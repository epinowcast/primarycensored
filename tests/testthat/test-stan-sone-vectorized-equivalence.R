skip_on_cran()

# primarycensored_sone_lpmf_vectorized() shares the analytical uniform
# primary terms between integer nodes. These tests check it against a frozen
# copy of the previous implementation (stan/reference-sone-lpmf-vectorized.stan)
# for values and gradients.
#
# Gradients need a compiled model, so this builds one whose target is a
# weighted sum of the PMF and the finite log PMF, and drives CmdStan's
# `log_prob` method directly. Values come from `generated quantities` under
# `fixed_param`. Both are printed at full precision.

equivalence_model <- function() {
  testthat::skip_if_not_installed("cmdstanr")
  testthat::skip_if(
    is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))
  )
  local_stan <- file.path("inst", "stan", "functions")
  stan_path <- if (dir.exists(local_stan)) local_stan else pcd_stan_path()
  functions <- pcd_load_stan_functions(stan_path = stan_path)
  reference <- readLines(
    testthat::test_path("stan", "reference-sone-lpmf-vectorized.stan")
  )
  code <- paste(
    c(
      "functions {", functions, reference, "}",
      "data {",
      "  int max_delay;",
      "  real L;",
      "  real D;",
      "  int dist_id;",
      "  int n_params;",
      "  real pwindow;",
      "  int primary_id;",
      "  int n_primary;",
      "  array[n_primary] real primary_params;",
      "  vector[max_delay + 1] w;",
      "  int use_ref;",
      "  int target_type;",
      "}",
      "parameters {",
      "  array[n_params] real params;",
      "}",
      "model {",
      "  if (target_type > 0) {",
      "    vector[max_delay + 1] lpmf = use_ref",
      "      ? ref_primarycensored_sone_lpmf_vectorized(",
      "          max_delay, L, D, dist_id, params, pwindow, primary_id,",
      "          primary_params)",
      "      : primarycensored_sone_lpmf_vectorized(",
      "          max_delay, L, D, dist_id, params, pwindow, primary_id,",
      "          primary_params);",
      "    for (i in 1:(max_delay + 1)) {",
      "      if (w[i] != 0) {",
      "        target += w[i] * (target_type == 1 ? exp(lpmf[i]) : lpmf[i]);",
      "      }",
      "    }",
      "  }",
      "}",
      "generated quantities {",
      "  vector[max_delay + 1] new_lpmf =",
      "    primarycensored_sone_lpmf_vectorized(",
      "      max_delay, L, D, dist_id, params, pwindow, primary_id,",
      "      primary_params);",
      "  vector[max_delay + 1] ref_lpmf =",
      "    ref_primarycensored_sone_lpmf_vectorized(",
      "      max_delay, L, D, dist_id, params, pwindow, primary_id,",
      "      primary_params);",
      "  vector[2] bounds = primarycensored_truncation_bounds(",
      "    L, D, dist_id, params, pwindow, primary_id, primary_params);",
      "  real log_mass = primarycensored_log_normalizer(",
      "    bounds[2], bounds[1], L);",
      "}"
    ),
    collapse = "\n"
  )
  path <- file.path(tempdir(), "pcd_sone_vectorized_equivalence.stan")
  writeLines(code, path)
  suppressMessages(suppressWarnings(cmdstanr::cmdstan_model(path)))
}

run_cmdstan <- function(model, args) {
  out <- file.path(tempdir(), "pcd_equivalence_out.csv")
  cmdstan_log <- suppressWarnings(system2(
    model$exe_file(),
    c(args, "output", paste0("file=", out), "sig_figs=18"),
    stdout = TRUE, stderr = TRUE
  ))
  if (!is.null(attr(cmdstan_log, "status"))) {
    stop(
      "CmdStan failed:\n", paste(cmdstan_log, collapse = "\n"),
      call. = FALSE
    )
  }
  csv <- readLines(out)
  utils::read.csv(text = csv[!startsWith(csv, "#")])
}

# Returns the new and reference log PMFs, and the gradients of a weighted
# sum of the PMF and of the log PMF, at one parameter set. The log PMF sum
# leaves out bins whose untruncated mass is below 1e-6: the gradient of a log
# PMF is 1 / mass times a difference of CDF gradients, so there rounding in
# either implementation is amplified past any fixed tolerance.
compare_at <- function(model, case, params) {
  n <- case$max_delay + 1
  stan_data <- c(
    case[c("max_delay", "L", "D", "dist_id", "pwindow", "primary_id")],
    list(
      n_params = length(params),
      n_primary = length(case$primary_params),
      primary_params = as.array(case$primary_params)
    )
  )
  init_file <- tempfile(fileext = ".json")
  cmdstanr::write_stan_json(list(params = as.array(params)), init_file)
  run <- function(args, w, use_ref = 0, target_type = 0) {
    data_file <- tempfile(fileext = ".json")
    cmdstanr::write_stan_json(
      c(stan_data, list(
        w = as.array(w), use_ref = use_ref,
        target_type = target_type
      )),
      data_file
    )
    run_cmdstan(model, c(
      args, "data", paste0("file=", data_file), paste0("init=", init_file)
    ))
  }
  grad <- function(w, use_ref, target_type) {
    unname(unlist(run(
      c(
        "method=log_prob", "jacobian=0",
        paste0("constrained_params=", init_file)
      ),
      w, use_ref, target_type
    )))
  }
  w <- seq(0.5, 1.5, length.out = n)
  gq <- run(c(
    "method=sample", "algorithm=fixed_param", "num_samples=1",
    "num_warmup=0"
  ), w)
  new_lpmf <- unname(unlist(gq[startsWith(names(gq), "new_lpmf")]))
  ref_lpmf <- unname(unlist(gq[startsWith(names(gq), "ref_lpmf")]))
  mass <- exp(ref_lpmf + gq$log_mass)
  w_log <- ifelse(!is.na(mass) & mass > 1e-6, w, 0)
  list(
    new = new_lpmf, ref = ref_lpmf,
    new_pmf_grad = grad(w, 0, 1), ref_pmf_grad = grad(w, 1, 1),
    new_lpmf_grad = grad(w_log, 0, 2), ref_lpmf_grad = grad(w_log, 1, 2)
  )
}

# Largest difference relative to max(1, |reference|). Matching non-finite
# entries (e.g. -Inf log PMFs) count as equal.
max_rel_diff <- function(new, ref) {
  both <- is.finite(new) & is.finite(ref)
  testthat::expect_identical(is.finite(new), is.finite(ref))
  testthat::expect_identical(new[!both], ref[!both])
  if (!any(both)) {
    return(0)
  }
  max(abs(new[both] - ref[both]) / pmax(1, abs(ref[both])))
}

cases <- function() {
  lognormal <- list(dist_id = 1, params = list(
    c(1.5, 0.5), c(0.2, 1.2), c(-1, 0.3), c(3.5, 0.05), c(4, 2)
  ))
  gamma_dist <- list(dist_id = 2, params = list(
    c(2, 0.5), c(0.5, 0.2), c(20, 4), c(1.2, 0.02), c(50, 0.5)
  ))
  weibull <- list(dist_id = 3, params = list(
    c(1.5, 5), c(0.7, 2), c(3, 30), c(5, 0.8)
  ))
  gengamma <- list(dist_id = 5, params = list(
    c(1.5, 3, 2), c(0.8, 1, 0.5)
  ))
  dists <- list(lognormal, gamma_dist, weibull, gengamma)
  settings <- list(
    list(max_delay = 20, L = 0, D = 21, pwindow = 1),
    list(max_delay = 40, L = 0, D = Inf, pwindow = 1),
    list(max_delay = 15, L = 0, D = 30, pwindow = 1),
    list(max_delay = 0, L = 0, D = 1, pwindow = 1),
    list(max_delay = 1, L = 0, D = Inf, pwindow = 1),
    list(max_delay = 20, L = 3, D = 21, pwindow = 1),
    list(max_delay = 20, L = 2.5, D = 25.5, pwindow = 1),
    list(max_delay = 20, L = -Inf, D = 21, pwindow = 1),
    list(max_delay = 20, L = 0, D = 21, pwindow = 2),
    list(max_delay = 20, L = 4, D = Inf, pwindow = 3),
    list(max_delay = 5, L = 0, D = 6, pwindow = 7),
    list(max_delay = 20, L = 0, D = 21, pwindow = 1.5)
  )
  out <- list()
  for (dist in dists) {
    for (setting in settings) {
      out[[length(out) + 1]] <- c(
        setting,
        list(
          dist_id = dist$dist_id, params = dist$params,
          primary_id = 1, primary_params = numeric(0)
        )
      )
    }
  }
  # Paths that are not optimised must be unchanged as well
  out[[length(out) + 1]] <- list(
    max_delay = 15, L = 0, D = 16, pwindow = 1, dist_id = 1,
    params = list(c(1.5, 0.5)), primary_id = 2, primary_params = 0.2
  )
  out[[length(out) + 1]] <- list(
    max_delay = 15, L = 0, D = 16, pwindow = 1, dist_id = 4,
    params = list(0.3), primary_id = 1, primary_params = numeric(0)
  )
  out
}

# Values are identical. Gradients differ only in the order adjoints are
# summed. Differences of nearly equal terms amplify that rounding, most in
# the log PMF, hence the looser gradient tolerances.
test_that(
  "primarycensored_sone_lpmf_vectorized matches the previous implementation
  in values and gradients",
  {
    model <- equivalence_model()
    max_diff <- c(value = 0, pmf_grad = 0, lpmf_grad = 0)
    for (case in cases()) {
      for (params in case$params) {
        res <- compare_at(model, case, params)
        info <- paste(
          "dist", case$dist_id, "params", toString(params),
          "max_delay", case$max_delay, "L", case$L, "D", case$D,
          "pwindow", case$pwindow, "primary", case$primary_id
        )
        expect_length(res$new, case$max_delay + 1)
        expect_length(res$new_pmf_grad, length(params) + 1)
        expect_length(res$new_lpmf_grad, length(params) + 1)
        diff <- c(
          value = max_rel_diff(res$new, res$ref),
          pmf_grad = max_rel_diff(res$new_pmf_grad, res$ref_pmf_grad),
          lpmf_grad = max_rel_diff(res$new_lpmf_grad, res$ref_lpmf_grad)
        )
        expect_lte(diff[["value"]], 1e-12, label = info)
        expect_lte(diff[["pmf_grad"]], 1e-11, label = info)
        expect_lte(diff[["lpmf_grad"]], 1e-10, label = info)
        max_diff <- pmax(max_diff, diff)
      }
    }
    message(
      "Maximum relative differences: ",
      toString(paste(names(max_diff), signif(max_diff, 3)))
    )
  }
)
