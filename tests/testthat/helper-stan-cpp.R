# Helpers for testing the C++ versions of Stan functions (see
# pcd_stan_header()) against their Stan bodies.

# Compiles `model_code` (a model without a functions block) with the Stan
# bodies (`stan`), and with the C++ functions declared and the header
# included (`cpp`). With `fallback = TRUE` it also builds the C++ version
# with the nested autodiff dP/da instead of Eigen's (`fallback`). Returns a
# named list of CmdStanModel objects.
compile_cpp_and_stan <- function(model_code, name, fallback = FALSE) {
  testthat::skip_if_not_installed("cmdstanr")
  testthat::skip_if(
    is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))
  )
  build <- function(cpp, variant) {
    # nolint next: object_usage_linter.
    functions <- pcd_load_stan_functions(wrap_in_block = TRUE, cpp = cpp)
    model_dir <- file.path(tempdir(), paste0(name, variant))
    dir.create(model_dir, showWarnings = FALSE)
    path <- file.path(model_dir, paste0(name, ".stan"))
    code <- paste0(functions, "\n", model_code)
    # Rewriting an unchanged file would make cmdstanr recompile the model
    code_lines <- strsplit(code, "\n", fixed = TRUE)[[1]]
    if (!file.exists(path) || !identical(readLines(path), code_lines)) {
      writeLines(code, path)
    }
    model_args <- list(path, dir = model_dir)
    if (cpp) {
      model_args$user_header <- pcd_stan_header() # nolint: object_usage_linter.
    }
    if (variant == "fallback") {
      model_args$cpp_options <- list(
        "CXXFLAGS += -DPRIMARYCENSORED_STAN_DPDA"
      )
    }
    suppressMessages(suppressWarnings(
      do.call(cmdstanr::cmdstan_model, model_args)
    ))
  }
  models <- list(stan = build(FALSE, "stan"), cpp = build(TRUE, "cpp"))
  if (fallback) {
    models$fallback <- build(TRUE, "fallback")
  }
  models
}

# Code for a probe model. Its parameters are on the unconstrained scale, so
# there is no Jacobian term. `out_code` fills `vector[N] out` from
# `dist_id`, `params` and `pwindow`, and the data declared in `data_code`.
# The target is the weighted sum of the finite elements of `out`, and `out`
# is also a generated quantity.
cpp_probe_code <- function(data_code, out_code) {
  paste0(
    "data {\n  int dist_id;\n  real pwindow;\n  int N;\n", data_code,
    "\n  vector[N] w;\n}\n",
    "parameters {\n  vector[2] theta;\n}\n",
    "transformed parameters {\n  array[2] real params;\n",
    "  params[1] = dist_id == 1 ? theta[1] : exp(theta[1]);\n",
    "  params[2] = exp(theta[2]);\n}\n",
    "model {\n  vector[N] out;\n", out_code, "\n",
    "  for (i in 1:N) {\n",
    "    if (!is_inf(out[i]) && !is_nan(out[i])) {\n",
    "      target += w[i] * out[i];\n    }\n  }\n}\n",
    "generated quantities {\n  vector[N] out;\n", out_code, "\n}\n"
  )
}

# Unconstrained parameters for the probe models.
cpp_probe_theta <- function(dist_id, params) {
  if (dist_id == 1L) c(params[1], log(params[2])) else log(params)
}

# Expects the C++ and Stan values of `out` to match, including which are
# infinite or NaN.
expect_cpp_values_match <- function(models, data, params, info) {
  init <- list(theta = cpp_probe_theta(data$dist_id, params))
  stan <- stan_gq(models$stan, data, init, "out")
  cpp <- stan_gq(models$cpp, data, init, "out")
  testthat::expect_identical(is.nan(cpp), is.nan(stan), info = info)
  testthat::expect_identical(
    is.infinite(cpp) & !is.nan(cpp), is.infinite(stan) & !is.nan(stan),
    info = info
  )
  finite <- is.finite(stan)
  testthat::expect_equal(
    cpp[finite], stan[finite],
    tolerance = 1e-13, info = info
  )
}

# Expects the C++ log density and gradient to match Stan's autodiff. The
# weights are positive, so the target does not cancel: far in the lower
# tail each gradient is large, and cancelling weights would turn small
# relative errors into a large relative error of the sum. Below a log CDF
# of about -700 the terms are close to the smallest doubles and both
# gradients are rounding error. A log PMF also uses the log CDF one delay
# below, so elements below -500 get zero weight. With
# `upper_tail = TRUE` (for log PMFs), elements past the mode with a log PMF
# below -10 also get zero weight. There the log PMF is a difference of log
# CDFs close to 0, so small differences in how Stan and the C++ code
# differentiate the special functions are amplified by about 1 / PMF, as
# described in issue 404. The C++ gradient is checked against a 50-digit
# reference there instead.
expect_cpp_gradient_matches <- function(models, data, params, info,
                                        model = "cpp", tolerance = 1e-8,
                                        upper_tail = FALSE) {
  upars <- list(theta = cpp_probe_theta(data$dist_id, params))
  out <- stan_gq(models$stan, data, upars, "out")
  data$w[!is.na(out) & out < -500] <- 0
  if (upper_tail) {
    past_mode <- seq_along(out) > which.max(out)
    data$w[past_mode & out < -10] <- 0
  }
  stan <- stan_log_prob(models$stan, data, upars)
  cpp <- stan_log_prob(models[[model]], data, upars)
  testthat::expect_equal(cpp$lp, stan$lp, tolerance = 1e-12, info = info)
  testthat::expect_equal(
    cpp$grad, stan$grad,
    tolerance = tolerance, info = info
  )
}

# Skips unless CmdStan's log_prob method can be run directly. It needs
# CmdStan 2.31 or later, and on Windows the binary needs the TBB library on
# the path, which cmdstanr only sets up for its own runs.
skip_if_no_log_prob <- function() {
  testthat::skip_on_os("windows")
  testthat::skip_if(
    cmdstanr::cmdstan_version(error_on_NA = FALSE) < "2.31"
  )
}

# Log density and its gradient at parameters `upars` (a named list), from
# CmdStan's log_prob method at full precision. The probe models have only
# unconstrained parameters, so these are on the unconstrained scale.
stan_log_prob <- function(model, data, upars) {
  run_dir <- tempfile("logprob")
  dir.create(run_dir)
  data_file <- file.path(run_dir, "data.json")
  pars_file <- file.path(run_dir, "pars.json")
  out_file <- file.path(run_dir, "out.csv")
  cmdstanr::write_stan_json(data, data_file)
  cmdstanr::write_stan_json(upars, pars_file)
  run_log <- suppressWarnings(system2(
    model$exe_file(),
    c(
      "log_prob", paste0("constrained_params=", pars_file),
      "data", paste0("file=", data_file),
      "output", paste0("file=", out_file), "sig_figs=18"
    ),
    stdout = TRUE, stderr = TRUE
  ))
  if (!is.null(attr(run_log, "status"))) {
    stop(paste(run_log, collapse = "\n"), call. = FALSE)
  }
  out <- utils::read.csv(out_file, comment.char = "#")
  list(lp = out$lp__[1], grad = unlist(out[1, -1], use.names = FALSE))
}

# Generated quantities named `name` at parameters `init`, at full precision.
stan_gq <- function(model, data, init, name) {
  fit <- model$sample(
    data = data, init = list(init), fixed_param = TRUE, chains = 1,
    iter_sampling = 1, sig_figs = 18, refresh = 0, show_messages = FALSE,
    show_exceptions = FALSE
  )
  as.vector(fit$draws(name, format = "draws_matrix"))
}
