skip_on_cran()

# Regression test for the NaN gradient reported in #333. `lognormal_lcdf`
# underflows to -inf deep in the lower tail of a narrow lognormal and its
# autodiff partial is then 0 / 0, which Stan's reverse pass chains into `mu`
# and `sigma` even where the term carries zero weight. The symptom is a
# finite log density with a non-finite gradient.
#
# Gradients are only observable from a compiled model, so this builds a
# minimal one whose whole target is `primarycensored_lcdf` and drives
# `diagnose test=gradient` directly rather than through `$diagnose()`, whose
# exit-status handling varies across cmdstanr versions.

gradient_probe_model <- function() {
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)))

  functions <- pcd_load_stan_functions(
    wrap_in_block = TRUE, write_to_file = FALSE
  )
  code <- paste0(
    functions, "\n",
    "data {\n",
    "  real d;\n",
    "  real pwindow;\n",
    "  int primary_id;\n",
    "  int n_primary;\n",
    "  array[n_primary] real primary_params;\n",
    "}\n",
    "parameters {\n",
    "  real mu;\n",
    "  real<lower=0> sigma;\n",
    "}\n",
    "model {\n",
    "  target += primarycensored_lcdf(\n",
    "    d | 1, {mu, sigma}, pwindow, 0, positive_infinity(),\n",
    "    primary_id, primary_params\n",
    "  );\n",
    "}\n"
  )
  path <- file.path(tempdir(), "pcd_lognormal_tail_gradient.stan")
  writeLines(code, path)
  suppressMessages(suppressWarnings(cmdstanr::cmdstan_model(path)))
}

# Returns the log density and the analytic and finite-difference gradients.
gradient_at <- function(model, d, sigma, mu = 1.8, pwindow = 1,
                        primary_id = 1, primary_params = numeric(0)) {
  data_file <- tempfile(fileext = ".json")
  cmdstanr::write_stan_json(
    list(
      d = d, pwindow = pwindow, primary_id = primary_id,
      n_primary = length(primary_params),
      primary_params = as.array(primary_params)
    ),
    data_file
  )
  init_file <- tempfile(fileext = ".json")
  cmdstanr::write_stan_json(list(mu = mu, sigma = sigma), init_file)

  out <- suppressWarnings(system2(
    model$exe_file(),
    c(
      "diagnose", "test=gradient",
      paste0("data file=", data_file),
      paste0("init=", init_file),
      "output", paste0("file=", tempfile(fileext = ".csv"))
    ),
    stdout = TRUE, stderr = TRUE
  ))

  rejected <- any(grepl("Rejecting initial value", out))
  not_finite <- any(grepl("Gradient evaluated at the initial value", out))
  rows <- grep("^\\s+\\d+\\s+", out, value = TRUE)
  parsed <- lapply(strsplit(trimws(rows), "\\s+"), as.numeric)

  list(
    rejected = rejected,
    gradient_not_finite = not_finite,
    gradient = vapply(parsed, function(x) x[3], numeric(1)),
    finite_diff = vapply(parsed, function(x) x[4], numeric(1))
  )
}

test_that("primarycensored_lcdf has finite gradients in the lower tail of a
   narrow lognormal", {
  model <- gradient_probe_model()

  # Each of these returned a finite log density with a NaN gradient before
  # the guard. `d = 1` with a uniform primary is included as a control: it
  # always worked, because `q = 0` took the existing `log(0)` guard.
  cases <- list(
    list(d = 1.0, sigma = 0.05, primary_id = 1, pp = numeric(0)),
    list(d = 1.5, sigma = 0.05, primary_id = 1, pp = numeric(0)),
    list(d = 2.0, sigma = 0.05, primary_id = 1, pp = numeric(0)),
    list(d = 2.0, sigma = 0.03, primary_id = 1, pp = numeric(0)),
    list(d = 3.0, sigma = 0.03, primary_id = 1, pp = numeric(0)),
    list(d = 1.0, sigma = 0.05, primary_id = 2, pp = 0.14),
    list(d = 1.0, sigma = 0.10, primary_id = 2, pp = 0.14)
  )

  for (case in cases) {
    label <- sprintf(
      "d = %g, sigma = %g, primary_id = %d", case$d, case$sigma,
      case$primary_id
    )
    res <- gradient_at(
      model, case$d, case$sigma,
      primary_id = case$primary_id, primary_params = case$pp
    )

    expect_false(res$gradient_not_finite, info = label)
    expect_false(res$rejected, info = label)
    expect_true(all(is.finite(res$gradient)), info = label)
    # The analytic gradient must agree with the finite difference.
    expect_equal(
      res$gradient, res$finite_diff,
      tolerance = 1e-4, info = label
    )
  }
})

test_that("a genuinely zero density is reported as log(0) rather than a
   non-finite gradient", {
  model <- gradient_probe_model()

  # Further into the tail the density really is zero. That must surface as
  # log(0), which points at the cause, not as a non-finite gradient.
  res <- gradient_at(model, d = 1.0, sigma = 0.03)

  expect_true(res$rejected)
  expect_false(res$gradient_not_finite)
})

test_that("the underflow guard leaves values in the normal range unchanged", {
  # The guard must only ever drop terms that are already unrepresentable, so
  # the analytical solution still has to match numeric integration wherever
  # it did before.
  for (sigma in c(0.05, 0.1, 0.5, 1)) {
    for (pwindow in c(1, 2)) {
      obj <- new_pcens(plnorm, dunif, list(), meanlog = 1.8, sdlog = sigma)
      q_values <- seq(0.5, 12, by = 0.5)

      analytic <- pcens_cdf(obj, q = q_values, pwindow = pwindow)
      numeric <- pcens_cdf(
        obj,
        q = q_values, pwindow = pwindow, use_numeric = TRUE
      )

      expect_equal(
        analytic, numeric,
        tolerance = 1e-6,
        info = sprintf("sigma = %g, pwindow = %g", sigma, pwindow)
      )
    }
  }
})
