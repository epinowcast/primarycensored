skip_on_cran()

# Gradient regression tests for #363. The log CDF of a generalised gamma in
# the lower tail is finite only if the incomplete gamma function is
# evaluated on the log scale, and its gradient with respect to `k` must stay
# accurate there. Gradients are only observable from a compiled model, so
# this builds a minimal one whose target is `primarycensored_lcdf` and runs
# `stan_gradient_at()` from helper-stan-gradient.R.

gengamma_probe_model <- function() {
  testthat::skip_if_not_installed("cmdstanr")
  testthat::skip_if(
    is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))
  )

  functions <- pcd_load_stan_functions(
    wrap_in_block = TRUE, write_to_file = FALSE
  )
  code <- paste0(
    functions, "\n",
    "data {\n",
    "  real d;\n",
    "  real pwindow;\n",
    "  real L;\n",
    "  real D;\n",
    "}\n",
    "parameters {\n",
    "  real<lower=0> shape;\n",
    "  real<lower=0> scale;\n",
    "  real<lower=0> k;\n",
    "}\n",
    "model {\n",
    "  target += primarycensored_lcdf(\n",
    "    d | 5, {shape, scale, k}, pwindow, L, D, 1, rep_array(0.0, 0)\n",
    "  );\n",
    "}\n"
  )
  path <- file.path(tempdir(), "pcd_gengamma_tail_gradient.stan")
  writeLines(code, path)
  suppressMessages(suppressWarnings(cmdstanr::cmdstan_model(path)))
}

gengamma_gradient_at <- function(model, d, params, pwindow = 1, L = 0,
                                 D = Inf) {
  stan_gradient_at( # nolint: object_usage_linter.
    model,
    data = list(d = d, pwindow = pwindow, L = L, D = D),
    init = list(shape = params[1], scale = params[2], k = params[3])
  )
}

expect_gradient_ok <- function(res, label) {
  testthat::expect_false(res$gradient_not_finite, info = label)
  testthat::expect_false(res$rejected, info = label)
  testthat::expect_length(res$gradient, 3)
  testthat::expect_true(all(is.finite(res$gradient)), info = label)
  # The analytic gradient must agree with the finite difference.
  testthat::expect_equal(
    res$gradient, res$finite_diff,
    tolerance = 1e-4, info = label
  )
}

test_that("primarycensored_lcdf has accurate finite gradients deep in the
   lower tail of a generalised gamma", {
  model <- gengamma_probe_model()
  # Values from -20 to -2400. The first two returned -inf or a wrong
  # gradient for k before the fix.
  cases <- list(
    list(d = 2, p = c(1, 5, 400)),
    list(d = 2, p = c(5, 5, 100)),
    list(d = 2, p = c(1.5, 3, 30)),
    list(d = 3, p = c(0.7, 4, 60)),
    list(d = 2, p = c(2, 6, 250)),
    list(d = 1, p = c(3, 10, 40))
  )
  for (case in cases) {
    for (pwindow in c(0.5, 1, 3)) {
      label <- sprintf(
        "d = %g, pwindow = %g, params = (%s)", case$d, pwindow,
        toString(case$p)
      )
      res <- gengamma_gradient_at(model, case$d, case$p, pwindow)
      expect_gradient_ok(res, label)
    }
  }
})

test_that("primarycensored_lcdf has finite gradients for extreme shapes", {
  model <- gengamma_probe_model()
  # k = 30000 puts x / (k + 1) between 0.5 and 0.9 with a log CDF far below
  # -600, which uses the wider series rule.
  res <- gengamma_gradient_at(
    model, 21000, c(1, 1, 30000),
    pwindow = 2000
  )
  expect_gradient_ok(res, "k = 30000")
  res <- gengamma_gradient_at(
    model, 0.6, c(12, 0.5, 1500),
    pwindow = 0.5
  )
  expect_gradient_ok(res, "shape = 12, k = 1500")
})

test_that("primarycensored_lcdf has finite gradients with truncation when
   both bounds are deep in the lower tail", {
  model <- gengamma_probe_model()
  cases <- list(
    list(d = 2, p = c(1, 5, 400), L = 0, D = 2.5),
    list(d = 2, p = c(5, 5, 100), L = 0, D = 2.5),
    list(d = 2, p = c(1, 5, 400), L = 1.5, D = 2.5)
  )
  for (case in cases) {
    label <- sprintf(
      "d = %g, L = %g, D = %g, params = (%s)", case$d, case$L, case$D,
      toString(case$p)
    )
    res <- gengamma_gradient_at(
      model, case$d, case$p,
      pwindow = 1, L = case$L, D = case$D
    )
    expect_gradient_ok(res, label)
  }
})

test_that("gradients are unchanged in the body of the distribution", {
  model <- gengamma_probe_model()
  for (d in c(0.5, 2, 6, 15)) {
    res <- gengamma_gradient_at(model, d, c(1.5, 3, 2), pwindow = 1)
    expect_gradient_ok(res, sprintf("d = %g", d))
  }
})
