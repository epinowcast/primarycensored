skip_on_cran()

# Regression test for the Weibull gradient dropped in the upper tail.
# Reverse-mode `gamma_p(a, x)` in Stan Math returns zero gradients when
# x / a > 10 (https://github.com/stan-dev/math/issues/2006). The Weibull
# uniform primary solution needs the lower incomplete gamma function at
# x = (t / scale)^shape, so without a workaround its gradient is wrong for
# large delays. At Weibull(1.5, 5) the cutoff is t > 32.6.
#
# Gradients are only observable from a compiled model, so this builds a
# minimal one whose whole target is `primarycensored_lpmf` and runs
# `stan_gradient_at()` from helper-stan-gradient.R.

weibull_gradient_probe_model <- function() {
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
    "  int d;\n",
    "  real pwindow;\n",
    "}\n",
    "parameters {\n",
    "  real<lower=0> shape;\n",
    "  real<lower=0> scale;\n",
    "}\n",
    "model {\n",
    "  target += primarycensored_lpmf(\n",
    "    d | 3, {shape, scale}, pwindow, d + 1, 0, positive_infinity(),\n",
    "    1, {0.0}\n",
    "  );\n",
    "}\n"
  )
  path <- file.path(tempdir(), "pcd_weibull_tail_gradient.stan")
  writeLines(code, path)
  suppressMessages(suppressWarnings(cmdstanr::cmdstan_model(path)))
}

test_that("Weibull uniform primary lpmf gradients match finite differences
   in the upper tail", {
  model <- weibull_gradient_probe_model()

  # d = 10 and 20 are below the gamma_p cutoff and act as controls.
  for (d in c(10, 20, 33, 35, 40)) {
    for (pwindow in c(1, 2)) {
      label <- sprintf("d = %d, pwindow = %g", d, pwindow)
      res <- stan_gradient_at( # nolint: object_usage_linter.
        model,
        data = list(d = d, pwindow = pwindow),
        init = list(shape = 1.5, scale = 5)
      )

      expect_false(res$rejected, info = label)
      expect_length(res$gradient, 2)
      expect_equal(
        res$gradient, res$finite_diff,
        tolerance = 1e-4, info = label
      )
    }
  }
})
