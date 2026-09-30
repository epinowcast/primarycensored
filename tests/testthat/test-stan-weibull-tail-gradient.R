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
# `stan_gradient_at()` from helper-stan-gradient.R. The parameters are on the
# log scale so there is no Jacobian term.
#
# CmdStan's own finite differences are not used as the reference. In the
# tail the analytical PMF is a difference of CDFs close to 1, so finite
# differences of it are only good to a few per cent. The reference instead
# integrates the difference of Weibull survival functions over the primary
# window, which does not cancel.

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
    "  real log_shape;\n",
    "  real log_scale;\n",
    "}\n",
    "model {\n",
    "  target += primarycensored_lpmf(\n",
    "    d | 3, {exp(log_shape), exp(log_scale)}, pwindow, d + 1, 0,\n",
    "    positive_infinity(), 1, {0.0}\n",
    "  );\n",
    "}\n"
  )
  path <- file.path(tempdir(), "pcd_weibull_tail_gradient.stan")
  writeLines(code, path)
  suppressMessages(suppressWarnings(cmdstanr::cmdstan_model(path)))
}

# Log PMF of a Weibull delay with a uniform primary event, computed from
# survival functions so that it stays accurate in the upper tail.
weibull_uniform_lpmf <- function(log_params, d, pwindow) {
  shape <- exp(log_params[1])
  scale <- exp(log_params[2])
  surv <- function(t) exp(-(pmax(t, 0) / scale)^shape)
  integrand <- function(u) surv(d - u) - surv(d + 1 - u)
  prob <- stats::integrate(
    integrand, 0, pwindow,
    rel.tol = 1e-12, abs.tol = 0
  )$value
  log(prob / pwindow)
}

# Central difference gradient of `weibull_uniform_lpmf()`.
weibull_uniform_lpmf_grad <- function(log_params, d, pwindow, h = 1e-5) {
  vapply(seq_along(log_params), function(i) {
    shift <- replace(numeric(length(log_params)), i, h)
    (weibull_uniform_lpmf(log_params + shift, d, pwindow) -
      weibull_uniform_lpmf(log_params - shift, d, pwindow)) / (2 * h)
  }, numeric(1))
}

test_that("Weibull uniform primary lpmf gradients are correct in the upper
   tail", {
  model <- weibull_gradient_probe_model()

  # Delays below the gamma_p cutoff act as controls. The cutoff is at
  # t > 32.6 for Weibull(1.5, 5) and t > 4.7 for Weibull(3, 2).
  cases <- list(
    list(shape = 1.5, scale = 5, d = c(10, 20, 33, 35, 40)),
    list(shape = 3, scale = 2, d = c(3, 5, 6))
  )

  for (case in cases) {
    log_params <- log(c(case$shape, case$scale))
    for (d in case$d) {
      for (pwindow in c(1, 2)) {
        label <- sprintf(
          "Weibull(%g, %g), d = %d, pwindow = %g",
          case$shape, case$scale, d, pwindow
        )
        res <- stan_gradient_at( # nolint: object_usage_linter.
          model,
          data = list(d = d, pwindow = pwindow),
          init = list(log_shape = log_params[1], log_scale = log_params[2])
        )

        expect_false(res$rejected, info = label)
        expect_length(res$gradient, 2)
        expect_equal(
          res$gradient,
          weibull_uniform_lpmf_grad(log_params, d, pwindow),
          tolerance = 1e-4, info = label
        )
      }
    }
  }
})
