skip_on_cran()

# Gradients are only observable from a compiled model, so these build
# minimal ones and use `stan_gradient_at()` from helper-stan-gradient.R.

gamma_probe_model <- function(target, parameters, name) {
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
    "parameters {\n", parameters, "}\n",
    "model {\n  target += ", target, ";\n}\n"
  )
  path <- file.path(tempdir(), paste0(name, ".stan"))
  writeLines(code, path)
  suppressMessages(suppressWarnings(cmdstanr::cmdstan_model(path)))
}

# Gamma delay with a uniform primary through the analytical solution
gamma_delay_probe_model <- function() {
  gamma_probe_model(
    paste0(
      "primarycensored_lcdf(\n",
      "    d | 2, {shape, rate}, pwindow, L, D, 1, rep_array(0.0, 0)\n",
      "  )"
    ),
    "  real<lower=0> shape;\n  real<lower=0> rate;\n",
    "pcd_gamma_tail_gradient"
  )
}

# An expression in (d, a) on its own, with d standing in for log x. The shape
# is unconstrained so that the gradient is with respect to the shape itself.
gamma_logx_probe_model <- function(target, name) {
  gamma_probe_model(target, "  real a;\n", name)
}

# Shape gradient of `target` at log x = d
gamma_logx_gradient_at <- function(model, d, a) {
  stan_gradient_at( # nolint: object_usage_linter.
    model,
    data = list(d = d, pwindow = 1, L = 0, D = Inf),
    init = list(a = a)
  )
}

gamma_gradient_at <- function(model, d, params, pwindow = 1, L = 0,
                              D = Inf) {
  stan_gradient_at( # nolint: object_usage_linter.
    model,
    data = list(d = d, pwindow = pwindow, L = L, D = D),
    init = list(shape = params[1], rate = params[2])
  )
}

expect_gamma_gradient_ok <- function(res, label) {
  testthat::expect_false(res$gradient_not_finite, info = label)
  testthat::expect_false(res$rejected, info = label)
  testthat::expect_length(res$gradient, 2)
  testthat::expect_true(all(is.finite(res$gradient)), info = label)
  # The analytic gradient must agree with the finite difference.
  testthat::expect_equal(
    res$gradient, res$finite_diff,
    tolerance = 1e-4, info = label
  )
}

# Derivative of log P(a, x) with respect to a, by a fifth order central
# difference of `pgamma()` with a step that shrinks with a
ref_dlogp_da <- function(x, a) {
  h <- a * min(1e-4, 0.01 / sqrt(a))
  f <- function(e) pgamma(x, a + e, log.p = TRUE)
  (-f(2 * h) + 8 * f(h) - 8 * f(-h) + f(-2 * h)) / (12 * h)
}

# Relative 1e-5, as CmdStan prints gradients to 6 significant figures
expect_gradient_close <- function(gradient, expected, label) {
  testthat::expect_length(gradient, 1)
  if (length(gradient) == 1) {
    testthat::expect_lt(
      abs(gradient - expected), 1e-5 * abs(expected) + 1e-10,
      label = label
    )
  }
}

test_that("primarycensored_lcdf has accurate finite gradients for a Gamma
   delay", {
  model <- gamma_delay_probe_model()
  # From the lower tail (log CDFs down to -2400) to the bulk
  cases <- list(
    list(d = 2, p = c(400, 1 / 5)),
    list(d = 2, p = c(100, 1 / 5)),
    list(d = 1.5, p = c(60, 1 / 4)),
    list(d = 3, p = c(250, 1 / 6)),
    list(d = 1, p = c(12, 1 / 4)),
    list(d = 1, p = c(40, 1 / 10)),
    list(d = 2, p = c(20, 1)),
    list(d = 10, p = c(100, 1)),
    list(d = 6, p = c(2.3, 0.5))
  )
  for (case in cases) {
    for (pwindow in c(0.5, 1, 3)) {
      label <- sprintf(
        "d = %g, pwindow = %g, params = (%s)", case$d, pwindow,
        toString(case$p)
      )
      res <- gamma_gradient_at(model, case$d, case$p, pwindow)
      expect_gamma_gradient_ok(res, label)
    }
  }
})

test_that("primarycensored_lcdf has finite gradients deep in the upper tail
   of a Gamma delay", {
  model <- gamma_delay_probe_model()
  # The CDF is 1 to double precision, so the gradient with respect to each
  # parameter is the 1 from the log Jacobian of its lower bound
  cases <- list(
    list(d = 60, p = c(2, 1), pwindow = 1),
    list(d = 120, p = c(20, 1), pwindow = 1),
    list(d = 1200, p = c(30, 0.1), pwindow = 1),
    list(d = 500, p = c(100, 1), pwindow = 1)
  )
  for (case in cases) {
    label <- sprintf(
      "d = %g, pwindow = %g, params = (%s)", case$d, case$pwindow,
      toString(case$p)
    )
    res <- gamma_gradient_at(model, case$d, case$p, case$pwindow)
    expect_false(res$rejected, info = label)
    expect_false(res$gradient_not_finite, info = label)
    expect_true(all(is.finite(res$gradient)), info = label)
    expect_equal(res$gradient, c(1, 1), tolerance = 1e-6, info = label)
  }
})

test_that("primarycensored_lcdf has finite gradients with truncation when
   both bounds are deep in the lower tail of a Gamma delay", {
  model <- gamma_delay_probe_model()
  cases <- list(
    list(d = 2, p = c(400, 1 / 5), L = 0, D = 2.5),
    list(d = 2, p = c(100, 1 / 5), L = 0, D = 2.5),
    list(d = 2, p = c(400, 1 / 5), L = 1.5, D = 2.5)
  )
  for (case in cases) {
    label <- sprintf(
      "d = %g, L = %g, D = %g, params = (%s)", case$d, case$L, case$D,
      toString(case$p)
    )
    res <- gamma_gradient_at(
      model, case$d, case$p,
      pwindow = 1, L = case$L, D = case$D
    )
    expect_false(res$gradient_not_finite, info = label)
    expect_false(res$rejected, info = label)
    # CmdStan's finite differences are noisy for log CDFs in the thousands
    expect_equal(
      res$gradient,
      ref_gamma_delay_gradient(case$d, case$p, 1, case$L, case$D),
      tolerance = 1e-5, info = label
    )
  }
})

test_that("gamma_lcdf_logx shape gradient is accurate in the lower tail and
   bulk", {
  model <- gamma_logx_probe_model(
    "gamma_lcdf_logx(d, a)", "pcd_gamma_lcdf_logx_gradient"
  )
  for (a in c(0.7, 1, 2, 5, 5.5, 9.5, 20, 100)) {
    for (frac in c(1e-4, 0.05, 0.2, 0.45, 0.55, 0.8, 1, 1.5, 3)) {
      x <- frac * (a + 1)
      res <- gamma_logx_gradient_at(model, log(x), a)
      label <- sprintf("a = %g, x over (a + 1) = %g", a, frac)
      expect_false(res$rejected, info = label)
      expect_false(res$gradient_not_finite, info = label)
      expect_gradient_close(res$gradient, ref_dlogp_da(x, a), label)
    }
  }
  # log P = a log x - lgamma(a + 1) as x tends to 0
  for (a in c(0.01, 1, 10, 400)) {
    res <- gamma_logx_gradient_at(model, -1000, a)
    expect_gradient_close(
      res$gradient, -1000 - digamma(a + 1), sprintf("a = %g", a)
    )
  }
})

test_that("summed shape gradient of the Gamma uniform primary PMF matches the
   reference", {
  # Issue #395: shape 20, rate 2, delays 0 to 19. The reference is a central
  # difference of the sum of the log PMFs integrated with `pgamma()`.
  model <- gamma_probe_model(
    paste0(
      "sum(primarycensored_sone_lpmf_vectorized(\n",
      "    19, 0, positive_infinity(), 2, {a, 2.0}, pwindow, 1,\n",
      "    rep_array(0.0, 0)\n",
      "  ))"
    ),
    "  real a;\n", "pcd_gamma_pmf_shape_gradient"
  )
  for (case in list(c(1, -5.194901), c(3, -6.466540))) {
    res <- stan_gradient_at( # nolint: object_usage_linter.
      model,
      data = list(d = 0, pwindow = case[1], L = 0, D = Inf),
      init = list(a = 20)
    )
    expect_gradient_close(
      res$gradient, case[2], sprintf("pwindow = %g", case[1])
    )
  }
})
