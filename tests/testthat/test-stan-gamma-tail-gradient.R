skip_on_cran()

# Gradient regression tests for #381. Stan's `gamma_lcdf` has an inaccurate
# or failing gradient with respect to its shape in the lower tail, and for a
# shape of about 1000 or more anywhere (nan, or an exception after 100000
# iterations). Gradients are only observable from a compiled model, so these
# build minimal ones and run `stan_gradient_at()` from
# helper-stan-gradient.R.

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

# A function of (log x, a) on its own, with d standing in for x. The shape is
# unconstrained so that the gradient is with respect to the shape itself.
gamma_logx_probe_model <- function(fn = "gamma_lcdf_logx") {
  gamma_probe_model(
    paste0(fn, "(log(d), a)"),
    "  real a;\n",
    paste0("pcd_", fn, "_gradient")
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

# Derivative of log P(a, x) (or log Q(a, x) for `lower = FALSE`) with
# respect to a, by a fifth order central difference of `pgamma()`, which is
# accurate to about 1e-10 relative. The step in log a shrinks with the
# shape, because a fixed step of 1e-4 is too coarse for a of 1e5 or more.
ref_dlogp_da <- function(x, a, lower = TRUE) {
  h <- a * min(1e-4, 0.01 / sqrt(a))
  f <- function(e) pgamma(x, a + e, lower.tail = lower, log.p = TRUE)
  (-f(2 * h) + 8 * f(h) - 8 * f(-h) + f(-2 * h)) / (12 * h)
}

# The error is checked against the size of the derivative, with a small
# absolute allowance for where it is below 1e-10. CmdStan prints gradients
# to 6 significant figures, so the tolerance is relative 1e-5.
expect_gradient_close <- function(gradient, expected, label) {
  testthat::expect_length(gradient, 1)
  if (length(gradient) == 1) {
    testthat::expect_lt(
      abs(gradient - expected), 1e-5 * abs(expected) + 1e-10,
      label = label
    )
  }
}

test_that("primarycensored_lcdf has accurate finite gradients deep in the
   lower tail of a Gamma delay", {
  model <- gamma_delay_probe_model()
  # Log CDFs from -20 to -2400. The first returned -inf before the fix.
  cases <- list(
    list(d = 2, p = c(400, 1 / 5)),
    list(d = 2, p = c(100, 1 / 5)),
    list(d = 1.5, p = c(60, 1 / 4)),
    list(d = 3, p = c(250, 1 / 6)),
    list(d = 1, p = c(12, 1 / 4)),
    list(d = 1, p = c(40, 1 / 10))
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

test_that("primarycensored_lcdf has finite gradients for a Gamma delay with
   a large shape", {
  model <- gamma_delay_probe_model()
  # Shapes of 1000 or more from the body into the upper tail, where
  # `gamma_lcdf` returned nan gradients or threw after 100000 iterations.
  cases <- list(
    list(d = 1400, p = c(1500, 1), pwindow = 10),
    list(d = 1500, p = c(1500, 1), pwindow = 10),
    list(d = 1530, p = c(1500, 1), pwindow = 10),
    list(d = 1700, p = c(1500, 1), pwindow = 10),
    list(d = 3000, p = c(1500, 1), pwindow = 10),
    list(d = 1050, p = c(1000, 1), pwindow = 1),
    list(d = 2900, p = c(3000, 1), pwindow = 10),
    list(d = 3150, p = c(3000, 1), pwindow = 10),
    list(d = 10400, p = c(10000, 1), pwindow = 10)
  )
  for (case in cases) {
    label <- sprintf(
      "d = %g, pwindow = %g, params = (%s)", case$d, case$pwindow,
      toString(case$p)
    )
    res <- gamma_gradient_at(model, case$d, case$p, case$pwindow)
    expect_gamma_gradient_ok(res, label)
  }
})

test_that("primarycensored_lcdf has finite gradients deep in the upper tail
   of a Gamma delay", {
  model <- gamma_delay_probe_model()
  # The CDF is 1 to double precision here, and Stan's gradient was nan. The
  # log CDF does not depend on the parameters, so the gradient with respect
  # to each is the 1 from the log Jacobian of its lower bound.
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
    # The log CDFs are in the thousands, so CmdStan's finite differences
    # carry rounding noise of 1e-4. Compare with a reference in R instead.
    expect_equal(
      res$gradient,
      ref_gamma_delay_gradient(case$d, case$p, 1, case$L, case$D),
      tolerance = 1e-5, info = label
    )
  }
})

test_that("gradients are unchanged in the body of a Gamma delay", {
  model <- gamma_delay_probe_model()
  for (d in c(0.5, 2, 6, 15)) {
    res <- gamma_gradient_at(model, d, c(2.3, 0.5), pwindow = 1)
    expect_gamma_gradient_ok(res, sprintf("d = %g", d))
  }
})

test_that("gamma_lcdf_logx gradient with respect to the shape is accurate", {
  model <- gamma_logx_probe_model()
  # Non-integer and integer shapes, from the lower tail to the upper tail
  for (a in c(0.7, 2, 5.5, 9.5, 10, 10.5, 20.5, 100, 700.5, 1000, 1500.5,
              3000, 10000)) {
    for (frac in c(0.3, 0.8, 0.95, 1, 1.02, 1.1, 1.5, 3)) {
      x <- frac * a
      res <- stan_gradient_at(
        model,
        data = list(d = x, pwindow = 1, L = 0, D = Inf),
        init = list(a = a)
      )
      label <- sprintf("a = %g, x over a = %g", a, frac)
      expect_false(res$rejected, info = label)
      expect_false(res$gradient_not_finite, info = label)
      expect_gradient_close(res$gradient, ref_dlogp_da(x, a), label)
    }
  }
})

test_that("gamma_lcdf_logx gradient is accurate for integer shapes", {
  model <- gamma_logx_probe_model()
  # The continued fraction has a zero numerator at step i = a
  for (a in c(10, 12, 25, 50, 100, 400)) {
    for (frac in c(1, 1.02, 1.1, 1.5, 3)) {
      x <- frac * (a + 1)
      res <- stan_gradient_at(
        model,
        data = list(d = x, pwindow = 1, L = 0, D = Inf),
        init = list(a = a)
      )
      expect_gradient_close(
        res$gradient, ref_dlogp_da(x, a),
        sprintf("a = %g, x over (a + 1) = %g", a, frac)
      )
    }
  }
})

test_that("gamma_lccdf_cf_logx gradient is accurate for integer shapes", {
  model <- gamma_logx_probe_model("gamma_lccdf_cf_logx")
  # The numerator of step i is i (i - a), which is zero at i = a. The
  # derivative with respect to a is only accurate there if the fraction has
  # converged by step a, which it has for a of 10 or more
  for (a in c(10, 11, 12, 15, 20, 25, 100)) {
    for (frac in c(1, 1.5, 3)) {
      x <- frac * (a + 1)
      res <- stan_gradient_at(
        model,
        data = list(d = x, pwindow = 1, L = 0, D = Inf),
        init = list(a = a)
      )
      expect_gradient_close(
        res$gradient, ref_dlogp_da(x, a, lower = FALSE),
        sprintf("a = %g, x over (a + 1) = %g", a, frac)
      )
    }
  }
})

test_that("gamma_lcdf_logx_pair gradients with respect to the shape are
   accurate", {
  for (component in 1:2) {
    model <- gamma_probe_model(
      paste0("gamma_lcdf_logx_pair(log(d), a)[", component, "]"),
      "  real a;\n",
      paste0("pcd_gamma_pair_gradient_", component)
    )
    # The second component is the CDF of a + 1
    shift <- component - 1
    for (a in c(0.3, 2, 5.5, 9.5, 10.5, 20.5, 100, 700.5, 3000, 10000)) {
      for (frac in c(1e-4, 0.05, 0.3, 0.6, 0.95, 1, 1.05, 1.5, 3)) {
        x <- frac * (a + 1)
        res <- stan_gradient_at(
          model,
          data = list(d = x, pwindow = 1, L = 0, D = Inf),
          init = list(a = a)
        )
        label <- sprintf(
          "component = %d, a = %g, x over (a + 1) = %g", component, a, frac
        )
        expect_false(res$rejected, info = label)
        expect_false(res$gradient_not_finite, info = label)
        expect_gradient_close(
          res$gradient,
          ref_dlogp_da(x, a + shift),
          label
        )
      }
    }
  }
})
