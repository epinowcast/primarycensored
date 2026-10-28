skip_on_cran()

# The reference is `pgamma(log.p = TRUE)`, with relative tolerance 1e-9
# unless stated.

# Relative error of a log CDF, with an absolute allowance of 1e-15 for where
# the CDF is within 1e-8 of 1
expect_lcdf_close <- function(actual, expected, label, tolerance = 1e-9) {
  testthat::expect_lte(
    abs(actual - expected), tolerance * abs(expected) + 1e-15,
    label = label
  )
}

# Shapes and rates that put the log CDF below -10 at the delays used
gamma_tail_cases <- data.frame(
  shape = c(400, 100, 60, 250, 40, 12),
  rate = c(1 / 5, 1 / 5, 1 / 4, 1 / 6, 1 / 10, 1 / 4)
)

test_that("dist_lcdf for a Gamma delay is finite and accurate deep in the
   lower tail", {
  for (i in seq_len(nrow(gamma_tail_cases))) {
    with(gamma_tail_cases[i, ], {
      for (y in c(0.05, 0.3, 1, 2, 3, 4.5)) {
        expect_equal(
          dist_lcdf(y, c(shape, rate), 2),
          pgamma(y * rate, shape, log.p = TRUE),
          tolerance = 1e-9,
          info = sprintf("y = %g, shape = %g, rate = %g", y, shape, rate)
        )
      }
    })
  }
  expect_identical(dist_lcdf(0, c(3, 1), 2), -Inf)
})

test_that("analytical gamma lcdf is finite and accurate in the lower tail", {
  cases <- list(
    # d > pwindow, so q > 0
    list(d = 2, pwindow = 1, p = c(400, 1 / 5)),
    list(d = 2, pwindow = 1, p = c(100, 1 / 5)),
    list(d = 1.5, pwindow = 0.5, p = c(60, 1 / 4)),
    list(d = 3, pwindow = 2, p = c(250, 1 / 6)),
    # d < pwindow, so q = 0
    list(d = 2, pwindow = 3, p = c(40, 1 / 10)),
    list(d = 1, pwindow = 4, p = c(12, 1 / 4)),
    # d equal to pwindow, q = 0
    list(d = 2, pwindow = 2, p = c(400, 1 / 5))
  )
  for (case in cases) {
    label <- sprintf(
      "d = %g, pwindow = %g, params = (%g, %g)",
      case$d, case$pwindow, case$p[1], case$p[2]
    )
    res <- primarycensored_analytical_lcdf(
      case$d, 2, case$p, case$pwindow, 0, Inf, 1, numeric(0)
    )
    expect_true(is.finite(res), info = label)
    expect_equal(
      res, ref_lcdf_unif_gamma(case$d, case$pwindow, case$p[1], case$p[2]),
      tolerance = 1e-8, info = label
    )
  }
})

test_that("gamma truncation normalisers stay finite when the CDF is deep in
   the lower tail", {
  params <- c(400, 1 / 5)
  pwindow <- 1
  d <- 2
  D <- 2.5
  ref <- function(x) ref_lcdf_unif_gamma(x, pwindow, 400, 1 / 5)
  expected <- ref(d) - ref(D)
  res <- primarycensored_analytical_lcdf(
    d, 2, params, pwindow, 0, D, 1, numeric(0)
  )
  expect_true(is.finite(res))
  expect_lt(res, 0)
  expect_equal(res, expected, tolerance = 1e-8)

  res_dispatch <- primarycensored_lcdf(
    d, 2, params, pwindow, 0, D, 1, numeric(0)
  )
  expect_equal(res_dispatch, expected, tolerance = 1e-8)

  # Lower and upper truncation together
  L <- 1.5
  log_diff <- function(a, b) a + log1p(-exp(b - a))
  expected_both <- log_diff(ref(d), ref(L)) - log_diff(ref(D), ref(L))
  res_both <- primarycensored_analytical_lcdf(
    d, 2, params, pwindow, L, D, 1, numeric(0)
  )
  expect_true(is.finite(res_both))
  expect_equal(res_both, expected_both, tolerance = 1e-8)
})

test_that("analytical gamma matches Stan's numerical path", {
  # Integrate the right hand side of `primarycensored_ode()` over
  # [d - pwindow, d]
  stan_numeric <- function(d, pwindow, params) {
    rhs <- function(t) {
      vapply(t, function(ti) {
        primarycensored_ode(ti, 0, params, c(d, pwindow), c(2L, 1L, 2L, 0L))
      }, numeric(1))
    }
    integrate(rhs, lower = d - pwindow, upper = d, rel.tol = 1e-12)$value
  }
  cases <- list(
    c(20, 1 / 5), c(12, 1 / 4), c(8, 1 / 3), c(40, 1 / 2), c(100, 1),
    c(2.5, 1), c(0.6, 2)
  )
  for (params in cases) {
    mean_delay <- params[1] / params[2]
    for (pwindow in c(0.5, 1, 3)) {
      for (d in mean_delay * c(0.3, 0.5, 0.7, 0.9, 1, 1.5, 3)) {
        label <- sprintf(
          "d = %g, pwindow = %g, params = (%s)", d, pwindow,
          toString(params)
        )
        analytic <- exp(primarycensored_analytical_lcdf(
          d, 2, params, pwindow, 0, Inf, 1, numeric(0)
        ))
        expect_equal(
          analytic, stan_numeric(d, pwindow, params),
          tolerance = 1e-6, info = label
        )
      }
    }
  }
})

test_that("Stan analytical gamma matches the R implementation", {
  for (shape in c(1, 2.5, 10, 30, 400)) {
    for (rate in c(0.2, 10)) {
      obj <- new_pcens(pgamma, dunif, list(), shape = shape, rate = rate)
      q <- shape / rate * c(0.3, 0.6, 0.9, 1, 1.1, 1.5, 3)
      for (pwindow in c(0.5, 1, 3)) {
        r_result <- pcens_cdf(obj, q = q, pwindow = pwindow)
        stan_result <- vapply(q, function(d) {
          exp(primarycensored_analytical_lcdf(
            d, 2, c(shape, rate), pwindow, 0, Inf, 1, numeric(0)
          ))
        }, numeric(1))
        expect_equal(
          stan_result, r_result,
          tolerance = 1e-6,
          info = sprintf(
            "shape = %g, rate = %g, pwindow = %g", shape, rate, pwindow
          )
        )
      }
    }
  }
})

test_that("gamma_lcdf_logx matches pgamma in the lower tail and bulk", {
  # frac = x / (a + 1) spans the series below 0.5 and gamma_lcdf() above
  for (a in c(0.01, 0.1, 1, 2, 3, 5, 9.5, 10, 20, 40, 100, 400)) {
    for (frac in c(
      1e-6, 1e-3, 0.1, 0.3, 0.49, 0.5, 0.51, 0.7, 0.9, 1, 1.2, 1.5, 3, 10
    )) {
      x <- frac * (a + 1)
      expect_lcdf_close(
        gamma_lcdf_logx(log(x), a), pgamma(x, a, log.p = TRUE),
        sprintf("a = %g, x over (a + 1) = %g", a, frac)
      )
    }
  }
  # Issue #363
  expect_equal(gamma_lcdf_logx(log(0.4), 400), -2367.416, tolerance = 1e-6)
})

test_that("gamma_lcdf_logx does not underflow when x does", {
  # exp(-1000) is 0 but log P is finite
  for (a in c(0.01, 1, 3, 30.5)) {
    expect_equal(
      gamma_lcdf_logx(-1000, a), a * -1000 - lgamma(a + 1),
      tolerance = 1e-12, info = sprintf("a = %g", a)
    )
  }
  expect_identical(gamma_lcdf_logx(-Inf, 3), -Inf)
  expect_identical(gamma_lcdf_logx(Inf, 3), 0)
})

test_that("the gamma log CDF helpers reject an invalid shape", {
  for (a in c(0, -0.5, -3, NaN, Inf)) {
    expect_error(gamma_lcdf_logx(log(2), a), "shape")
    expect_error(dist_lcdf(2, c(a, 1), 2), "shape")
  }
  expect_error(dist_lcdf(2, c(1.5, -1), 2))
})

test_that("the gamma log CDFs reject a rate that is not positive and finite", {
  for (rate in c(0, -1, Inf, -Inf, NaN)) {
    expect_error(dist_lcdf(1, c(2, rate), 2), "rate")
    expect_error(
      primarycensored_lcdf(1, 2L, c(2, rate), 1, 0, Inf, 1L, numeric(0)),
      "rate"
    )
    expect_error(
      primarycensored_gamma_uniform_terms(1, c(2, rate)), "rate"
    )
  }
})
