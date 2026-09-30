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
  expect_equal(
    dist_lcdf(2, c(400, 0.2), 2),
    pgamma(0.4, 400, log.p = TRUE),
    tolerance = 1e-9
  )
  for (i in seq_len(nrow(gamma_tail_cases))) {
    with(gamma_tail_cases[i, ], {
      for (y in c(0.05, 0.3, 1, 2, 3, 4.5)) {
        expect_equal(
          dist_lcdf(y, c(shape, rate), 2),
          ref_lgamma_delay(y, shape, rate),
          tolerance = 1e-9,
          info = sprintf("y = %g, shape = %g, rate = %g", y, shape, rate)
        )
      }
    })
  }
  expect_identical(dist_lcdf(0, c(3, 1), 2), -Inf)
})

test_that("dist_lcdf for a Gamma delay is accurate in the body and the
   upper tail", {
  # Shapes either side of 10, where the evaluation rule changes
  for (shape in c(0.3, 2, 9.5, 10, 10.5, 75, 400, 3000)) {
    mu <- shape / 1.5
    for (frac in c(0.3, 0.7, 0.95, 1, 1.05, 1.3, 2, 5, 50)) {
      expect_lcdf_close(
        dist_lcdf(frac * mu, c(shape, 1.5), 2),
        pgamma(frac * mu * 1.5, shape, log.p = TRUE),
        sprintf("shape = %g, y over mean = %g", shape, frac)
      )
    }
  }
})

test_that("gamma uniform terms are finite and accurate in the lower tail", {
  for (i in seq_len(nrow(gamma_tail_cases))) {
    with(gamma_tail_cases[i, ], {
      for (t in c(0.05, 0.3, 1, 2, 3, 4.5)) {
        terms <- primarycensored_gamma_uniform_terms(t, c(shape, rate))
        expected <- c(
          log(t) + ref_lgamma_delay(t, shape, rate),
          log(shape / rate) + pgamma(t * rate, shape + 1, log.p = TRUE)
        )
        expect_true(all(is.finite(terms)))
        expect_equal(
          terms, expected,
          tolerance = 1e-9,
          info = sprintf("t = %g, shape = %g, rate = %g", t, shape, rate)
        )
      }
    })
  }
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

test_that("Stan analytical gamma matches the R implementation for shapes
   of 10 or more", {
  for (shape in c(10, 30, 400, 1500)) {
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

test_that("Stan analytical gamma matches rprimarycensored samples", {
  set.seed(381)
  n <- 2e5
  for (case in list(list(30, 1), list(400, 2), list(2000, 1))) {
    shape <- case[[1]]
    rate <- case[[2]]
    pwindow <- 2
    samples <- rprimarycensored(
      n, rgamma,
      pwindow = pwindow, swindow = 0, shape = shape, rate = rate
    )
    probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
    for (d in unname(quantile(samples, probs))) {
      empirical <- mean(samples <= d)
      analytic <- exp(primarycensored_lcdf(
        d, 2L, c(shape, rate), pwindow, 0, Inf, 1L, numeric(0)
      ))
      expect_lt(
        abs(analytic - empirical),
        4 * sqrt(empirical * (1 - empirical) / n),
        label = sprintf("shape = %g, rate = %g, d = %g", shape, rate, d)
      )
    }
  }
})

test_that("analytical gamma lcdf is accurate for large shapes", {
  cases <- list(
    list(d = 1400, pwindow = 10, p = c(1500, 1)),
    list(d = 1510, pwindow = 10, p = c(1500, 1)),
    list(d = 1600, pwindow = 10, p = c(1500, 1)),
    list(d = 2900, pwindow = 10, p = c(3000, 1)),
    list(d = 3100, pwindow = 5, p = c(3000, 1))
  )
  for (case in cases) {
    label <- sprintf(
      "d = %g, pwindow = %g, params = (%g, %g)",
      case$d, case$pwindow, case$p[1], case$p[2]
    )
    res <- primarycensored_analytical_lcdf(
      case$d, 2, case$p, case$pwindow, 0, Inf, 1, numeric(0)
    )
    expect_equal(
      res, ref_lcdf_unif_gamma(case$d, case$pwindow, case$p[1], case$p[2]),
      tolerance = 1e-7, info = label
    )
  }
})

test_that("gamma_lcdf_logx matches pgamma across the evaluation rules", {
  # frac = x / (a + 1) spans the lower tail, the body and the upper tail
  for (a in c(0.1, 1, 5, 9.5, 9.999, 10, 10.001, 40, 400, 4000, 30000, 1e5)) {
    for (frac in c(
      1e-6, 1e-3, 0.1, 0.3, 0.5, 0.7, 0.89, 0.9, 0.99, 1, 1.01, 1.05, 1.2,
      1.5, 3, 10, 100
    )) {
      x <- frac * (a + 1)
      expect_lcdf_close(
        gamma_lcdf_logx(log(x), a),
        pgamma(x, a, log.p = TRUE),
        sprintf("a = %g, x over (a + 1) = %g", a, frac)
      )
    }
  }
})

test_that("gamma_lcdf_logx is accurate either side of the rule changes", {
  # x = a + 1 switches the series and the continued fraction, and a = 10
  # switches `gamma_lcdf`
  for (a in c(1.5, 9.99, 10.01, 150.5)) {
    for (x in c(a + 1 - 1e-7, a + 1 + 1e-7)) {
      expect_lcdf_close(
        gamma_lcdf_logx(log(x), a),
        pgamma(x, a, log.p = TRUE),
        sprintf("a = %g, x = %.9g", a, x),
        tolerance = 1e-10
      )
    }
  }
})

test_that("gamma_lcdf_logx is exact for integer shapes", {
  # The continued fraction terminates at i = a
  for (a in c(2, 20, 21, 100, 1000)) {
    for (frac in c(1.01, 1.5, 4)) {
      x <- frac * (a + 1)
      expect_equal(
        gamma_lcdf_logx(log(x), a),
        pgamma(x, a, log.p = TRUE),
        tolerance = 1e-9
      )
    }
  }
})

test_that("gamma_lcdf_logx does not underflow when x does", {
  # exp(-800) underflows to 0 but the log CDF is finite
  expect_equal(
    gamma_lcdf_logx(-800, 30.5), 30.5 * -800 - lgamma(31.5),
    tolerance = 1e-12
  )
  expect_identical(gamma_lcdf_logx(-Inf, 30.5), -Inf)
})

test_that("gamma_lcdf_logx_pair matches pgamma for a and a + 1", {
  # frac = x / (a + 1) covers the series, the continued fraction and
  # `gamma_lcdf`
  for (a in c(0.05, 0.5, 1, 5, 9.5, 9.999, 10, 10.001, 40, 400, 4000, 3e4)) {
    for (frac in c(
      1e-8, 1e-4, 0.01, 0.1, 0.3, 0.49, 0.5, 0.51, 0.7, 0.89, 0.9, 0.99, 1,
      1.01, 1.2, 1.5, 3, 10, 100
    )) {
      x <- frac * (a + 1)
      pair <- gamma_lcdf_logx_pair(log(x), a)
      label <- sprintf("a = %g, x over (a + 1) = %g", a, frac)
      expect_length(pair, 2)
      expect_lcdf_close(
        pair[1], pgamma(x, a, log.p = TRUE), label
      )
      expect_lcdf_close(
        pair[2], pgamma(x, a + 1, log.p = TRUE), label
      )
    }
  }
})

test_that("gamma_lcdf_logx_pair does not underflow when x does", {
  # x = exp(-800) is 0 but both log CDFs are finite
  for (a in c(0.5, 9.5, 30.5)) {
    expect_equal(
      gamma_lcdf_logx_pair(-800, a),
      c(a * -800 - lgamma(a + 1), (a + 1) * -800 - lgamma(a + 2)),
      tolerance = 1e-12
    )
  }
  expect_identical(gamma_lcdf_logx_pair(-Inf, 3), c(-Inf, -Inf))
  expect_identical(gamma_lcdf_logx_pair(Inf, 3), c(0, 0))
})

test_that("gamma_lcdf_logx_pair keeps P(a + 1) accurate when it is far
   below P(a)", {
  # x much less than a + 1, where the recursion would cancel
  for (a in c(0.2, 3, 8)) {
    for (x in c(1e-12, 1e-9, 1e-6, 1e-3)) {
      expect_equal(
        gamma_lcdf_logx_pair(log(x), a)[2],
        pgamma(x, a + 1, log.p = TRUE),
        tolerance = 1e-12,
        info = sprintf("a = %g, x = %g", a, x)
      )
    }
  }
})

test_that("gamma_lcdf_logx and gamma_lcdf_logx_pair reject an invalid shape", {
  for (a in c(0, -0.5, -3)) {
    expect_error(gamma_lcdf_logx(log(2), a), "shape")
    expect_error(gamma_lcdf_logx_pair(log(2), a), "shape")
    expect_error(dist_lcdf(2, c(a, 1), 2), "shape")
  }
  expect_error(gamma_lcdf_logx(log(2), NaN), "shape")
  expect_error(gamma_lcdf_logx(log(2), Inf), "shape")
  expect_error(gamma_lcdf_logx_pair(log(2), Inf), "shape")
  expect_error(dist_lcdf(2, c(Inf, 1), 2), "shape")
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

test_that("gamma_log_lead_logx agrees with the Poisson log PMF for large
   shapes", {
  # dpois() is accurate near x = a
  for (a in c(50, 99, 100, 150, 1e3, 1e4, 1e5, 1e6, 1e7)) {
    for (r in c(1e-6, 0.01, 0.3, 0.5, 0.6, 0.9, 0.999, 1, 1.001, 1.5,
                1.99, 2.1, 5, 30)) {
      log_x <- log(a * r)
      expected <- dpois(a, exp(log_x), log = TRUE)
      expect_lt(
        abs(gamma_log_lead_logx(log_x, a) - expected),
        1e-12 * max(1, abs(expected)),
        label = sprintf("a = %g, x over a = %g", a, r)
      )
    }
  }
})

test_that("primarycensored_lcdf is accurate for a large Gamma shape when
   the delay is long compared with the primary window", {
  # Rounding error is amplified by about d / pwindow
  cases <- list(
    list(d = 1e7, shape = 1e5, rate = 0.01, pwindow = 0.01, tol = 1e-5),
    list(d = 1e6, shape = 1e4, rate = 0.01, pwindow = 0.01, tol = 1e-6),
    list(d = 2e6, shape = 1e5, rate = 0.05, pwindow = 0.5, tol = 1e-7),
    list(d = 1e5, shape = 1e5, rate = 1, pwindow = 1, tol = 1e-8)
  )
  for (case in cases) {
    expected <- ref_lcdf_unif_gamma(
      case$d, case$pwindow, case$shape, case$rate
    )
    actual <- primarycensored_lcdf(
      case$d, 2L, c(case$shape, case$rate), case$pwindow, 0, Inf, 1L,
      numeric(0)
    )
    expect_lt(
      abs(actual - expected), case$tol,
      label = sprintf(
        "d = %g, shape = %g, rate = %g, pwindow = %g", case$d, case$shape,
        case$rate, case$pwindow
      )
    )
  }
})
