skip_on_cran()

# Regression tests for #381. The Gamma delay (dist_id 2) called
# `gamma_lcdf`, which underflows to -inf deep in the lower tail. The
# uniform terms then subtracted two -inf values with `log_diff_exp`, giving
# NaN. The same fix as for the generalised gamma in #363 applies, via
# `gamma_lcdf_logx()`.
#
# The reference is R's `pgamma(log.p = TRUE)`, which is accurate in the
# tails, so the tolerance is relative 1e-9 unless stated.

# log F_T(t) for a Gamma delay via R
ref_lgamma_delay <- function(t, shape, rate) {
  pgamma(t * rate, shape = shape, log.p = TRUE)
}

# log of the uniform primary event censored CDF, from
# F_{S+}(d) = (1 / w) int_{max(d - w, 0)}^{d} F_T(u) du. The integrand is
# scaled by its maximum at u = d so that nothing underflows. The range is
# cut where it has fallen by a factor of exp(-40), which only matters in the
# lower tail.
ref_lcdf_unif_gamma <- function(d, pwindow, shape, rate) {
  q_lo <- max(d - pwindow, 0)
  log_max <- ref_lgamma_delay(d, shape, rate)
  target <- log_max - 40
  lower <- q_lo
  if (ref_lgamma_delay(q_lo, shape, rate) < target) {
    lower <- uniroot(
      function(u) ref_lgamma_delay(u, shape, rate) - target,
      lower = q_lo, upper = d, tol = 1e-14
    )$root
  }
  scaled <- integrate(
    function(u) exp(ref_lgamma_delay(u, shape, rate) - log_max),
    lower = lower, upper = d, rel.tol = 1e-13, subdivisions = 1000L
  )$value
  log_max + log(scaled) - log(pwindow)
}

# Shape, rate and a delay that puts the log CDF well below -10, the
# smallest at -2400
gamma_tail_cases <- data.frame(
  shape = c(400, 100, 60, 250, 40, 12),
  rate = c(1 / 5, 1 / 5, 1 / 4, 1 / 6, 1 / 10, 1 / 4)
)

test_that("dist_lcdf for a Gamma delay is finite and accurate deep in the
   lower tail", {
  # log F_T(2) for shape 400 and rate 1 / 5 is about -2367
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

test_that("dist_lcdf for a Gamma delay is unchanged in the body and the
   upper tail", {
  # shape above and below the point where the evaluation rule changes
  for (shape in c(0.3, 2, 19.5, 20, 20.5, 75, 400, 3000)) {
    mu <- shape / 1.5
    for (frac in c(0.3, 0.7, 0.95, 1, 1.05, 1.3, 2, 5, 50)) {
      expect_equal(
        dist_lcdf(frac * mu, c(shape, 1.5), 2),
        pgamma(frac * mu * 1.5, shape, log.p = TRUE),
        tolerance = 1e-9,
        info = sprintf("shape = %g, y over mean = %g", shape, frac)
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
  # The Stan numerical path integrates exp(dist_lcdf(t)) / pwindow over
  # [d - pwindow, d] with `primarycensored_ode()`. Integrate that same
  # right hand side so the analytical path is checked against the delay CDF
  # that the numerical path uses, from deep in the lower tail to the upper
  # tail.
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

test_that("analytical gamma lcdf is accurate for large shapes", {
  # Shapes of 1000 or more, where the CDF is a narrow step
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
  # a spans the body, the point where the rule changes at 20, and the large
  # shapes. frac = x / (a + 1) covers the lower tail, the body and the upper
  # tail on both sides of x = a + 1.
  for (a in c(0.1, 1, 5, 19.5, 19.999, 20, 20.001, 40, 400, 4000, 30000, 1e5)) {
    for (frac in c(
      1e-6, 1e-3, 0.1, 0.3, 0.5, 0.7, 0.89, 0.9, 0.99, 1, 1.01, 1.05, 1.2,
      1.5, 3, 10, 100
    )) {
      x <- frac * (a + 1)
      expect_equal(
        gamma_lcdf_logx(log(x), a),
        pgamma(x, a, log.p = TRUE),
        tolerance = 1e-9,
        info = sprintf("a = %g, x over (a + 1) = %g", a, frac)
      )
    }
  }
})

test_that("gamma_lcdf_logx is continuous across the rule changes", {
  # x = a + 1 switches between the series and the continued fraction, and
  # a = 20 between Stan's gamma_lcdf and the series and fraction
  for (a in c(1.5, 19.99, 20.01, 150.5)) {
    for (x in c(a + 1 - 1e-7, a + 1 + 1e-7)) {
      expect_equal(
        gamma_lcdf_logx(log(x), a),
        pgamma(x, a, log.p = TRUE),
        tolerance = 1e-10
      )
    }
  }
  x <- 30
  expect_equal(
    gamma_lcdf_logx(log(x), 19.9999999),
    gamma_lcdf_logx(log(x), 20.0000001),
    tolerance = 1e-7
  )
})

test_that("gamma_lcdf_logx is exact for integer shapes", {
  # The continued fraction terminates at i = a, so its derivative needs the
  # rest of the fraction. The value is checked here, gradients elsewhere.
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
  # exp(-800) underflows to 0 here, but the log CDF is finite. The
  # relative error of 1e-12 is the rounding error of the leading term.
  expect_equal(
    gamma_lcdf_logx(-800, 30.5), 30.5 * -800 - lgamma(31.5),
    tolerance = 1e-12
  )
  expect_identical(gamma_lcdf_logx(-Inf, 30.5), -Inf)
})
