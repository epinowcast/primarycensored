skip_on_cran()

# Regression tests for #363. `gamma_lcdf` underflows to -inf deep in the
# lower tail of the generalised gamma, which gave NaN or +Inf log densities
# once the truncation normaliser was subtracted. Stan's partial of the
# regularised incomplete gamma with respect to its shape is also inaccurate
# well before that point.
#
# The reference is R's `pgamma(log.p = TRUE)`, which is accurate in the
# tails, so the tolerance is relative 1e-9 unless stated.

# log F_T(t) for the generalised gamma (Stacy parameterisation) via R
ref_lgengamma <- function(t, shape, scale, k) {
  pgamma((t / scale)^shape, shape = k, log.p = TRUE)
}

# log of the uniform primary event censored CDF, from
# F_{S+}(d) = (1 / w) int_{max(d - w, 0)}^{d} F_T(u) du. The integrand is
# scaled by its maximum at u = d so that nothing underflows. The range is
# cut where it has fallen by a factor of exp(-40), which only matters in the
# lower tail.
ref_lcdf_unif <- function(d, pwindow, shape, scale, k) {
  q_lo <- max(d - pwindow, 0)
  log_max <- ref_lgengamma(d, shape, scale, k)
  target <- log_max - 40
  lower <- q_lo
  if (ref_lgengamma(q_lo, shape, scale, k) < target) {
    lower <- uniroot(
      function(u) ref_lgengamma(u, shape, scale, k) - target,
      lower = q_lo, upper = d, tol = 1e-14
    )$root
  }
  scaled <- integrate(
    function(u) exp(ref_lgengamma(u, shape, scale, k) - log_max),
    lower = lower, upper = d, rel.tol = 1e-13, subdivisions = 1000L
  )$value
  log_max + log(scaled) - log(pwindow)
}

tail_cases <- data.frame(
  shape = c(1, 5, 1.5, 0.7, 2, 3),
  scale = c(5, 5, 3, 4, 6, 10),
  k = c(400, 100, 30, 60, 250, 40)
)

test_that("gengamma_lcdf is finite and accurate deep in the lower tail", {
  # Examples from #363
  expect_equal(
    gengamma_lcdf(2, 1, 5, 400),
    pgamma(0.4, 400, log.p = TRUE),
    tolerance = 1e-9
  )
  expect_equal(gengamma_lcdf(2, 1, 5, 400), -2367.4, tolerance = 1e-4)
  expect_equal(gengamma_lcdf(2, 5, 5, 100), -821.9, tolerance = 1e-4)

  for (i in seq_len(nrow(tail_cases))) {
    with(tail_cases[i, ], {
      for (y in c(0.05, 0.3, 1, 2, 3, 4.5)) {
        expect_equal(
          gengamma_lcdf(y, shape, scale, k),
          ref_lgengamma(y, shape, scale, k),
          tolerance = 1e-9,
          info = sprintf(
            "y = %g, shape = %g, scale = %g, k = %g", y, shape, scale, k
          )
        )
      }
    })
  }
})

test_that("gengamma_lcdf agrees with flexsurv in the tails", {
  skip_if_not_installed("flexsurv")
  for (y in c(0.1, 1, 2, 4, 8, 20)) {
    expect_equal(
      gengamma_lcdf(y, 1.5, 3, 12),
      flexsurv::pgengamma.orig(y,
        shape = 1.5, scale = 3, k = 12,
        log.p = TRUE
      ),
      tolerance = 1e-9
    )
  }
})

test_that("gengamma_lcdf does not underflow when the power does", {
  # (y / scale)^shape underflows to 0 here, but the log CDF is finite
  y <- 1e-7
  shape <- 50
  scale <- 1
  k <- 2
  log_x <- shape * log(y / scale)
  expect_identical((y / scale)^shape, 0)
  # P(k, x) = x^k e^{-x} / Gamma(k + 1) * (1 + x / (k + 1) + ...)
  expected <- k * log_x - lgamma(k + 1)
  expect_equal(gengamma_lcdf(y, shape, scale, k), expected, tolerance = 1e-12)
})

test_that("gamma_lcdf_logx matches pgamma across the series rules", {
  # a spans the body and the extreme shapes where only the series is finite.
  # frac = x / (a + 1) covers both sides of each switch in the rule.
  for (a in c(0.1, 1, 5, 40, 400, 4000, 8000, 30000)) {
    for (frac in c(
      1e-6, 1e-3, 0.1, 0.3, 0.49, 0.5, 0.51, 0.7, 0.89, 0.9,
      0.95, 1.2, 3
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
  expect_identical(gamma_lcdf_logx(-Inf, 2), -Inf)
  # x underflows to 0 but the log CDF is finite
  expect_equal(
    gamma_lcdf_logx(-800, 3), 3 * -800 - lgamma(4),
    tolerance = 1e-12
  )
})

test_that("gengamma_lcdf is continuous across the tail rule", {
  # Values either side of any switch between evaluation rules must agree
  for (k in c(3, 25, 150)) {
    for (frac in c(0.3, 0.45, 0.5, 0.55, 0.7)) {
      x <- frac * (k + 1)
      eps <- 1e-7
      lo <- gengamma_lcdf((x * (1 - eps))^(1 / 1.3) * 2, 1.3, 2, k)
      hi <- gengamma_lcdf((x * (1 + eps))^(1 / 1.3) * 2, 1.3, 2, k)
      expect_lt(abs(hi - lo), 1e-4)
      expect_equal(
        gengamma_lcdf(x^(1 / 1.3) * 2, 1.3, 2, k),
        pgamma(x, k, log.p = TRUE),
        tolerance = 1e-9
      )
    }
  }
})

test_that("gengamma_lcdf is -inf at zero", {
  expect_identical(gengamma_lcdf(0, 1.5, 2, 3), -Inf)
})

test_that("analytical generalised gamma lcdf is finite and accurate in the
   lower tail", {
  cases <- list(
    # d > pwindow, so q > 0
    list(d = 2, pwindow = 1, p = c(1, 5, 400)),
    list(d = 2, pwindow = 1, p = c(5, 5, 100)),
    list(d = 1.5, pwindow = 0.5, p = c(1.5, 3, 30)),
    list(d = 3, pwindow = 2, p = c(0.7, 4, 60)),
    # d < pwindow, so q = 0
    list(d = 2, pwindow = 3, p = c(2, 6, 250)),
    list(d = 1, pwindow = 4, p = c(3, 10, 40)),
    # d equal to pwindow, q = 0
    list(d = 2, pwindow = 2, p = c(1, 5, 400))
  )
  for (case in cases) {
    shape <- case$p[1]
    scale <- case$p[2]
    k <- case$p[3]
    label <- sprintf(
      "d = %g, pwindow = %g, params = (%g, %g, %g)",
      case$d, case$pwindow, shape, scale, k
    )
    res <- primarycensored_analytical_lcdf(
      case$d, 5, case$p, case$pwindow, 0, Inf, 1, numeric(0)
    )
    expect_true(is.finite(res), info = label)
    expect_equal(
      res, ref_lcdf_unif(case$d, case$pwindow, shape, scale, k),
      tolerance = 1e-8, info = label
    )
  }
})

test_that("truncation normalisers stay finite when the CDF is deep in the
   lower tail", {
  params <- c(1, 5, 400)
  pwindow <- 1
  d <- 2
  # Upper truncation D also deep in the tail: log_cdf_D was -inf, so the
  # normalised result was NaN or +Inf.
  D <- 2.5
  expected <- ref_lcdf_unif(d, pwindow, 1, 5, 400) -
    ref_lcdf_unif(D, pwindow, 1, 5, 400)
  res <- primarycensored_analytical_lcdf(
    d, 5, params, pwindow, 0, D, 1, numeric(0)
  )
  expect_true(is.finite(res))
  expect_lt(res, 0)
  expect_equal(res, expected, tolerance = 1e-8)

  res_dispatch <- primarycensored_lcdf(
    d, 5, params, pwindow, 0, D, 1, numeric(0)
  )
  expect_equal(res_dispatch, expected, tolerance = 1e-8)

  # Lower and upper truncation together
  L <- 1.5
  lo <- ref_lcdf_unif(L, pwindow, 1, 5, 400)
  hi <- ref_lcdf_unif(D, pwindow, 1, 5, 400)
  mid <- ref_lcdf_unif(d, pwindow, 1, 5, 400)
  # Log of (F(d) - F(L)) / (F(D) - F(L))
  log_diff <- function(a, b) a + log1p(-exp(b - a))
  expected_both <- log_diff(mid, lo) - log_diff(hi, lo)
  res_both <- primarycensored_analytical_lcdf(
    d, 5, params, pwindow, L, D, 1, numeric(0)
  )
  expect_true(is.finite(res_both))
  expect_equal(res_both, expected_both, tolerance = 1e-8)
})

test_that("analytical generalised gamma still matches the numeric path
   where the CDF is representable", {
  for (params in list(c(2, 5, 20), c(1, 4, 12), c(1.5, 3, 8))) {
    for (pwindow in c(0.5, 1, 3)) {
      for (d in c(0.4, 1, 2.5, 6, 12, 25)) {
        label <- sprintf(
          "d = %g, pwindow = %g, params = (%s)", d, pwindow,
          toString(params)
        )
        analytic <- primarycensored_analytical_lcdf(
          d, 5, params, pwindow, 0, Inf, 1, numeric(0)
        )
        # The numeric path integrates the CDF, so only compare where the
        # result is not vanishingly small
        if (analytic > -30) {
          numeric <- log(primarycensored_cdf(
            d, 5, params, pwindow, 0, Inf, 1, numeric(0)
          ))
          expect_equal(analytic, numeric, tolerance = 1e-6, info = label)
        }
      }
    }
  }
})
