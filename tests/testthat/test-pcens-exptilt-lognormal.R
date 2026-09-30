# The lognormal delay with an exponentially tilted primary. The tilt
# transform has no closed form, so it is evaluated by quadrature for a tilt
# rho = -xi > 0 and by a series of partial moments for rho < 0. These tests
# check the transform against a reference integral, and the CDF against the
# reference integral of [exptilt_reference()] and the numerical path.

cases <- exptilt_lnorm_cases()

lnorm_object <- function(case, rho) {
  exptilt_object(exptilt_lnorm_family(case), rho)
}

test_that("lognormal delays with a tilted primary dispatch to the transform", {
  obj <- lnorm_object(cases[[1]], 0.2)
  expect_s3_class(obj, "pcens_plnorm_dexpgrowth")
  expect_false(
    is.null(
      utils::getS3method(
        "pcens_cdf", "pcens_plnorm_dexpgrowth",
        optional = TRUE
      )
    )
  )
  expect_identical(.pcens_tilt_lower(obj), 0)
  expect_true(.pcens_tilt_available(obj, -0.2))
  expect_true(.pcens_tilt_available(obj, 0))
  expect_true(.pcens_tilt_available(obj, 0.7))
})

test_that("the lognormal tilt is unavailable where the mode overflows", {
  obj <- new_pcens(
    plnorm, dexpgrowth, list(r = 1e300),
    meanlog = 650, sdlog = 1
  )
  expect_false(.pcens_tilt_available(obj, -1e300))
  obj <- new_pcens(plnorm, dexpgrowth, list(r = 1), meanlog = 0, sdlog = 0)
  expect_false(.pcens_tilt_available(obj, -1))
  obj <- new_pcens(plnorm, dexpgrowth, list(r = 1), meanlog = NA, sdlog = 1)
  expect_false(.pcens_tilt_available(obj, -1))
})

test_that("the lognormal transform matches a reference integral", {
  z <- c(-12, -6, -3, -1, 0, 0.5, 1, 2, 3)
  xis <- c(-5, -1, -0.25, -0.01, -1e-6, 0, 1e-3, 0.1, 0.5, 1)
  for (case in cases) {
    obj <- lnorm_object(case, 0.1)
    t <- exp(case$meanlog + case$sdlog * z)
    for (xi in xis) {
      info <- sprintf(
        "meanlog %g, sdlog %g, xi %g", case$meanlog, case$sdlog, xi
      )
      lower <- lnorm_tilt_reference(t, case$meanlog, case$sdlog, xi)
      actual <- .pcens_tilt_transform(obj, t, xi)
      # Values that underflow to a few hundred are not compared
      keep <- lower > -700
      expect_lt(max(abs(actual[keep] - lower[keep])), 1e-9, label = info)
      expect_identical(actual[!is.finite(lower)], lower[!is.finite(lower)])
      if (xi <= 0) {
        upper <- lnorm_tilt_reference(t, case$meanlog, case$sdlog, xi, TRUE)
        actual <- .pcens_tilt_transform(obj, t, xi, upper = TRUE)
        keep <- upper > -700
        expect_lt(max(abs(actual[keep] - upper[keep])), 1e-9, label = info)
      } else {
        # The total diverges for a positive tilt
        expect_identical(
          .pcens_tilt_transform(obj, t, xi, upper = TRUE),
          rep(Inf, length(t))
        )
      }
    }
  }
})

test_that("the pair of transforms is the lower and upper transform", {
  obj <- lnorm_object(cases[[2]], 0.2)
  t <- c(-1, 0, 0.05, 0.8, 2, 7, 40)
  for (xi in c(-1, -0.2, 0, 0.3)) {
    pair <- .pcens_tilt_pair(obj, t, xi)
    expect_identical(colnames(pair), c("lower", "upper"))
    expect_identical(pair[, 1], .pcens_tilt_transform(obj, t, xi))
    expect_identical(pair[, 2], .pcens_tilt_transform(obj, t, xi, upper = TRUE))
  }
  # The default method of the generic gives the same from the two transforms
  gamma_obj <- exptilt_object(exptilt_families()[[3]], 0.2)
  pair <- .pcens_tilt_pair(gamma_obj, t, -0.2)
  expect_identical(pair[, 1], .pcens_tilt_transform(gamma_obj, t, -0.2))
  expect_identical(
    pair[, 2], .pcens_tilt_transform(gamma_obj, t, -0.2, upper = TRUE)
  )
})

test_that("the lower and upper lognormal transforms sum to one total", {
  for (case in cases) {
    obj <- lnorm_object(case, 0.2)
    t <- exp(case$meanlog + case$sdlog * c(-6, -2, -0.3, 0, 0.4, 1.5, 4))
    for (xi in c(-3, -0.4, -1e-3)) {
      pair <- .pcens_tilt_pair(obj, t, xi)
      total <- exp(pair[, 1]) + exp(pair[, 2])
      expect_equal(total, rep(total[1], length(t)), tolerance = 1e-9)
    }
  }
})

test_that("the lognormal transform at xi = 0 is the CDF", {
  obj <- lnorm_object(cases[[4]], 0.2)
  t <- c(1e-6, 0.3, 1, 5, 80)
  pair <- .pcens_tilt_pair(obj, t, 0)
  expect_equal(pair[, 1], plnorm(t, -1, 1.5, log.p = TRUE), tolerance = 1e-14)
  expect_equal(
    pair[, 2], plnorm(t, -1, 1.5, lower.tail = FALSE, log.p = TRUE),
    tolerance = 1e-14
  )
})

test_that("the lognormal transform is zero or total below the support", {
  obj <- lnorm_object(cases[[1]], 0.2)
  expect_identical(.pcens_tilt_transform(obj, c(-2, 0), -0.3), c(-Inf, -Inf))
  expect_identical(.pcens_tilt_transform(obj, c(-2, 0), 0.3), c(-Inf, -Inf))
  upper <- .pcens_tilt_transform(obj, c(-2, 0, 1e-300), -0.3, upper = TRUE)
  expect_identical(upper[1], upper[2])
  expect_equal(upper[3], upper[1], tolerance = 1e-12)
  expected <- lnorm_tilt_integral(-45, 45, 1.6, 0.5, -0.3)
  expect_equal(upper[1], expected, tolerance = 1e-10)
})

test_that("the lognormal moments match numerical integrals", {
  for (case in cases) {
    obj <- lnorm_object(case, 0.1)
    t <- exp(case$meanlog + case$sdlog * c(-3, -1, 0, 1, 2.5))
    moments <- .pcens_tilt_moments(obj, t)
    f <- function(u) dlnorm(u, case$meanlog, case$sdlog)
    for (k in 1:2) {
      expected <- vapply(t, function(tt) {
        stats::integrate(
          function(u) (tt - u)^k * f(u), 0, tt,
          rel.tol = 1e-12, abs.tol = 0, subdivisions = 500L
        )$value
      }, numeric(1))
      expect_equal(exp(moments[, k]), expected, tolerance = 1e-7)
    }
  }
  obj <- lnorm_object(cases[[1]], 0.1)
  moments <- .pcens_tilt_moments(obj, c(-1, 0))
  expect_identical(moments[, 1], c(-Inf, -Inf))
  expect_identical(moments[, 2], c(-Inf, -Inf))
})

test_that("the lognormal series agrees with quadrature for a small tilt", {
  # The two methods cover different signs, so compare the transform across
  # zero by the first order change in the tilt
  obj <- lnorm_object(cases[[1]], 0.1)
  t <- c(0.5, 3, 9, 30)
  lower_minus <- .pcens_tilt_transform(obj, t, -1e-7)
  lower_plus <- .pcens_tilt_transform(obj, t, 1e-7)
  expect_equal(lower_minus, lower_plus, tolerance = 1e-5)
  expect_equal(
    0.5 * (lower_minus + lower_plus), plnorm(t, 1.6, 0.5, log.p = TRUE),
    tolerance = 1e-8
  )
})

test_that("the lognormal series needs many terms for a large tilt times t", {
  obj <- lnorm_object(cases[[2]], 0.1)
  t <- c(1, 50, 400)
  lower <- .pcens_tilt_transform(obj, t, 0.5)
  expected <- lnorm_tilt_reference(t, 0, 1, 0.5)
  expect_lt(max(abs(lower - expected)), 1e-9)
})

test_that("the lognormal CDF matches a reference integral", {
  pwindows <- c(0.5, 1, 2, 7)
  rhos <- c(
    -1, -0.3, -0.05, -3e-4, -1e-4, -1e-5, -1e-8, 1e-8, 1e-5, 1e-4, 3e-4,
    0.05, 0.3, 1, 4
  )
  for (case in cases) {
    family <- exptilt_lnorm_family(case)
    for (pwindow in pwindows) {
      q <- sort(c(
        1e-6, 1e-3, 0.3 * pwindow, pwindow - 1e-3, pwindow, pwindow + 1e-3,
        2, 3, 6, 12, 25, 60
      ))
      for (rho in rhos) {
        obj <- exptilt_object(family, rho)
        expected <- exptilt_reference(
          q, pwindow, rho, exptilt_lnorm_cdf(case)
        )
        actual <- .pcens_cdf_exptilt(obj, q, pwindow)
        # Deep in the lower tail the direct form cancels, by about a factor
        # 1 / (rho q) times the gap between q and the mean of the delays
        # below q. The error is below 1e-7 in relative terms.
        expect_lt(
          max_rel_diff(actual, expected), 1e-7,
          label = exptilt_label(family, pwindow, rho)
        )
      }
    }
  }
})

test_that("the lognormal CDF agrees with use_numeric = TRUE", {
  q <- c(0.05, 0.5, 1.5, 3, 6, 12, 20)
  for (case in cases[1:4]) {
    for (pwindow in c(1, 2, 7)) {
      for (rho in c(-1, -0.5, -1e-8, 1e-8, 0.5, 1)) {
        obj <- lnorm_object(case, rho)
        analytic <- .pcens_cdf_exptilt(obj, q, pwindow)
        numeric <- pcens_cdf(obj, q, pwindow, use_numeric = TRUE)
        expect_equal(
          analytic, numeric,
          tolerance = 1e-6,
          info = sprintf(
            "meanlog %g, sdlog %g, pwindow = %g, r = %g",
            case$meanlog, case$sdlog, pwindow, rho
          )
        )
      }
    }
  }
})

test_that("the lognormal CDF is continuous in the tilt through zero", {
  for (case in cases) {
    family <- exptilt_lnorm_family(case)
    for (pwindow in c(0.5, 2, 7)) {
      q <- c(1e-3, 0.3 * pwindow, pwindow, 3, 6, 12)
      uniform <- exptilt_reference(
        q, pwindow, 0, exptilt_lnorm_cdf(case)
      )
      obj_zero <- exptilt_object(family, 0)
      expect_lt(
        max_rel_diff(.pcens_cdf_exptilt(obj_zero, q, pwindow), uniform), 1e-9
      )
      for (sign in c(-1, 1)) {
        for (rho in sign * c(1e-12, 1e-9, 1e-7)) {
          obj <- exptilt_object(family, rho)
          expect_lt(
            max_rel_diff(.pcens_cdf_exptilt(obj, q, pwindow), uniform), 1e-6,
            label = exptilt_label(family, pwindow, rho)
          )
        }
      }
    }
  }
})

test_that("the lognormal CDF has no jump where the small tilt form ends", {
  for (case in cases) {
    family <- exptilt_lnorm_family(case)
    for (pwindow in c(0.5, 2, 7)) {
      q <- c(1e-3, 0.3 * pwindow, pwindow, 3, 6, 12, 25)
      for (sign in c(-1, 1)) {
        below <- exptilt_object(family, sign * 0.9999e-4 / pwindow)
        above <- exptilt_object(family, sign * 1.0001e-4 / pwindow)
        expect_lt(
          max_rel_diff(
            .pcens_cdf_exptilt(below, q, pwindow),
            .pcens_cdf_exptilt(above, q, pwindow)
          ),
          1e-7,
          label = exptilt_label(family, pwindow, sign * 1e-4 / pwindow)
        )
      }
    }
  }
})

test_that("the lognormal CDF has no jump where series and quadrature meet", {
  # The sign of the tilt chooses the method of the transform. At a window
  # where |rho| w is 1e-3 the small tilt forms are not used.
  for (case in cases) {
    family <- exptilt_lnorm_family(case)
    pwindow <- 2
    q <- c(0.3, 1, 2, 3, 6, 12, 25)
    expected <- exptilt_reference(q, pwindow, 0, exptilt_lnorm_cdf(case))
    for (rho in c(-1e-3, 1e-3)) {
      obj <- exptilt_object(family, rho)
      expect_lt(
        max_rel_diff(.pcens_cdf_exptilt(obj, q, pwindow), expected), 5e-3,
        label = exptilt_label(family, pwindow, rho)
      )
    }
    # Both signs are close to the uniform window to first order in rho
    plus <- .pcens_cdf_exptilt(exptilt_object(family, 1e-3), q, pwindow)
    minus <- .pcens_cdf_exptilt(exptilt_object(family, -1e-3), q, pwindow)
    expect_lt(max_rel_diff(0.5 * (plus + minus), expected), 1e-5)
  }
})

test_that("the lognormal CDF is accurate far from the origin", {
  # The small window form must not be used where |r| q^2 / w is large
  pwindow <- 1
  for (m in c(1e5, 1e6, 1e7, 1e8)) {
    q <- m * seq(0.5, 2, length.out = 12)
    delay <- list(
      pdist = plnorm, args = list(meanlog = log(m), sdlog = 0.5)
    )
    for (rho in c(3e-5, -3e-5, 1e-6)) {
      obj <- exptilt_object(delay, rho)
      expected <- exptilt_reference(q, pwindow, rho, exptilt_cdf(delay))
      expect_lt(
        max_rel_diff(pcens_cdf(obj, q, pwindow), expected), 1e-7,
        label = sprintf("m = %g, r = %g", m, rho)
      )
    }
  }
  # With a window that is small next to the delay the CDF is the delay CDF
  delay <- list(pdist = plnorm, args = list(meanlog = 30, sdlog = 1))
  q <- exp(30) * seq(0.3, 3, length.out = 12)
  obj <- exptilt_object(delay, 1e-6)
  expect_lt(
    max_rel_diff(pcens_cdf(obj, q, 1), plnorm(q, 30, 1)), 1e-6
  )
})

test_that("the lognormal series needs to fit the terms it is given", {
  obj <- lnorm_object(cases[[2]], 0.1)
  expect_true(all(.pcens_tilt_fits(obj, 0.5, c(1, 1e4))))
  expect_true(all(.pcens_tilt_fits(obj, -0.5, c(1, 1e8))))
  # A positive tilt needs about xi t + 9 sqrt(xi t) + 30 terms
  expect_identical(
    .pcens_tilt_fits(obj, 1, c(1.8e4, 1.9e4, 1e6)), c(TRUE, FALSE, FALSE)
  )
  expect_identical(.pcens_tilt_fits(obj, 1, c(-3, 0)), c(TRUE, TRUE))
  # Delays without such a limit always fit
  expect_true(all(.pcens_tilt_fits(
    exptilt_object(exptilt_families()[[3]], 0.2), 0.5, 1e9
  )))
})

test_that("the lognormal CDF uses the numerical method where the series is
  too long", {
  # xi q above about 1.9e4 needs too many terms. The quantiles that fit use
  # the series and the rest use the numerical method, so nothing errors.
  delay <- list(pdist = plnorm, args = list(meanlog = 6, sdlog = 0.5))
  obj <- exptilt_object(delay, -20)
  q <- c(seq(200, 900, length.out = 6), seq(950, 1100, length.out = 6))
  expected <- exptilt_reference(q, 1, -20, exptilt_cdf(delay))
  actual <- expect_no_error(pcens_cdf(obj, q, 1))
  expect_true(all(actual >= 0 & actual <= 1))
  # The numerical method has the tolerances of stats::integrate()
  expect_lt(max(abs(actual - expected)), 1e-3)
  # The quantiles that fit keep the accuracy of the series
  expect_lt(max_rel_diff(actual[1:6], expected[1:6]), 1e-7)
  expect_no_error(dprimarycensored(
    q, plnorm, pwindow = 1, dprimary = dexpgrowth,
    primary_args = list(r = -20), meanlog = 6, sdlog = 0.5
  ))
})

test_that("the lognormal CDF is accurate for q near zero", {
  for (case in cases[c(1, 2, 4)]) {
    cdf <- exptilt_lnorm_cdf(case)
    family <- exptilt_lnorm_family(case)
    for (rho in c(-0.05, 0.5, 1)) {
      obj <- exptilt_object(family, rho)
      q <- c(1e-9, 1e-7, 1e-5, 1e-4, 1e-3)
      expected <- exptilt_reference(q, 2, rho, cdf)
      expect_lt(
        max_rel_diff(.pcens_cdf_exptilt(obj, q, 2), expected), 1e-7,
        label = exptilt_label(family, 2, rho)
      )
    }
  }
})

test_that("the lognormal CDF handles boundary values of q", {
  obj <- lnorm_object(cases[[1]], 0.3)
  result <- .pcens_cdf_exptilt(obj, c(-Inf, Inf, NA), 2)
  expect_identical(result, c(0, 1, NA_real_))
  expect_identical(.pcens_cdf_exptilt(obj, numeric(0), 2), numeric(0))
  expect_equal(
    .pcens_cdf_exptilt(obj, c(1e4, 1e6), 2), c(1, 1),
    tolerance = 1e-12
  )
  expect_identical(.pcens_cdf_exptilt(obj, c(-3, -1e-12, 0), 2), c(0, 0, 0))
  # Endpoints shared between q and q - pwindow
  endpoints <- .exptilt_endpoints(1:10, 3, lower = 0)
  expect_identical(endpoints, c(0, 1:10))
})

test_that("the lognormal CDF is a CDF for either sign of the tilt", {
  for (case in cases) {
    for (rho in c(-0.5, 1e-8, 0.3, 1)) {
      obj <- lnorm_object(case, rho)
      q <- seq(-2, 60, by = 0.25)
      result <- .pcens_cdf_exptilt(obj, q, 3)
      expect_true(all(result >= 0 & result <= 1))
      # Rounding can make the upper tail decrease by about 1e-14
      expect_gte(min(diff(result)), -1e-12)
    }
  }
})

test_that("large tilts do not overflow for the lognormal", {
  for (case in cases[1:3]) {
    for (rho in c(-5, 20, 100)) {
      obj <- lnorm_object(case, rho)
      q <- c(0.5, 2, 5, 30)
      result <- .pcens_cdf_exptilt(obj, q, 3)
      expect_true(all(is.finite(result)))
      expect_true(all(result >= 0 & result <= 1))
    }
  }
})

test_that("the lognormal pmf matches differences of the reference CDF", {
  for (case in cases[1:4]) {
    cdf <- exptilt_lnorm_cdf(case)
    for (rho in c(-0.5, 1e-8, 0.4)) {
      obj <- lnorm_object(case, rho)
      ref_cdf <- exptilt_reference(0:12, 3, rho, cdf)
      expected <- diff(ref_cdf)
      actual <- pcens_pmf(obj, 0:11, pwindow = 3)
      expect_equal(actual, expected, tolerance = 1e-8)
    }
  }
})

test_that("use_numeric = TRUE uses the default method for the lognormal", {
  obj <- lnorm_object(cases[[1]], 0.5)
  expect_identical(
    pcens_cdf(obj, c(1, 4), 2, use_numeric = TRUE),
    pcens_cdf.default(obj, c(1, 4), 2)
  )
})

test_that("pprimarycensored uses the lognormal transform with truncation", {
  plnorm_obj <- function(q, ...) {
    pprimarycensored(
      q, plnorm,
      pwindow = 2, dprimary = dexpgrowth, primary_args = list(r = 0.3),
      meanlog = 1.6, sdlog = 0.5, ...
    )
  }
  ref <- function(x) {
    exptilt_reference(x, 2, 0.3, function(u) plnorm(u, 1.6, 0.5))
  }
  q <- c(0.5, 2, 4, 7, 12, 20)
  expect_equal(plnorm_obj(q, D = 25), ref(q) / ref(25), tolerance = 1e-7)
})

test_that("the lognormal transform agrees with the delays it generalises", {
  # A lognormal with a tiny sdlog is a point mass at exp(meanlog). The
  # transform is then a step of size exp(xi * exp(meanlog)) at that point
  obj <- new_pcens(
    plnorm, dexpgrowth, list(r = 0.3),
    meanlog = log(2), sdlog = 0.02
  )
  t <- c(1, 3)
  lower <- .pcens_tilt_transform(obj, t, -0.3)
  expect_equal(exp(lower[1]), 0, tolerance = 1e-12)
  expect_equal(exp(lower[2]), exp(-0.3 * 2), tolerance = 2e-3)
})

test_that("pcens_cdf uses the numerical method for a few quantiles", {
  # The quadrature has a fixed cost that the numerical method beats for
  # fewer than 10 quantiles, see `.lnorm_exptilt_min_q`
  obj <- lnorm_object(cases[[1]], 0.3)
  few <- c(0.5, 2, 4, 8)
  expect_identical(
    pcens_cdf(obj, few, 2), pcens_cdf.default(obj, few, 2)
  )
  many <- seq(0.5, 20, length.out = 12)
  expect_identical(
    pcens_cdf(obj, many, 2), .pcens_cdf_exptilt(obj, many, 2)
  )
  expect_equal(
    pcens_cdf(obj, many, 2), pcens_cdf(obj, many, 2, use_numeric = TRUE),
    tolerance = 1e-6
  )
  # The other delays use the closed forms for any number of quantiles
  gamma_obj <- exptilt_object(exptilt_families()[[3]], 0.3)
  expect_identical(
    pcens_cdf(gamma_obj, few, 2), .pcens_cdf_exptilt(gamma_obj, few, 2)
  )
})

test_that("a missing tilt is an error for the lognormal with few quantiles", {
  obj <- new_pcens(plnorm, dexpgrowth, list(), meanlog = 1, sdlog = 0.5)
  expect_error(
    pcens_cdf(obj, 1, 2), "r parameter is required for the exponential growth"
  )
})

test_that("the lognormal series stops where it needs too many terms", {
  expect_error(
    .lnorm_tilt_series(1e6, 0, 1, 1), "needs more than 20000 terms"
  )
  # A large tilt times the point still gives a finite transform
  lower <- .lnorm_tilt_series(c(10, 2000), 0, 1, 1)
  expect_true(all(is.finite(lower)))
  expect_lt(lower[1], lower[2])
})

test_that("default parameters of the lognormal are as in plnorm", {
  obj <- new_pcens(plnorm, dexpgrowth, list(r = 0.2))
  q <- seq(0.3, 12, length.out = 14)
  expect_equal(
    pcens_cdf(obj, q, 1), pcens_cdf(obj, q, 1, use_numeric = TRUE),
    tolerance = 1e-6
  )
})
