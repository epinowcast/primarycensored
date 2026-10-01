families <- exptilt_families()
pwindows <- c(0.5, 1, 2, 7)
# The tilts in the acceptance criteria plus the small tilt regimes
rhos <- c(
  -1, -0.5, -0.05, -3e-4, -1e-4, -1e-5, -1e-8, 1e-8, 1e-5, 1e-4, 3e-4, 0.05,
  0.5, 1
)

test_that("exponentially tilted primaries dispatch to analytic methods", {
  for (family in families) {
    obj <- exptilt_object(family, 0.2)
    expect_s3_class(obj, "pcens")
  }
  expect_s3_class(
    exptilt_object(families[[1]], 0.2), "pcens_pexp_dexpgrowth"
  )
  expect_s3_class(
    exptilt_object(families[[3]], 0.2), "pcens_pgamma_dexpgrowth"
  )
  expect_s3_class(
    exptilt_object(families[[6]], 0.2), "pcens_pnorm_dexpgrowth"
  )
  for (cls in c(
    "pcens_pexp_dexpgrowth", "pcens_pgamma_dexpgrowth",
    "pcens_pnorm_dexpgrowth"
  )) {
    expect_false(
      is.null(utils::getS3method("pcens_cdf", cls, optional = TRUE)),
      info = cls
    )
  }
})

test_that("tilt transforms match numerical integration of the density", {
  ts <- c(1e-3, 0.4, 1, 3.5, 12, 40)
  xis <- c(0, -1, -0.25, 0.1, 0.25)
  for (family in families) {
    obj <- exptilt_object(family, 0.2)
    ddist <- switch(family$label,
      "exponential rate 2" = function(u) dexp(u, 2),
      "exponential rate 0.3" = function(u) dexp(u, 0.3),
      "gamma shape 0.6" = function(u) dgamma(u, 0.6, 1.3),
      "gamma shape 2.5" = function(u) dgamma(u, 2.5, scale = 2.5),
      "gamma shape 20" = function(u) dgamma(u, 20, 4),
      "normal mean 3" = function(u) dnorm(u, 3, 2),
      "normal mean -1" = function(u) dnorm(u, -1, 3)
    )
    # The integrand of the normal is not finite at -Inf
    lower_end <- if (family$positive) 0 else -60
    for (xi in xis) {
      if (family$positive && family$rate - xi <= 0) {
        expect_false(.pcens_tilt_available(obj, xi), info = family$label)
        next
      }
      expect_true(.pcens_tilt_available(obj, xi), info = family$label)
      lower <- .pcens_tilt_transform(obj, ts, xi)
      upper <- .pcens_tilt_transform(obj, ts, xi, upper = TRUE)
      expected <- vapply(ts, function(t) {
        stats::integrate(
          function(u) exp(xi * u) * ddist(u), lower_end, t,
          rel.tol = 1e-12, abs.tol = 0
        )$value
      }, numeric(1))
      expect_equal(exp(lower), expected, tolerance = 1e-9, info = family$label)
      total <- exp(lower) + exp(upper)
      # The transform over the whole support is the same at every t
      expect_equal(total, rep(total[1], length(ts)), tolerance = 1e-12)
    }
  }
})

test_that("tilt transforms are 0 or total below the support of positive
  delays", {
  obj <- exptilt_object(families[[3]], 0.2)
  expect_identical(.pcens_tilt_transform(obj, c(-2, 0), -0.3), c(-Inf, -Inf))
  total <- (1.3 / (1.3 + 0.3))^0.6
  expect_equal(
    exp(.pcens_tilt_transform(obj, c(-2, 0), -0.3, upper = TRUE)),
    rep(total, 2),
    tolerance = 1e-14
  )
})

test_that("the analytic CDF matches a reference integral", {
  for (family in families) {
    for (pwindow in pwindows) {
      q <- sort(c(
        1e-6, 1e-3, 0.3 * pwindow, pwindow - 1e-3, pwindow, pwindow + 1e-3,
        2, 3, 6, 12, 25,
        if (!family$positive) c(-10, -3, -0.5)
      ))
      for (rho in rhos) {
        if (family$positive && !exptilt_admissible(family, rho) &&
          abs(rho) * pwindow >= 1e-2) {
          next
        }
        obj <- exptilt_object(family, rho)
        expected <- exptilt_reference(q, pwindow, rho, exptilt_cdf(family))
        actual <- pcens_cdf(obj, q, pwindow)
        expect_lt(
          max_rel_diff(actual, expected), 1e-7,
          label = exptilt_label(family, pwindow, rho)
        )
      }
    }
  }
})

test_that("the analytic CDF agrees with use_numeric = TRUE", {
  q <- c(0.05, 0.5, 1.5, 3, 6, 12, 20)
  for (family in families) {
    for (pwindow in c(1, 2, 7)) {
      for (rho in c(-1, -0.5, -1e-8, 1e-8, 0.5, 1)) {
        obj <- exptilt_object(family, rho)
        analytic <- pcens_cdf(obj, q, pwindow)
        numeric <- pcens_cdf(obj, q, pwindow, use_numeric = TRUE)
        expect_equal(
          analytic, numeric,
          tolerance = 1e-6,
          info = sprintf(
            "%s, pwindow = %g, r = %g", family$label, pwindow, rho
          )
        )
      }
    }
  }
})

test_that("use_numeric = TRUE uses the default method", {
  obj <- exptilt_object(families[[3]], 0.5)
  expect_identical(
    pcens_cdf(obj, c(1, 4), 2, use_numeric = TRUE),
    pcens_cdf.default(obj, c(1, 4), 2)
  )
})

test_that("delay families without a tilt form use the numerical method", {
  obj <- new_pcens(
    pdist = pweibull, dprimary = dexpgrowth,
    primary_args = list(r = 0.2), shape = 2, scale = 3
  )
  expect_false(.pcens_tilt_available(obj, -0.2))
  q <- c(0.5, 2, 5)
  expect_identical(pcens_cdf(obj, q, 2), pcens_cdf.default(obj, q, 2))
})

test_that("inadmissible tilts use the numerical method", {
  for (family in families[c(2, 3, 4)]) {
    for (rho in c(-1, -0.5)) {
      if (exptilt_admissible(family, rho)) {
        next
      }
      obj <- exptilt_object(family, rho)
      expect_false(.pcens_tilt_available(obj, -rho))
      q <- c(0.5, 2, 5, 10)
      expect_identical(
        pcens_cdf(obj, q, 2), pcens_cdf.default(obj, q, 2),
        info = sprintf("%s, r = %g", family$label, rho)
      )
    }
  }
})

test_that("the analytic CDF is continuous in the tilt through zero", {
  for (family in families) {
    for (pwindow in c(0.5, 2, 7)) {
      q <- c(1e-3, 0.3 * pwindow, pwindow, 3, 6, 12)
      obj_zero <- exptilt_object(family, 0)
      # A tilt of exactly zero is the uniform window
      uniform <- exptilt_reference(q, pwindow, 0, exptilt_cdf(family))
      expect_lt(max_rel_diff(pcens_cdf(obj_zero, q, pwindow), uniform), 1e-9)
      for (sign in c(-1, 1)) {
        for (rho in sign * c(1e-12, 1e-9, 1e-7)) {
          obj <- exptilt_object(family, rho)
          expect_lt(
            max_rel_diff(pcens_cdf(obj, q, pwindow), uniform),
            1e-6,
            label = exptilt_label(family, pwindow, rho)
          )
        }
      }
    }
  }
})

test_that("the analytic CDF has no jump where the small tilt form ends", {
  for (family in families) {
    limit <- if (family$positive) 1e-2 else 1e-3
    for (pwindow in c(0.5, 2, 7)) {
      q <- c(1e-3, 0.3 * pwindow, pwindow, 3, 6, 12, 25)
      for (sign in c(-1, 1)) {
        below <- exptilt_object(family, sign * 0.999999 * limit / pwindow)
        above <- exptilt_object(family, sign * 1.000001 * limit / pwindow)
        expect_lt(
          max_rel_diff(
            pcens_cdf(below, q, pwindow), pcens_cdf(above, q, pwindow)
          ),
          1e-7,
          label = exptilt_label(family, pwindow, sign * limit / pwindow)
        )
      }
    }
  }
})

test_that("the analytic CDF is accurate for gamma delays with large shapes", {
  # Mean 10, so that the lower tail and the bulk are at q = 6 to 10
  q <- c(6, 8, 9.5, 10, 10.5)
  for (shape in c(500, 1000, 5000)) {
    for (pwindow in c(1, 7)) {
      for (rho in c(-1e-3, -1e-4, -1e-5, 2e-5, 1.5e-4, 1e-3, 1e-2)) {
        obj <- new_pcens(
          pdist = pgamma, dprimary = dexpgrowth, primary_args = list(r = rho),
          shape = shape, rate = shape / 10
        )
        expected <- exptilt_gamma_log_reference(
          q, pwindow, rho, shape, shape / 10
        )
        expect_lt(
          max(abs(expm1(log(pcens_cdf(obj, q, pwindow)) - expected))), 1e-6,
          label = sprintf(
            "shape = %g, pwindow = %g, r = %g", shape, pwindow, rho
          )
        )
      }
    }
  }
})

test_that("the analytic CDF is accurate in the lower tail of a normal delay
  for small tilts", {
  grid <- exptilt_normal_tail_grid()
  expected <- exptilt_normal_tail_reference(grid)
  actual <- vapply(seq_len(nrow(grid)), function(i) {
    obj <- new_pcens(
      pdist = pnorm, dprimary = dexpgrowth,
      primary_args = list(r = grid$rho[i]), mean = -4, sd = 0.3
    )
    log(pcens_cdf(obj, grid$d[i], grid$pwindow[i]))
  }, numeric(1))
  expect_lt(max(abs(expm1(actual - expected))), 3e-7)
})

test_that("the tilt moments are the moments of the delay about t", {
  t <- c(0.05, 0.7, 2, 6, 15)
  for (family in families) {
    obj <- exptilt_object(family, 0.1)
    moments <- .pcens_tilt_moments(obj, t)
    expect_identical(colnames(moments), c("G1", "G2", "G3"))
    ddist <- if (identical(family$pdist, pexp)) {
      dexp
    } else if (identical(family$pdist, pgamma)) {
      dgamma
    } else {
      dnorm
    }
    density <- function(x) do.call(ddist, c(list(x), family$args))
    lower <- if (family$positive) 0 else -Inf
    expected <- vapply(1:3, function(k) {
      vapply(t, function(tt) {
        stats::integrate(
          function(u) (tt - u)^k * density(u), lower, tt,
          rel.tol = 1e-12, abs.tol = 0
        )$value
      }, numeric(1))
    }, numeric(length(t)))
    expect_equal(
      exp(moments), expected, tolerance = 1e-8, ignore_attr = TRUE,
      info = family$label
    )
  }
})

test_that("the analytic CDF is accurate for q near zero", {
  for (family in families[c(2, 3, 4, 5)]) {
    cdf <- exptilt_cdf(family)
    for (rho in c(-0.05, 0.5, 1)) {
      obj <- exptilt_object(family, rho)
      q <- c(1e-9, 1e-7, 1e-5, 1e-4, 1e-3)
      expected <- exptilt_reference(q, 2, rho, cdf)
      expect_lt(
        max_rel_diff(pcens_cdf(obj, q, 2), expected), 1e-7,
        label = sprintf("%s, r = %g", family$label, rho)
      )
    }
  }
})

test_that("the analytic CDF handles boundary values of q", {
  for (family in families) {
    obj <- exptilt_object(family, 0.3)
    result <- pcens_cdf(obj, c(-Inf, Inf, NA), 2)
    expect_identical(result, c(0, 1, NA_real_))
    expect_identical(pcens_cdf(obj, numeric(0), 2), numeric(0))
    result <- pcens_cdf(obj, c(1e4, 1e6), 2)
    expect_equal(result, c(1, 1), tolerance = 1e-12)
    if (family$positive) {
      expect_identical(pcens_cdf(obj, c(-3, -1e-12, 0), 2), c(0, 0, 0))
    } else {
      expect_equal(
        pcens_cdf(obj, -50, 2), 0,
        tolerance = 1e-12
      )
    }
  }
})

test_that("the analytic CDF is a CDF", {
  for (family in families) {
    for (rho in c(-0.5, 1e-8, 0.3, 1)) {
      if (family$positive && !exptilt_admissible(family, rho)) {
        next
      }
      obj <- exptilt_object(family, rho)
      q <- seq(-2, 60, by = 0.25)
      result <- pcens_cdf(obj, q, 3)
      expect_true(all(result >= 0 & result <= 1))
      # Rounding can make the upper tail decrease by about 1e-14
      expect_gte(min(diff(result)), -1e-12)
    }
  }
})

test_that("large tilts do not overflow", {
  for (family in families) {
    for (rho in c(-20, 20, 100)) {
      if (family$positive && !exptilt_admissible(family, rho)) {
        next
      }
      obj <- exptilt_object(family, rho)
      q <- c(0.5, 2, 5, 30)
      result <- pcens_cdf(obj, q, 3)
      expect_true(all(is.finite(result)))
      expect_true(all(result >= 0 & result <= 1))
    }
  }
})

test_that("the pmf matches differences of the reference CDF", {
  for (family in families) {
    cdf <- exptilt_cdf(family)
    for (rho in c(-0.5, 1e-8, 0.4)) {
      if (family$positive && !exptilt_admissible(family, rho)) {
        next
      }
      obj <- exptilt_object(family, rho)
      ref_cdf <- exptilt_reference(0:12, 3, rho, cdf)
      expected <- diff(ref_cdf)
      actual <- pcens_pmf(obj, 0:11, pwindow = 3)
      expect_equal(actual, expected, tolerance = 1e-8)
    }
  }
})

test_that("pprimarycensored uses the analytic methods with truncation", {
  pnorm_obj <- function(q, ...) {
    pprimarycensored(
      q, pnorm,
      pwindow = 2, dprimary = dexpgrowth, primary_args = list(r = 0.3),
      mean = 3, sd = 2, ...
    )
  }
  cdf <- function(x) pnorm(x, 3, 2)
  ref <- function(x) exptilt_reference(x, 2, 0.3, cdf)
  q <- c(-1, 0.5, 2, 4, 7)
  expected <- (ref(q) - ref(-2)) / (ref(9) - ref(-2))
  expect_equal(pnorm_obj(q, L = -2, D = 9), expected, tolerance = 1e-8)
  expect_equal(pnorm_obj(q, D = 9), ref(q) / ref(9), tolerance = 1e-8)
})

test_that("a missing or invalid tilt is an error", {
  obj <- new_pcens(pgamma, dexpgrowth, list(), shape = 2, rate = 1)
  expect_error(
    pcens_cdf(obj, 1, 2), "r parameter is required for the exponential growth"
  )
  obj <- new_pcens(
    pgamma, dexpgrowth, list(r = c(0.1, 0.2)),
    shape = 2, rate = 1
  )
  expect_error(pcens_cdf(obj, 1, 2), "single finite")
  obj <- new_pcens(
    pexp, dexpgrowth, list(r = NA_real_),
    rate = 1
  )
  expect_error(pcens_cdf(obj, 1, 2), "single finite")
})

test_that("missing delay parameters are errors", {
  obj <- new_pcens(pgamma, dexpgrowth, list(r = 0.2), rate = 1)
  expect_error(
    pcens_cdf(obj, 1, 2), "shape parameter is required for Gamma"
  )
})

test_that("delay parameters default as in the stats functions", {
  obj <- new_pcens(pexp, dexpgrowth, list(r = 0.2))
  expect_equal(
    pcens_cdf(obj, c(0.5, 2), 1),
    pcens_cdf(obj, c(0.5, 2), 1, use_numeric = TRUE),
    tolerance = 1e-6
  )
  # A gamma with only a shape has rate 1, as in pgamma()
  obj <- new_pcens(pgamma, dexpgrowth, list(r = 0.3), shape = 2)
  expect_equal(
    pcens_cdf(obj, c(1, 5), 2),
    c(0.04114682, 0.89135648),
    tolerance = 1e-7
  )
  expect_equal(
    pcens_cdf(obj, c(1, 5), 2),
    pcens_cdf(obj, c(1, 5), 2, use_numeric = TRUE),
    tolerance = 1e-6
  )
  expect_equal(
    pprimarycensored(
      c(1, 5), pgamma,
      pwindow = 2, dprimary = dexpgrowth, primary_args = list(r = 0.3),
      shape = 2
    ),
    c(0.04114682, 0.89135648),
    tolerance = 1e-7
  )
  obj_rate <- new_pcens(
    pgamma, dexpgrowth, list(r = 0.3),
    shape = 2, rate = 1
  )
  expect_identical(
    pcens_cdf(obj, c(1, 5), 2), pcens_cdf(obj_rate, c(1, 5), 2)
  )
  # A rate that would not admit the tilt is checked against the default
  obj <- new_pcens(pgamma, dexpgrowth, list(r = -1.5), shape = 2)
  expect_equal(
    pcens_cdf(obj, c(1, 5), 2),
    pcens_cdf(obj, c(1, 5), 2, use_numeric = TRUE),
    tolerance = 1e-6
  )
  obj <- new_pcens(pnorm, dexpgrowth, list(r = 0.2))
  expect_equal(
    pcens_cdf(obj, c(-0.5, 2), 1),
    pcens_cdf(obj, c(-0.5, 2), 1, use_numeric = TRUE),
    tolerance = 1e-6
  )
})

test_that("endpoints are shared between neighbouring delays", {
  endpoints <- .exptilt_endpoints(1:10, 3, lower = 0)
  expect_identical(endpoints, c(0, 1:10))
  endpoints <- .exptilt_endpoints(1:10, 3, lower = -Inf)
  expect_identical(endpoints, as.numeric(-2:10))
  endpoints <- .exptilt_endpoints(c(0.5, 2.5, 4), 2, lower = 0)
  expect_identical(endpoints, c(0, 0.5, 2, 2.5, 4))
})

test_that("the analytic CDF matches rprimarycensored samples", {
  set.seed(101)
  n <- 20000
  pwindow <- 2
  for (family in families[c(1, 3, 4, 6, 7)]) {
    for (rho in c(-0.2, 0.5)) {
      obj <- exptilt_object(family, rho)
      samples <- do.call(
        rprimarycensored,
        c(
          list(
            n = n, rdist = family$rdist, pwindow = pwindow, swindow = 0,
            rprimary = rexpgrowth, rprimary_args = list(r = rho)
          ),
          family$args
        )
      )
      ks <- stats::ks.test(samples, function(x) pcens_cdf(obj, x, pwindow))
      expect_gt(ks$p.value, 1e-3, label = exptilt_label(family, pwindow, rho))
    }
  }
})

test_that("pcens_cdf is unchanged for the uniform primary", {
  obj <- new_pcens(pgamma, dunif, list(), shape = 2, rate = 0.5)
  expect_s3_class(obj, "pcens_pgamma_dunif")
  expect_equal(
    pcens_cdf(obj, c(1, 4), 2),
    pcens_cdf(obj, c(1, 4), 2, use_numeric = TRUE),
    tolerance = 1e-6
  )
})
