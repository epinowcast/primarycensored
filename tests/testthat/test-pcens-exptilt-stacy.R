skip_if_not_installed("flexsurv")

# Weibull and generalised gamma delays with an exponentially tilted primary.

stacy_families <- exptilt_stacy_families()
pwindows <- c(0.5, 2, 7)
rhos <- c(-1, -0.3, -1e-3, -1e-5, 1e-5, 1e-3, 0.05, 0.3, 1)

test_that("the series transform matches quadrature of the density", {
  ts <- c(1e-3, 0.4, 3.5, 12, 40)
  for (family in stacy_families) {
    obj <- exptilt_object(family, 0.2)
    for (xi in c(-1, -0.25, -0.05, 0.05, 0.25)) {
      actual <- .pcens_tilt_transform(obj, ts, xi)
      keep <- !is.na(actual)
      expect_true(any(keep))
      expect_equal(
        c(actual)[keep], log(stacy_transform_reference(family, ts, xi))[keep],
        tolerance = 1e-8, info = paste(family$label, xi)
      )
    }
  }
})

test_that("the series is NaN where its terms cancel or do not converge", {
  obj <- exptilt_object(stacy_families[[2]], 1)
  expect_false(anyNA(.pcens_tilt_transform(obj, c(0.4, 3), -1)))
  expect_true(is.na(.pcens_tilt_transform(obj, 30, -1)))
  expect_true(all(is.na(.pcens_tilt_transform(obj, c(1000, 2000), 0.5))))
  # The default integration has a relative tolerance of about 1e-4
  expect_equal(
    pcens_cdf(obj, c(0.5, 30, 2000), 2),
    pcens_cdf(obj, c(0.5, 30, 2000), 2, use_numeric = TRUE),
    tolerance = 1e-4
  )
})

test_that("the xi = 0 transform is the delay CDF and its survival function", {
  ts <- c(-1, 0, 1e-3, 0.4, 12, 300)
  for (family in stacy_families) {
    obj <- exptilt_object(family, 0.2)
    expect_equal(
      c(.pcens_tilt_transform(obj, ts, 0)), log(exptilt_cdf(family)(ts)),
      tolerance = 1e-13
    )
    # The survival function underflows before its log does
    expected <- do.call(
      family$pdist,
      c(list(ts), family$args, list(lower.tail = FALSE, log.p = TRUE))
    )
    upper <- c(.pcens_tilt_transform(obj, ts, 0, upper = TRUE))
    expect_equal(
      upper[is.finite(expected)], expected[is.finite(expected)],
      tolerance = 1e-9
    )
    lower <- c(.pcens_tilt_transform(obj, c(-2, 0), -0.3))
    expect_identical(lower, c(-Inf, -Inf))
    expect_true(all(is.na(.pcens_tilt_transform(obj, 1, -0.3, upper = TRUE))))
  }
})

test_that("the partial moments match quadrature", {
  ts <- c(1e-2, 0.3, 4, 20)
  for (family in stacy_families) {
    obj <- exptilt_object(family, 0.2)
    moments <- .pcens_tilt_moments(obj, c(-1, 0, ts))
    expect_identical(unname(moments[1:2, ]), matrix(-Inf, 2, 3))
    for (j in 1:3) {
      expected <- vapply(ts, function(tt) {
        stats::integrate(
          function(u) {
            (tt - u)^j * do.call(family$ddist, c(list(u), family$args))
          },
          0, tt,
          rel.tol = 1e-12, abs.tol = 0
        )$value
      }, numeric(1))
      expect_equal(
        unname(moments[-(1:2), j]), log(expected),
        tolerance = 1e-8, info = family$label
      )
    }
  }
})

test_that("the analytic CDF matches a reference integral", {
  # The small tilt forms have a truncation error of up to about 3e-8
  for (family in stacy_families) {
    for (pwindow in pwindows) {
      q <- c(
        1e-6, 1e-3, 0.3 * pwindow, pwindow - 1e-3, pwindow, pwindow + 1e-3,
        3, 6, 12, 25
      )
      for (rho in rhos) {
        expected <- exptilt_reference(q, pwindow, rho, exptilt_cdf(family))
        actual <- pcens_cdf(exptilt_object(family, rho), q, pwindow)
        expect_lt(
          max_rel_diff(actual, expected), 1e-7,
          label = exptilt_label(family, pwindow, rho)
        )
      }
    }
  }
})

test_that("the analytic CDF agrees with use_numeric = TRUE", {
  # The default integration has a relative tolerance of about 1e-4
  q <- c(0.05, 0.5, 1.5, 3, 6, 12, 20)
  for (family in stacy_families) {
    for (pwindow in c(1, 7)) {
      for (rho in c(-1, -0.5, -1e-8, 1e-8, 0.5, 1)) {
        obj <- exptilt_object(family, rho)
        expect_equal(
          pcens_cdf(obj, q, pwindow),
          pcens_cdf(obj, q, pwindow, use_numeric = TRUE),
          tolerance = 1e-4, info = exptilt_label(family, pwindow, rho)
        )
      }
    }
  }
})

test_that("the analytic CDF is continuous through zero tilt and the small
  tilt threshold", {
  for (family in stacy_families) {
    q <- c(1e-3, 0.6, 2, 6, 12)
    uniform <- exptilt_reference(q, 2, 0, exptilt_cdf(family))
    for (rho in c(0, 1e-12, -1e-9, 1e-7)) {
      expect_lt(
        max_rel_diff(pcens_cdf(exptilt_object(family, rho), q, 2), uniform),
        1e-6
      )
    }
    for (sign in c(-1, 1)) {
      expect_lt(
        max_rel_diff(
          pcens_cdf(exptilt_object(family, sign * 0.9999e-4 / 2), q, 2),
          pcens_cdf(exptilt_object(family, sign * 1.0001e-4 / 2), q, 2)
        ),
        1e-7
      )
    }
  }
})

test_that("the analytic CDF is accurate for q near zero and is a CDF", {
  for (family in stacy_families) {
    for (rho in c(-0.05, 0.5, 1)) {
      obj <- exptilt_object(family, rho)
      q <- c(1e-9, 1e-7, 1e-5, 1e-3)
      expect_lt(
        max_rel_diff(
          pcens_cdf(obj, q, 2),
          exptilt_reference(q, 2, rho, exptilt_cdf(family))
        ),
        1e-7
      )
      grid <- pcens_cdf(obj, seq(-2, 60, by = 0.25), 3)
      expect_true(all(grid >= 0 & grid <= 1))
      expect_gte(min(diff(grid)), -1e-6)
    }
    obj <- exptilt_object(family, 0.3)
    expect_identical(pcens_cdf(obj, c(-Inf, Inf, NA), 2), c(0, 1, NA_real_))
    expect_identical(pcens_cdf(obj, c(-3, 0), 2), c(0, 0))
  }
})

test_that("the PMF matches differences of the reference CDF", {
  for (family in stacy_families) {
    for (rho in c(-0.5, 1e-8, 0.4)) {
      expect_equal(
        pcens_pmf(exptilt_object(family, rho), 0:11, pwindow = 3),
        diff(exptilt_reference(0:12, 3, rho, exptilt_cdf(family))),
        tolerance = 1e-7
      )
    }
  }
})

test_that("the upper tail PMF matches a survival based reference", {
  x <- 0:70
  for (family in stacy_families) {
    for (rho in c(0.3, 0.05, -0.3)) {
      for (pwindow in c(1, 3)) {
        expected <- exptilt_pmf_reference(family, x, pwindow, rho)
        actual <- pcens_pmf(exptilt_object(family, rho), x, pwindow = pwindow)
        keep <- expected > 1e-8
        expect_lt(
          max_rel_diff(actual[keep], expected[keep]), 1e-6,
          label = exptilt_label(family, pwindow, rho)
        )
      }
    }
  }
})

test_that("a short tailed delay uses the direct form in the upper tail", {
  family <- stacy_families[[6]]
  obj <- exptilt_object(family, 2)
  q <- c(0.5, 1.5, 2.5, 3)
  expected <- exptilt_reference(q, 5, 2, exptilt_cdf(family))
  expect_lt(max_rel_diff(pcens_cdf(obj, q, 5), expected), 1e-7)
})

test_that("the transform is evaluated once per unique endpoint", {
  obj <- exptilt_object(stacy_families[[1]], 0.3)
  calls <- new.env()
  calls$n <- numeric(0)
  testthat::local_mocked_bindings(
    .pcens_tilt_transform = function(object, t, xi, upper = FALSE) {
      if (xi == -0.3 && !upper) {
        calls$n <- c(calls$n, length(t))
      }
      .stacy_tilt_transform(object, t, xi, upper)
    }
  )
  pcens_cdf(obj, 1:10, 3)
  # The endpoints are 0, which stands for every point at or below 0, and 1:10
  expect_identical(calls$n, 11)
})

test_that("the Prentice parameterisation maps to the Stacy one", {
  mu <- 0.4
  sigma <- 0.7
  Q <- 1.3
  prentice <- new_pcens(
    flexsurv::pgengamma, dexpgrowth, list(r = 0.3),
    mu = mu, sigma = sigma, Q = Q
  )
  stacy <- new_pcens(
    flexsurv::pgengamma.orig, dexpgrowth, list(r = 0.3),
    shape = Q / sigma, scale = exp(mu) * Q^(2 * sigma / Q), k = Q^-2
  )
  q <- c(0.2, 1, 3, 8, 20)
  expect_identical(pcens_cdf(prentice, q, 2), pcens_cdf(stacy, q, 2))
  expect_equal(
    pcens_cdf(prentice, q, 2), pcens_cdf(prentice, q, 2, use_numeric = TRUE),
    tolerance = 1e-6
  )
  # Q <= 0 is not in the Stacy family and uses the numerical method
  for (Q in c(0, -0.5)) {
    obj <- new_pcens(
      flexsurv::pgengamma, dexpgrowth, list(r = 0.3),
      mu = mu, sigma = sigma, Q = Q
    )
    expect_false(.pcens_tilt_available(obj, -0.3))
    expect_identical(pcens_cdf(obj, q, 2), pcens_cdf.default(obj, q, 2))
  }
})

test_that("the analytic CDF matches the empirical CDF of rprimarycensored", {
  set.seed(2024)
  n <- 1e5
  pwindow <- 2
  q <- c(1, 2.5, 5, 8, 14)
  cases <- list(
    stacy_families[[1]], stacy_families[[4]],
    list(
      pdist = flexsurv::pgengamma, rdist = flexsurv::rgengamma,
      args = list(mu = 1.2, sigma = 0.6, Q = 1.5)
    )
  )
  for (case in cases) {
    for (rho in c(-0.5, 0.3)) {
      obj <- do.call(
        new_pcens,
        c(
          list(
            pdist = case$pdist, dprimary = dexpgrowth,
            primary_args = list(r = rho)
          ),
          case$args
        )
      )
      samples <- do.call(
        rprimarycensored,
        c(
          list(
            n, case$rdist, pwindow = pwindow, swindow = 0,
            rprimary = rexpgrowth, rprimary_args = list(r = rho)
          ),
          case$args
        )
      )
      expected <- pcens_cdf(obj, q, pwindow)
      empirical <- vapply(q, function(qq) mean(samples <= qq), numeric(1))
      body <- expected > 0.01 & expected < 0.99
      z <- (empirical - expected)[body] /
        sqrt(expected * (1 - expected) / n)[body]
      expect_lt(max(abs(z)), 4.5, label = paste(class(obj)[1], "r =", rho))
    }
  }
})
