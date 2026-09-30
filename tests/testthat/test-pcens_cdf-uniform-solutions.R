# Uniform primary analytical solutions for the exponential, normal,
# chi-square and beta delays. Each is checked against a near double
# precision quadrature of the delay CDF over the primary window (see
# helper-uniform-reference.R) to a relative 1e-9, and against the package
# numerical path with `use_numeric = TRUE` to 1e-6 on a grid where its
# default stats::integrate tolerance is enough. The grids include windows
# wider than q, q near 0, q below 0 for the normal, and both tails.

delays_positive <- c(0.001, 0.01, 0.1, 0.5, 1, 2, 5, 10, 30, 100)

test_that("the new delays have analytical pcens_cdf methods", {
  classes <- c(
    "pcens_pexp_dunif", "pcens_pnorm_dunif", "pcens_pchisq_dunif",
    "pcens_pbeta_dunif"
  )
  for (class in classes) {
    expect_false(
      is.null(utils::getS3method("pcens_cdf", class, optional = TRUE)),
      info = class
    )
  }
  expect_s3_class(new_pcens(pexp, dunif, list()), classes[[1]])
  expect_s3_class(new_pcens(pnorm, dunif, list()), classes[[2]])
  expect_s3_class(new_pcens(pchisq, dunif, list(), df = 1), classes[[3]])
  expect_s3_class(
    new_pcens(pbeta, dunif, list(), shape1 = 1, shape2 = 1), classes[[4]]
  )
})

test_that("pcens_cdf for the exponential matches quadrature and the numerical
  path", {
  for (rate in c(0.001, 0.05, 0.5, 3, 20)) {
    obj <- new_pcens(pexp, dunif, list(), rate = rate)
    for (pwindow in c(0.5, 1, 3, 10)) {
      info <- sprintf("rate = %g, pwindow = %g", rate, pwindow)
      analytic <- pcens_cdf(obj, delays_positive, pwindow)
      reference <- reference_uniform_cdf(
        function(t) pexp(t, rate), delays_positive, pwindow, 0
      )
      expect_rel_equal(analytic, reference, info = info)
    }
  }
  obj <- new_pcens(pexp, dunif, list(), rate = 0.5)
  q <- seq(0, 30, by = 1)
  for (pwindow in c(1, 2, 5, 10)) {
    expect_equal(
      pcens_cdf(obj, q, pwindow),
      pcens_cdf(obj, q, pwindow, use_numeric = TRUE),
      tolerance = 1e-6, info = paste("pwindow =", pwindow)
    )
  }
})

test_that("pcens_cdf for the exponential agrees with the gamma with shape 1", {
  q <- c(0.001, 0.5, 2, 10, 40)
  exp_obj <- new_pcens(pexp, dunif, list(), rate = 0.3)
  gamma_obj <- new_pcens(pgamma, dunif, list(), shape = 1, rate = 0.3)
  for (pwindow in c(0.5, 2, 7)) {
    expect_equal(
      pcens_cdf(exp_obj, q, pwindow), pcens_cdf(gamma_obj, q, pwindow),
      tolerance = 1e-10
    )
  }
})

test_that("pcens_cdf for the exponential uses pexp's default rate", {
  q <- c(0.5, 2, 10)
  expect_identical(
    pcens_cdf(new_pcens(pexp, dunif, list()), q, 2),
    pcens_cdf(new_pcens(pexp, dunif, list(), rate = 1), q, 2)
  )
})

test_that("pcens_cdf for the normal matches quadrature and the numerical
  path, including negative q and both tails", {
  delays <- c(-40, -20, -10, -5, -2, -1, -0.5, 0, 0.5, 1, 2, 5, 10, 20, 40)
  for (mu in c(-3, 0, 2, 10)) {
    for (sigma in c(0.1, 1, 3)) {
      obj <- new_pcens(pnorm, dunif, list(), mean = mu, sd = sigma)
      for (pwindow in c(0.5, 1, 3, 10)) {
        info <- sprintf(
          "mean = %g, sd = %g, pwindow = %g", mu, sigma, pwindow
        )
        analytic <- pcens_cdf(obj, delays, pwindow)
        reference <- reference_uniform_cdf(
          function(t) pnorm(t, mu, sigma), delays, pwindow
        )
        expect_rel_equal(analytic, reference, info = info)
      }
    }
  }
  obj <- new_pcens(pnorm, dunif, list(), mean = 2, sd = 1.5)
  q <- seq(-6, 12, by = 1)
  for (pwindow in c(1, 2, 5)) {
    expect_equal(
      pcens_cdf(obj, q, pwindow),
      pcens_cdf(obj, q, pwindow, use_numeric = TRUE),
      tolerance = 1e-6, info = paste("pwindow =", pwindow)
    )
  }
})

test_that("pcens_cdf for the normal is continuous across the tail switch", {
  # The lower tail uses an asymptotic series below z = -10. The value from
  # either side of the switch must agree to rounding.
  mu <- 0
  sigma <- 1
  obj <- new_pcens(pnorm, dunif, list(), mean = mu, sd = sigma)
  pwindow <- 1
  below <- pcens_cdf(obj, -10 + 1e-9 - 1e-7, pwindow)
  above <- pcens_cdf(obj, -10 + 1e-9 + 1e-7, pwindow)
  expect_equal(below, above, tolerance = 1e-5)
  expect_rel_equal(
    .norm_shortfall(c(-10.000001, -9.999999)),
    vapply(c(-10.000001, -9.999999), function(z) {
      stats::integrate(
        pnorm, -Inf, z, rel.tol = 1.2e-14, abs.tol = 0, subdivisions = 500L
      )$value
    }, numeric(1)),
    tolerance = 1e-12
  )
})

test_that("pcens_cdf for the normal uses pnorm's defaults and handles
  infinite q", {
  q <- c(-2, 0, 3)
  expect_identical(
    pcens_cdf(new_pcens(pnorm, dunif, list()), q, 2),
    pcens_cdf(new_pcens(pnorm, dunif, list(), mean = 0, sd = 1), q, 2)
  )
  obj <- new_pcens(pnorm, dunif, list(), mean = 1, sd = 2)
  expect_identical(pcens_cdf(obj, c(-Inf, Inf), 1), c(0, 1))
})

test_that("pcens_cdf for the chi-square is the gamma with shape df / 2 and
  scale 2", {
  q <- c(0.001, 0.5, 2, 10, 40)
  for (df in c(0.5, 1, 3, 10, 40)) {
    obj <- new_pcens(pchisq, dunif, list(), df = df)
    gamma_obj <- new_pcens(pgamma, dunif, list(), shape = df / 2, scale = 2)
    for (pwindow in c(0.5, 3)) {
      info <- sprintf("df = %g, pwindow = %g", df, pwindow)
      analytic <- pcens_cdf(obj, q, pwindow)
      expect_identical(analytic, pcens_cdf(gamma_obj, q, pwindow), info = info)
      reference <- reference_uniform_cdf(
        function(t) pchisq(t, df), q, pwindow, 0
      )
      expect_rel_equal(analytic, reference, info = info)
    }
  }
  obj <- new_pcens(pchisq, dunif, list(), df = 4)
  q <- seq(0, 30, by = 1)
  expect_equal(
    pcens_cdf(obj, q, 2), pcens_cdf(obj, q, 2, use_numeric = TRUE),
    tolerance = 1e-6
  )
})

test_that("pcens_cdf for the non-central chi-square and beta falls back to the
  numerical path", {
  q <- c(0.5, 2, 5)
  chisq_obj <- new_pcens(pchisq, dunif, list(), df = 3, ncp = 1)
  expect_identical(
    pcens_cdf(chisq_obj, q, 1),
    pcens_cdf(chisq_obj, q, 1, use_numeric = TRUE)
  )
  beta_obj <- new_pcens(pbeta, dunif, list(), shape1 = 2, shape2 = 3, ncp = 1)
  expect_identical(
    pcens_cdf(beta_obj, c(0.2, 0.5), 0.3),
    pcens_cdf(beta_obj, c(0.2, 0.5), 0.3, use_numeric = TRUE)
  )
  # A zero ncp is the central distribution
  zero_obj <- new_pcens(pchisq, dunif, list(), df = 3, ncp = 0)
  expect_identical(
    pcens_cdf(zero_obj, q, 1),
    pcens_cdf(new_pcens(pchisq, dunif, list(), df = 3), q, 1)
  )
})

test_that("pcens_cdf for the beta matches quadrature and the numerical path
  for windows above and below the support", {
  delays <- c(0.001, 0.1, 0.3, 0.5, 0.9, 1, 1.5, 2, 3, 10)
  for (a in c(0.5, 1, 2, 5)) {
    for (b in c(0.5, 1, 3, 10)) {
      obj <- new_pcens(pbeta, dunif, list(), shape1 = a, shape2 = b)
      for (pwindow in c(0.3, 1, 3)) {
        info <- sprintf("a = %g, b = %g, pwindow = %g", a, b, pwindow)
        analytic <- pcens_cdf(obj, delays, pwindow)
        reference <- reference_uniform_cdf(
          function(t) pbeta(t, a, b), delays, pwindow, c(0, 1)
        )
        expect_rel_equal(analytic, reference, info = info)
      }
    }
  }
  obj <- new_pcens(pbeta, dunif, list(), shape1 = 2, shape2 = 3)
  q <- seq(0, 3, by = 0.25)
  expect_equal(
    pcens_cdf(obj, q, 0.5), pcens_cdf(obj, q, 0.5, use_numeric = TRUE),
    tolerance = 1e-6
  )
})

test_that("the new pcens_cdf methods error for missing parameters", {
  expect_error(
    pcens_cdf(new_pcens(pchisq, dunif, list()), 1, 1),
    "df parameter is required for Chi-square distribution"
  )
  expect_error(
    pcens_cdf(new_pcens(pbeta, dunif, list(), shape2 = 2), 1, 1),
    "shape1 parameter is required for Beta distribution"
  )
  expect_error(
    pcens_cdf(new_pcens(pbeta, dunif, list(), shape1 = 2), 1, 1),
    "shape2 parameter is required for Beta distribution"
  )
})

test_that("the new pcens_cdf methods return values in [0, 1] and increase in
  q", {
  objs <- list(
    list(obj = new_pcens(pexp, dunif, list(), rate = 0.4), q = 0:30),
    list(
      obj = new_pcens(pnorm, dunif, list(), mean = 3, sd = 2),
      q = seq(-10, 20, by = 0.5)
    ),
    list(obj = new_pcens(pchisq, dunif, list(), df = 5), q = 0:40),
    list(
      obj = new_pcens(pbeta, dunif, list(), shape1 = 2, shape2 = 2),
      q = seq(0, 3, by = 0.1)
    )
  )
  for (case in objs) {
    for (pwindow in c(0.5, 2)) {
      result <- pcens_cdf(case$obj, case$q, pwindow)
      expect_true(all(result >= 0 & result <= 1))
      expect_true(all(diff(result) >= -1e-12))
    }
  }
})

test_that("pprimarycensored and dprimarycensored work end to end for the new
  delays", {
  expect_identical(
    pprimarycensored(c(2, 5, 9), pexp, pwindow = 2, rate = 0.3),
    pcens_cdf(new_pcens(pexp, dunif, list(), rate = 0.3), c(2, 5, 9), 2)
  )
  pmf <- dprimarycensored(
    0:15, pnorm, pwindow = 2, swindow = 1, L = 0, D = 16, mean = 6, sd = 2
  )
  expect_equal(sum(pmf), 1, tolerance = 1e-10)
  # The window is not clipped at 0, so a normal delay has mass below 0
  lower <- pprimarycensored(-1, pnorm, pwindow = 2, mean = 2, sd = 2)
  expect_gt(lower, 0)
})
