# Analytical uniform primary CDFs against quadrature, see
# helper-uniform-solutions-reference.R

test_that("the new delays have analytical pcens_cdf methods", {
  expect_s3_class(new_pcens(pexp, dunif, list()), "pcens_pexp_dunif")
  expect_s3_class(new_pcens(pnorm, dunif, list()), "pcens_pnorm_dunif")
  expect_s3_class(
    new_pcens(pchisq, dunif, list(), df = 1), "pcens_pchisq_dunif"
  )
  expect_s3_class(
    new_pcens(pbeta, dunif, list(), shape1 = 1, shape2 = 1),
    "pcens_pbeta_dunif"
  )
  for (class in c("pexp", "pnorm", "pchisq", "pbeta")) {
    method <- paste0("pcens_", class, "_dunif")
    expect_false(
      is.null(utils::getS3method("pcens_cdf", method, optional = TRUE)),
      info = method
    )
  }
})

test_that("the new pcens_cdf methods match quadrature", {
  delays_positive <- c(0.001, 0.01, 0.1, 0.5, 1, 2, 5, 10, 30, 100)
  # Each case is a delay with a grid of its arguments, the delays and
  # windows at which to evaluate it, and the kinks of F_T in t
  cases <- list(
    list(
      pdist = pexp, grid = expand.grid(rate = c(0.001, 0.05, 0.5, 3, 20)),
      delays = delays_positive, pwindows = c(0.5, 1, 3, 10), kinks = 0
    ),
    list(
      pdist = pnorm,
      grid = expand.grid(mean = c(-3, 0, 2, 10), sd = c(0.1, 1, 3)),
      delays = c(-40, -20, -10, -5, -2, -1, -0.5, 0, 0.5, 1, 2, 5, 10, 20, 40),
      pwindows = c(0.5, 1, 3, 10), kinks = numeric(0)
    ),
    list(
      pdist = pchisq, grid = expand.grid(df = c(0.5, 1, 3, 10, 40)),
      delays = c(0.001, 0.5, 2, 10, 40), pwindows = c(0.5, 3), kinks = 0
    ),
    list(
      pdist = pbeta,
      grid = expand.grid(shape1 = c(0.5, 1, 2, 5), shape2 = c(0.5, 1, 3, 10)),
      delays = c(0.001, 0.1, 0.3, 0.5, 0.9, 1, 1.5, 2, 3, 10),
      pwindows = c(0.3, 1, 3), kinks = c(0, 1)
    )
  )
  for (case in cases) {
    for (i in seq_len(nrow(case$grid))) {
      args <- as.list(case$grid[i, , drop = FALSE])
      obj <- do.call(new_pcens, c(list(case$pdist, dunif, list()), args))
      cdf <- function(t) do.call(case$pdist, c(list(t), args))
      for (pwindow in case$pwindows) {
        info <- paste(
          toString(paste(names(args), args, sep = " = ")), "pwindow", pwindow
        )
        expect_rel_equal(
          pcens_cdf(obj, case$delays, pwindow),
          reference_uniform_cdf(cdf, case$delays, pwindow, case$kinks),
          info = info
        )
      }
    }
  }
})

test_that("the new pcens_cdf methods match the numerical path, lie in [0, 1]
  and increase in q", {
  cases <- list(
    list(obj = new_pcens(pexp, dunif, list(), rate = 0.4), q = 0:30),
    list(
      obj = new_pcens(pnorm, dunif, list(), mean = 3, sd = 2),
      q = seq(-10, 20, by = 0.5)
    ),
    list(obj = new_pcens(pchisq, dunif, list(), df = 5), q = 0:40),
    list(
      obj = new_pcens(pbeta, dunif, list(), shape1 = 2, shape2 = 3),
      q = seq(0, 3, by = 0.1)
    )
  )
  for (case in cases) {
    for (pwindow in c(0.5, 2)) {
      info <- paste(class(case$obj)[1], "pwindow", pwindow)
      result <- pcens_cdf(case$obj, case$q, pwindow)
      expect_equal(
        result, pcens_cdf(case$obj, case$q, pwindow, use_numeric = TRUE),
        tolerance = 1e-6, info = info
      )
      expect_true(all(result >= 0), info = info)
      expect_true(all(result <= 1), info = info)
      expect_true(all(diff(result) >= -1e-12), info = info)
    }
  }
})

test_that("the exponential and chi-square methods agree with the gamma", {
  q <- c(0.001, 0.5, 2, 10, 40)
  for (pwindow in c(0.5, 2, 7)) {
    expect_equal(
      pcens_cdf(new_pcens(pexp, dunif, list(), rate = 0.3), q, pwindow),
      pcens_cdf(
        new_pcens(pgamma, dunif, list(), shape = 1, rate = 0.3), q, pwindow
      ),
      tolerance = 1e-10
    )
  }
  for (df in c(0.5, 3, 40)) {
    expect_identical(
      pcens_cdf(new_pcens(pchisq, dunif, list(), df = df), q, 3),
      pcens_cdf(
        new_pcens(pgamma, dunif, list(), shape = df / 2, scale = 2), q, 3
      ),
      info = paste("df", df)
    )
  }
})

test_that("the exponential and normal methods use the pdist defaults", {
  q <- c(-2, 0.5, 3, 10)
  expect_identical(
    pcens_cdf(new_pcens(pexp, dunif, list()), q, 2),
    pcens_cdf(new_pcens(pexp, dunif, list(), rate = 1), q, 2)
  )
  expect_identical(
    pcens_cdf(new_pcens(pnorm, dunif, list()), q, 2),
    pcens_cdf(new_pcens(pnorm, dunif, list(), mean = 0, sd = 1), q, 2)
  )
})

test_that("the normal method handles infinite q and is continuous across the
  switch to the asymptotic series", {
  obj <- new_pcens(pnorm, dunif, list(), mean = 1, sd = 2)
  expect_identical(pcens_cdf(obj, c(-Inf, Inf), 1), c(0, 1))
  obj <- new_pcens(pnorm, dunif, list(), mean = 0, sd = 1)
  expect_equal(
    pcens_cdf(obj, -10 + 1e-9 - 1e-7, 1), pcens_cdf(obj, -10 + 1e-9 + 1e-7, 1),
    tolerance = 1e-5
  )
  z <- c(-10.000001, -9.999999)
  expect_rel_equal(
    .norm_shortfall(z),
    vapply(z, function(zi) {
      stats::integrate(
        pnorm, -Inf, zi,
        rel.tol = 1.2e-14, abs.tol = 0, subdivisions = 500L
      )$value
    }, numeric(1)),
    tolerance = 1e-12
  )
})

test_that("the new pcens_cdf methods use the numerical path for arguments
  they do not handle, rather than ignoring them", {
  q <- c(0.5, 2, 6)
  # A non-zero ncp, lower.tail and log.p change what pdist returns
  objs <- list(
    new_pcens(pexp, dunif, list(), rate = 0.4, lower.tail = FALSE),
    new_pcens(pnorm, dunif, list(), mean = 2, log.p = TRUE),
    new_pcens(pchisq, dunif, list(), df = 3, lower.tail = FALSE),
    new_pcens(pchisq, dunif, list(), df = 3, ncp = 1),
    new_pcens(pbeta, dunif, list(), shape1 = 2, shape2 = 3, lower.tail = FALSE),
    new_pcens(pbeta, dunif, list(), shape1 = 2, shape2 = 3, ncp = 1)
  )
  for (obj in objs) {
    expect_identical(
      pcens_cdf(obj, q, 0.5), pcens_cdf(obj, q, 0.5, use_numeric = TRUE),
      info = toString(names(obj$args))
    )
  }
  expect_gt(
    max(abs(
      pcens_cdf(objs[[1]], q, 2) -
        pcens_cdf(new_pcens(pexp, dunif, list(), rate = 0.4), q, 2)
    )),
    0.1
  )
  # An unknown argument is an error in pdist
  exp_scale <- new_pcens(pexp, dunif, list(), scale = 2)
  expect_error(pcens_cdf(exp_scale, q, 2), "unused argument")
  # A zero ncp is the central distribution
  expect_identical(
    pcens_cdf(new_pcens(pchisq, dunif, list(), df = 3, ncp = 0), q, 1),
    pcens_cdf(new_pcens(pchisq, dunif, list(), df = 3), q, 1)
  )
})

test_that("the new pcens_cdf methods error for missing parameters", {
  expect_error(
    pcens_cdf(new_pcens(pchisq, dunif, list()), 1, 1),
    "df parameter is required for Chi-square distribution"
  )
  for (obj in list(
    new_pcens(pbeta, dunif, list(), shape2 = 2),
    new_pcens(pbeta, dunif, list(), shape1 = 2)
  )) {
    expect_error(
      pcens_cdf(obj, 1, 1),
      "shape1 and shape2 parameters are required for Beta distribution"
    )
  }
})

test_that("pprimarycensored and dprimarycensored work end to end for the new
  delays", {
  expect_identical(
    pprimarycensored(c(2, 5, 9), pexp, pwindow = 2, rate = 0.3),
    pcens_cdf(new_pcens(pexp, dunif, list(), rate = 0.3), c(2, 5, 9), 2)
  )
  pmf <- dprimarycensored(
    0:15, pnorm,
    pwindow = 2, swindow = 1, L = 0, D = 16, mean = 6, sd = 2
  )
  expect_equal(sum(pmf), 1, tolerance = 1e-10)
  # The window is not clipped at 0, so a normal delay has mass below 0
  lower <- pprimarycensored(-1, pnorm, pwindow = 2, mean = 2, sd = 2)
  expect_gt(lower, 0)
})

test_that("pprimarycensored matches the empirical CDF of rprimarycensored
  samples for the new delays", {
  set.seed(20260930)
  n <- 1e5
  probs <- c(0.05, 0.3, 0.6, 0.9)
  cases <- list(
    list(pdist = pexp, rdist = rexp, args = list(rate = 0.4)),
    list(pdist = pnorm, rdist = rnorm, args = list(mean = 2, sd = 1.5)),
    list(pdist = pchisq, rdist = rchisq, args = list(df = 4)),
    list(pdist = pbeta, rdist = rbeta, args = list(shape1 = 2, shape2 = 3))
  )
  for (case in cases) {
    for (pwindow in c(0.5, 2)) {
      samples <- do.call(
        rprimarycensored,
        c(list(n, case$rdist, pwindow = pwindow, swindow = 0), case$args)
      )
      delays <- unname(stats::quantile(samples, probs))
      actual <- do.call(
        pprimarycensored,
        c(list(delays, case$pdist, pwindow = pwindow), case$args)
      )
      # Within 4 standard errors of the probability at its own quantile
      expect_lt(
        max(abs(actual - probs) / sqrt(probs * (1 - probs) / n)), 4,
        label = sprintf("pwindow %g", pwindow)
      )
    }
  }
})
