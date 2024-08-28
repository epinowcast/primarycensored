skip_on_cran()
skip_on_os("windows")
skip_on_os("mac")
skip_on_local()

test_that("Stan expgrowth functions match R implementations", {
  min <- 0
  max <- 10
  r <- 0.5
  x <- seq(min, max, length.out = 100)

  r_pdf <- dexpgrowth(x, min, max, r)
  r_cdf <- pexpgrowth(x, min, max, r)

  stan_pdf <- sapply(x, expgrowth_pdf, min = min, max = max, r = r)
  stan_cdf <- sapply(x, expgrowth_cdf, min = min, max = max, r = r)

  expect_equal(r_pdf, stan_pdf, tolerance = 1e-6)
  expect_equal(r_cdf, stan_cdf, tolerance = 1e-6)
})

test_that("Stan expgrowth_rng matches R rexpgrowth", {
  min <- 0
  max <- 10
  r <- 0.5
  n <- 10000

  set.seed(123)
  r_samples <- rexpgrowth(n, min, max, r)

  set.seed(123)
  stan_samples <- replicate(n, expgrowth_rng(min, max, r))

  ks_test <- ks.test(r_samples, stan_samples)
  expect_gt(ks_test$p.value, 0.05)
})

test_that("Stan expgrowth functions handle edge cases", {
  x <- 5
  min <- 0
  max <- 10
  r <- 1e-11

  r_pdf <- dexpgrowth(x, min, max, r)
  r_cdf <- pexpgrowth(x, min, max, r)

  stan_pdf <- expgrowth_pdf(x, min, max, r)
  stan_cdf <- expgrowth_cdf(x, min, max, r)

  expect_equal(r_pdf, stan_pdf, tolerance = 1e-6)
  expect_equal(r_cdf, stan_cdf, tolerance = 1e-6)

  x_below <- -1
  x_above <- 11

  expect_identical(dexpgrowth(x_below, min, max, r), 0)
  expect_identical(dexpgrowth(x_above, min, max, r), 0)
  expect_identical(pexpgrowth(x_below, min, max, r), 0)
  expect_identical(pexpgrowth(x_above, min, max, r), 1)

  expect_identical(expgrowth_pdf(x_below, min, max, r), 0)
  expect_identical(expgrowth_pdf(x_above, min, max, r), 0)
  expect_identical(expgrowth_cdf(x_below, min, max, r), 0)
  expect_identical(expgrowth_cdf(x_above, min, max, r), 1)
})

test_that("Stan expgrowth log functions match R implementations", {
  min <- 0
  max <- 10
  r <- 0.5
  x <- seq(min, max, length.out = 100)

  r_log_pdf <- dexpgrowth(x, min, max, r, log = TRUE)
  r_log_cdf <- pexpgrowth(x, min, max, r, log.p = TRUE)

  stan_log_pdf <- sapply(x, expgrowth_lpdf, min = min, max = max, r = r)
  stan_log_cdf <- sapply(x, expgrowth_lcdf, min = min, max = max, r = r)

  expect_equal(r_log_pdf, stan_log_pdf, tolerance = 1e-6)
  expect_equal(r_log_cdf, stan_log_cdf, tolerance = 1e-6)
})
