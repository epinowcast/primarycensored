skip_on_cran()
if (on_ci()) {
  skip_on_os("windows")
  skip_on_os("mac")
}


test_that("Stan expgrowth_pdf matches R dexpgrowth", {
  x <- seq(0, 1, by = 0.1)
  min <- 0
  max <- 1
  r <- 0.5

  stan_pdf <- sapply(x, expgrowth_pdf, min, max, r)
  r_pdf <- dexpgrowth(x, min, max, r)

  expect_equal(stan_pdf, r_pdf, tolerance = 1e-4)
})

test_that("Stan expgrowth_lpdf matches R dexpgrowth with log = TRUE", {
  x <- seq(0, 1, by = 0.1)
  min <- 0
  max <- 1
  r <- 0.5

  stan_lpdf <- sapply(x, expgrowth_lpdf, min, max, r)
  r_lpdf <- dexpgrowth(x, min, max, r, log = TRUE)

  expect_equal(stan_lpdf, r_lpdf, tolerance = 1e-4)
})

test_that("Stan expgrowth_cdf matches R pexpgrowth", {
  x <- seq(0, 1, by = 0.1)
  min <- 0
  max <- 1
  r <- 0.5

  stan_cdf <- sapply(x, expgrowth_cdf, min, max, r)
  r_cdf <- pexpgrowth(x, min, max, r)

  expect_equal(stan_cdf, r_cdf, tolerance = 1e-4)
})

test_that("Stan expgrowth_lcdf matches R pexpgrowth with log.p = TRUE", {
  x <- seq(0, 1, by = 0.1)
  min <- 0
  max <- 1
  r <- 0.5

  stan_lcdf <- sapply(x, expgrowth_lcdf, min, max, r)
  r_lcdf <- pexpgrowth(x, min, max, r, log.p = TRUE)

  expect_equal(stan_lcdf, r_lcdf, tolerance = 1e-6)
})

test_that("Stan expgrowth_rng matches distribution of R rexpgrowth", {
  n <- 10000
  min <- 0
  max <- 1
  r <- 0.5

  set.seed(123)
  stan_samples <- replicate(n, expgrowth_rng(min, max, r))
  set.seed(123)
  r_samples <- rexpgrowth(n, min, max, r)

  expect_equal(mean(stan_samples), mean(r_samples), tolerance = 1e-2)
  expect_equal(sd(stan_samples), sd(r_samples), tolerance = 1e-2)
  expect_equal(quantile(stan_samples), quantile(r_samples), tolerance = 1e-2)
})
