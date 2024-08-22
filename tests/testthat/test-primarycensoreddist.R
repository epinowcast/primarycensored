library(testthat)
library(primarycensoreddist)

test_that("pprimarycensoreddist produces valid CDF", {
  # Use exponential distribution with rate 1
  dist_func <- pexp
  pwindow <- 1
  D <- Inf
  q <- seq(0, 5, by = 0.5)

  cdf <- pprimarycensoreddist(q, dist_func, pwindow = pwindow, D = D, rate = 1)

  # Check that CDF is monotonically increasing
  expect_true(all(diff(cdf) >= 0))

  # Check that CDF starts at 0 and approaches 1
  expect_equal(cdf[1], 0)
  expect_true(cdf[length(cdf)] > 0.99)

  # Check that CDF is between 0 and 1
  expect_true(all(cdf >= 0 & cdf <= 1))
})

test_that("dprimarycensoreddist produces valid PMF", {
  # Use exponential distribution with rate 1
  dist_func <- pexp
  pwindow <- 1
  D <- Inf
  swindow <- 0.5
  x <- seq(0, 5, by = swindow)

  pmf <- dprimarycensoreddist(x, dist_func, pwindow = pwindow, swindow = swindow, D = D, rate = 1)

  # Check that PMF is non-negative
  expect_true(all(pmf >= 0))

  # Check that PMF sums to approximately 1
  expect_equal(sum(pmf), 1, tolerance = 1e-6)
})

test_that("pprimarycensoreddist and dprimarycensoreddist are consistent", {
  # Use exponential distribution with rate 1
  dist_func <- pexp
  pwindow <- 1
  D <- Inf
  swindow <- 0.5
  x <- seq(0, 5, by = swindow)

  cdf <- pprimarycensoreddist(x, dist_func, pwindow = pwindow, D = D, rate = 1)
  pmf <- dprimarycensoreddist(x, dist_func, pwindow = pwindow, swindow = swindow, D = D, rate = 1)

  # Check that the difference of CDF equals PMF
  cdf_diff <- c(cdf[1], diff(cdf))
  expect_equal(cdf_diff, pmf, tolerance = 1e-6)
})

test_that("pprimarycensoreddist handles truncation correctly", {
  dist_func <- pexp
  pwindow <- 1
  D <- 3
  q <- seq(0, 5, by = 0.5)

  cdf_truncated <- pprimarycensoreddist(q, dist_func, pwindow = pwindow, D = D, rate = 1)
  cdf_untruncated <- pprimarycensoreddist(q, dist_func, pwindow = pwindow, D = Inf, rate = 1)

  # Check that truncated CDF reaches 1 at D
  expect_equal(cdf_truncated[q >= D], 1)

  # Check that truncated CDF is greater than or equal to untruncated CDF
  expect_true(all(cdf_truncated >= cdf_untruncated))
})

test_that("primary distribution affects the result", {
  dist_func <- pexp
  pwindow <- 1
  D <- Inf
  q <- seq(0, 5, by = 0.5)

  cdf_unif <- pprimarycensoreddist(q, dist_func, pwindow = pwindow, D = D, rate = 1)
  cdf_exp <- pprimarycensoreddist(q, dist_func, pwindow = pwindow, D = D, rate = 1,
                                  primary_dist = exp_primary_dist,
                                  primary_args = list(rate = 2))

  # Check that the CDFs are different
  expect_false(all(cdf_unif == cdf_exp))
})