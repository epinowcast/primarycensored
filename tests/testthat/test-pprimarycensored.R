test_that("pprimarycensored returns 0 for non-positive quantiles", {
  pwindow <- 1
  D <- 10
  cdf <- ppcens(c(-1, 0), plnorm, pwindow, D = D, meanlog = 1, sdlog = 1)
  expect_identical(cdf, c(0, 0))
})

test_that("pprimarycensored approaches 1 for large quantiles", {
  pwindow <- 1
  D <- Inf
  cdf <- ppcens(1000, plnorm, pwindow, D = D, meanlog = 1, sdlog = 1)
  expect_equal(cdf, 1, tolerance = 1e-6)
})

test_that("pprimarycensored is monotonically increasing", {
  pwindow <- 1
  D <- 10
  q <- seq(0, D, by = 0.5)
  cdf <- ppcens(q, plnorm, pwindow, D = D, meanlog = 1, sdlog = 1)
  expect_true(all(diff(cdf) >= 0))
})

test_that("pprimarycensored handles finite D correctly", {
  pwindow <- 1
  D <- 10
  cdf <- ppcens(D, plnorm, pwindow, D = D, meanlog = 1, sdlog = 1)
  expect_equal(cdf, 1, tolerance = 1e-6)
})

test_that("pprimarycensored handles custom primary distributions", {
  pwindow <- 5
  D <- 20
  cdf_uniform <- ppcens(
    c(1, 5, 10), plnorm, pwindow,
    D = D,
    meanlog = 1, sdlog = 1
  )
  cdf_expgrowth <- ppcens(
    c(1, 5, 10), plnorm, pwindow,
    D = D,
    dprimary = dexpgrowth, dprimary_args = list(r = 0.2),
    meanlog = 1, sdlog = 1
  )
  expect_false(all(cdf_uniform == cdf_expgrowth))
})

test_that("pprimarycensored is consistent with dprimarycensored", {
  pwindow <- 1
  D <- 10
  q <- 0:9
  cdf <- ppcens(q + 1, plnorm, pwindow, D = D, meanlog = 1, sdlog = 1)
  pmf <- dpcens(q, plnorm, pwindow, D = D, meanlog = 1, sdlog = 1)
  cdf_from_pmf <- cumsum(pmf)
  expect_equal(cdf, cdf_from_pmf, tolerance = 1e-6)
})
