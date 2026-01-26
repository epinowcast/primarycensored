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
    c(1, 5, 10),
    plnorm,
    pwindow,
    D = D,
    meanlog = 1,
    sdlog = 1
  )
  cdf_expgrowth <- ppcens(
    c(1, 5, 10),
    plnorm,
    pwindow,
    D = D,
    dprimary = dexpgrowth,
    dprimary_args = list(r = 0.2),
    meanlog = 1,
    sdlog = 1
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

# Left truncation (L parameter) tests

test_that("pprimarycensored returns 0 for quantiles <= L", {
  pwindow <- 1
  D <- 10
  L <- 2
  cdf <- ppcens(
    c(0, 1, 2), plnorm, pwindow,
    D = D, L = L, meanlog = 1, sdlog = 1
  )
  expect_identical(cdf, c(0, 0, 0))
})
test_that("pprimarycensored returns 1 for quantiles >= D with L > 0", {
  pwindow <- 1
  D <- 10
  L <- 2
  cdf <- ppcens(
    c(10, 15), plnorm, pwindow,
    D = D, L = L, meanlog = 1, sdlog = 1
  )
  expect_identical(cdf, c(1, 1))
})

test_that("pprimarycensored is monotonically increasing with L > 0", {
  pwindow <- 1
  D <- 10
  L <- 2
  q <- seq(L, D, by = 0.5)
  cdf <- ppcens(q, plnorm, pwindow, D = D, L = L, meanlog = 1, sdlog = 1)
  expect_true(all(diff(cdf) >= 0))
})

test_that("pprimarycensored errors when L >= D", {
  expect_error(
    ppcens(5, plnorm, pwindow = 1, D = 10, L = 10, meanlog = 1, sdlog = 1),
    "L must be less than D"
  )
  expect_error(
    ppcens(5, plnorm, pwindow = 1, D = 10, L = 15, meanlog = 1, sdlog = 1),
    "L must be less than D"
  )
})

test_that("pprimarycensored errors when L < 0", {
  expect_error(
    ppcens(5, plnorm, pwindow = 1, D = 10, L = -1, meanlog = 1, sdlog = 1),
    "L must be non-negative"
  )
})

test_that("pprimarycensored with L = 0 matches default behaviour", {
  pwindow <- 1
  D <- 10
  q <- seq(0, D, by = 0.5)
  cdf_default <- ppcens(q, plnorm, pwindow, D = D, meanlog = 1, sdlog = 1)
  cdf_explicit <- ppcens(
    q, plnorm, pwindow,
    D = D, L = 0, meanlog = 1, sdlog = 1
  )
  expect_identical(cdf_default, cdf_explicit)
})

test_that("pprimarycensored with L > 0 shifts distribution correctly", {
  pwindow <- 1
  D <- 10
  L <- 2

  # CDF at L should be 0
  cdf_at_L <- ppcens(
    L, plnorm, pwindow,
    D = D, L = L, meanlog = 1, sdlog = 1
  )
  expect_identical(cdf_at_L, 0)

  # CDF at D should be 1
  cdf_at_D <- ppcens(
    D, plnorm, pwindow,
    D = D, L = L, meanlog = 1, sdlog = 1
  )
  expect_identical(cdf_at_D, 1)

  # CDF between L and D should be between 0 and 1
  mid_cdf <- ppcens(6, plnorm, pwindow, D = D, L = L, meanlog = 1, sdlog = 1)
  expect_gt(mid_cdf, 0)
  expect_lt(mid_cdf, 1)
})

test_that("pprimarycensored works with L > 0 and D = Inf", {
  pwindow <- 1
  L <- 2

  # CDF at L should be 0
  cdf_at_L <- ppcens(L, plnorm, pwindow, D = Inf, L = L, meanlog = 1, sdlog = 1)
  expect_identical(cdf_at_L, 0)

  # CDF above L should be positive
  cdf_above_L <- ppcens(5, plnorm, pwindow, D = Inf, L = L, meanlog = 1, sdlog = 1)
  expect_gt(cdf_above_L, 0)
  expect_lt(cdf_above_L, 1)

  # CDF should approach 1 for large values
  cdf_large <- ppcens(50, plnorm, pwindow, D = Inf, L = L, meanlog = 1, sdlog = 1)
  expect_gt(cdf_large, 0.99)
})
