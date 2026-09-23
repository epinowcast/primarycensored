test_that("pcens_pmf matches dprimarycensored", {
  obj <- new_pcens(
    pweibull, dexpgrowth, list(r = 0.3), shape = 1.5, scale = 2
  )
  x <- 0:9
  expected <- dprimarycensored(
    x, pweibull,
    pwindow = 2, D = 10, dprimary = dexpgrowth,
    primary_args = list(r = 0.3), shape = 1.5, scale = 2
  )
  expect_identical(pcens_pmf(obj, x, pwindow = 2, D = 10), expected)
})

test_that("pcens_pmf sums to one over [L, D)", {
  obj <- new_pcens(pgamma, dunif, list(), shape = 2, scale = 1.5)
  expect_equal(
    sum(pcens_pmf(obj, 0:14, pwindow = 1, D = 15)), 1,
    tolerance = 1e-10
  )
  expect_equal(
    sum(pcens_pmf(obj, 2:14, pwindow = 1, L = 2, D = 15)), 1,
    tolerance = 1e-10
  )
})

test_that("pcens_pmf clips the secondary window at D", {
  obj <- new_pcens(pgamma, dunif, list(), shape = 2, scale = 1.5)
  expect_message(
    pcens_pmf(obj, 9, pwindow = 1, swindow = 2, D = 10),
    "clipping the upper end"
  )
  pmf <- suppressMessages(
    pcens_pmf(obj, 9, pwindow = 1, swindow = 2, D = 10)
  )
  cdf <- pcens_cdf(obj, c(9, 10), pwindow = 1)
  expect_equal(pmf, diff(cdf) / cdf[2], tolerance = 1e-10)
})

test_that("pcens_pmf errors for x outside [L, D)", {
  obj <- new_pcens(pgamma, dunif, list(), shape = 2, scale = 1.5)
  expect_error(pcens_pmf(obj, 0:3, pwindow = 1, L = 1), "below L")
  expect_error(pcens_pmf(obj, 0:10, pwindow = 1, D = 10), "strictly less")
  expect_error(pcens_pmf(obj, 1, pwindow = 1, L = 5, D = 2), "less than D")
  expect_error(pcens_pmf(list(), 1, pwindow = 1), "pcens")
})

test_that("pcens_pmf gives zero mass at infinite x", {
  # Analytical CDF methods return NaN at Inf, so pcens_pmf() sets the CDF to
  # 0 and 1 at -Inf and Inf as pprimarycensored() does.
  obj <- new_pcens(pgamma, dunif, list(), shape = 2, scale = 1.5)
  pmf <- pcens_pmf(obj, c(-Inf, 0, Inf), pwindow = 1)
  expect_identical(pmf[c(1, 3)], c(0, 0))
  expect_identical(pmf[2], pcens_cdf(obj, 1, pwindow = 1))
})

test_that("pcens_pmf dispatches to the default method", {
  obj <- new_pcens(pgamma, dunif, list(), shape = 2, scale = 1.5)
  expect_true(is.function(pcens_pmf.default))
  expect_identical(
    pcens_pmf(obj, 0:5, pwindow = 1, D = 10),
    pcens_pmf.default(obj, 0:5, pwindow = 1, D = 10)
  )
  # A method for a more specific class takes precedence over the default
  pcens_pmf.pcens_test <- function(object, x, pwindow, ...) {
    rep(-1, length(x))
  }
  class(obj) <- c("pcens_test", class(obj))
  expect_identical(pcens_pmf(obj, 0:2, pwindow = 1), c(-1, -1, -1))
})
