test_that("dprimarycensoreddist sums to 1 for discrete values", {
  pwindow <- 1
  D <- 10
  pmf <- dpcens(
    0:(D - 1), plnorm, pwindow,
    D = D, meanlog = 1, sdlog = 1
  )
  expect_equal(sum(pmf), 1, tolerance = 1e-6)
})

test_that("dprimarycensoreddist handles log probabilities", {
  pwindow <- 1
  D <- 10
  pmf <- dpcens(
    0:(D - 1), plnorm, pwindow,
    D = D, meanlog = 1, sdlog = 1
  )
  log_pmf <- dpcens(
    0:(D - 1), plnorm, pwindow,
    D = D, meanlog = 1, sdlog = 1, log = TRUE
  )
  expect_equal(exp(log_pmf), pmf, tolerance = 1e-6)
})

test_that("dprimarycensoreddist handles non-finite D", {
  pwindow <- 1
  D <- Inf
  pmf <- dpcens(
    0:100, plnorm, pwindow,
    D = D, meanlog = 1, sdlog = 1
  )
  expect_lt(sum(pmf), 1)
  expect_equal(sum(pmf), 1, tolerance = 0.01)
})

test_that("dprimarycensoreddist matches difference of pprimarycensoreddist", {
  x <- c(1, 2, 3)
  pwindow <- 5
  swindow <- 0.5
  D <- 10

  pmf <- dpcens(
    x, plnorm, pwindow, swindow, D,
    meanlog = 0, sdlog = 1
  )
  cdf_diff <- sapply(x, function(xi) {
    ppcens(
      xi + swindow, plnorm, pwindow, D,
      meanlog = 0, sdlog = 1
    ) -
      ppcens(
        xi, plnorm, pwindow, D,
        meanlog = 0, sdlog = 1
      )
  })

  expect_equal(pmf, cdf_diff, tolerance = 1e-6)
})

test_that("dprimarycensoreddist throws an error for invalid upper truncation point", {
  d <- 10
  pwindow <- 1
  swindow <- 1
  D <- 10

  expect_error(
    dpcens(
      d, plnorm, pwindow, swindow, D,
      meanlog = 0, sdlog = 1
    ),
    "Upper truncation point is greater than D"
  )
})
