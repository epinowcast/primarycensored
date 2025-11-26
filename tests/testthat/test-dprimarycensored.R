test_that("dprimarycensored sums to 1 for discrete values", {
  pwindow <- 1
  D <- 10
  pmf <- dpcens(
    0:(D - 1), plnorm, pwindow,
    D = D, meanlog = 1, sdlog = 1
  )
  expect_equal(sum(pmf), 1, tolerance = 1e-6)
})

test_that("dprimarycensored handles log probabilities", {
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

test_that("dprimarycensored handles non-finite D", {
  pwindow <- 1
  D <- Inf
  pmf <- dpcens(
    0:100, plnorm, pwindow,
    D = D, meanlog = 1, sdlog = 1
  )
  expect_lt(sum(pmf), 1)
  expect_equal(sum(pmf), 1, tolerance = 0.01)
})

test_that("dprimarycensored matches difference of pprimarycensored", {
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

test_that(
  "dprimarycensored throws an error for invalid upper truncation point",
  {
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
  }
)

test_that("dprimarycensored returns 0 for negative d", {
  d <- -1
  pwindow <- 1
  swindow <- 0.5
  D <- 10

  expect_identical(
    dpcens(
      d, plnorm, pwindow, swindow, D,
      meanlog = 0, sdlog = 1
    ), 0
  )
  expect_identical(
    dpcens(
      c(8, d), plnorm, pwindow, swindow, D,
      meanlog = 0, sdlog = 1
    )[2], 0
  )
})

test_that("dprimarycensored returns non-negative values", {
  # Test case from issue #238: lognormal with specific parameters

  # that previously produced negative values due to floating-point precision
  pmf <- dpcens(
    x = seq(0, 29), plnorm, pwindow = 1, swindow = 1, D = 30,
    meanlog = 0.55, sdlog = 0.27
  )
  expect_true(all(pmf >= 0), info = "PMF should never be negative")

  # Also test with infinite D

  pmf_inf <- dpcens(
    x = seq(0, 29), plnorm, pwindow = 1, swindow = 1, D = Inf,
    meanlog = 0.55, sdlog = 0.27
  )
  expect_true(
    all(pmf_inf >= 0),
    info = "PMF with D=Inf should never be negative"
  )

  # Test with other distributions
  pmf_gamma <- dpcens(
    x = seq(0, 29), pgamma, pwindow = 1, swindow = 1, D = 30,
    shape = 2, rate = 0.5
  )
  expect_true(all(pmf_gamma >= 0), info = "Gamma PMF should never be negative")
})
