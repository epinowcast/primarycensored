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
    x, plnorm,
    pwindow = pwindow, swindow = swindow, D = D,
    meanlog = 0, sdlog = 1
  )
  cdf_diff <- sapply(x, function(xi) {
    ppcens(
      xi + swindow, plnorm,
      pwindow = pwindow, D = D,
      meanlog = 0, sdlog = 1
    ) -
      ppcens(
        xi, plnorm,
        pwindow = pwindow, D = D,
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
        d, plnorm,
        pwindow = pwindow, swindow = swindow, D = D,
        meanlog = 0, sdlog = 1
      ),
      "Upper truncation point is greater than D"
    )
  }
)

test_that("dprimarycensored errors for negative d", {
  d <- -1
  pwindow <- 1
  swindow <- 0.5
  D <- 10

  expect_error(
    dpcens(
      d, plnorm,
      pwindow = pwindow, swindow = swindow, D = D,
      meanlog = 0, sdlog = 1
    ),
    "values of x are below L"
  )
  expect_error(
    dpcens(
      c(8, d), plnorm,
      pwindow = pwindow, swindow = swindow, D = D,
      meanlog = 0, sdlog = 1
    ),
    "values of x are below L"
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

  # Test with exponential growth primary distribution
  pmf_expgrowth <- dpcens(
    x = seq(0, 29), plnorm, pwindow = 1, swindow = 1, D = 30,
    dprimary = dexpgrowth,
    dprimary_args = list(r = 0.2),
    meanlog = 0.55, sdlog = 0.27
  )
  expect_true(
    all(pmf_expgrowth >= 0),
    info = "PMF with expgrowth primary should never be negative"
  )
})

# Left truncation (L parameter) tests

test_that("dprimarycensored sums to 1 over [L, D) with L > 0", {
  pwindow <- 1
  D <- 10
  L <- 2
  pmf <- dpcens(
    L:(D - 1), plnorm, pwindow,
    D = D, L = L, meanlog = 1, sdlog = 1
  )
  expect_equal(sum(pmf), 1, tolerance = 1e-6)
})

test_that("dprimarycensored errors for x < L", {
  pwindow <- 1
  D <- 10
  L <- 2

  expect_error(
    dpcens(
      0:(L - 1), plnorm, pwindow,
      D = D, L = L, meanlog = 1, sdlog = 1
    ),
    "values of x are below L"
  )
})

test_that("dprimarycensored errors when L >= D", {
  expect_error(
    dpcens(5, plnorm, pwindow = 1, D = 10, L = 10, meanlog = 1, sdlog = 1),
    "L must be less than D"
  )
  expect_error(
    dpcens(5, plnorm, pwindow = 1, D = 10, L = 15, meanlog = 1, sdlog = 1),
    "L must be less than D"
  )
})

test_that("dprimarycensored errors when L < 0", {
  expect_error(
    dpcens(5, plnorm, pwindow = 1, D = 10, L = -1, meanlog = 1, sdlog = 1),
    "L must be non-negative"
  )
})

test_that("dprimarycensored with L = 0 matches default behaviour", {
  pwindow <- 1
  D <- 10
  pmf_default <- dpcens(
    0:(D - 1), plnorm, pwindow,
    D = D, meanlog = 1, sdlog = 1
  )
  pmf_explicit <- dpcens(
    0:(D - 1), plnorm, pwindow,
    D = D, L = 0, meanlog = 1, sdlog = 1
  )
  expect_identical(pmf_default, pmf_explicit)
})

test_that("dprimarycensored is consistent with pprimarycensored for L > 0", {
  pwindow <- 1
  D <- 10
  L <- 2
  x <- L:(D - 2)
  pmf <- dpcens(x, plnorm, pwindow, D = D, L = L, meanlog = 1, sdlog = 1)
  cdf_diff <- sapply(x, function(xi) {
    ppcens(
      xi + 1, plnorm, pwindow,
      D = D, L = L,
      meanlog = 1, sdlog = 1
    ) -
      ppcens(
        xi, plnorm, pwindow,
        D = D, L = L,
        meanlog = 1, sdlog = 1
      )
  })
  expect_equal(pmf, cdf_diff, tolerance = 1e-6)
})

test_that("dprimarycensored handles log probabilities with L > 0", {
  pwindow <- 1
  D <- 10
  L <- 2
  pmf <- dpcens(
    L:(D - 1), plnorm, pwindow,
    D = D, L = L, meanlog = 1, sdlog = 1
  )
  log_pmf <- dpcens(
    L:(D - 1), plnorm, pwindow,
    D = D, L = L, meanlog = 1, sdlog = 1, log = TRUE
  )
  expect_equal(exp(log_pmf), pmf, tolerance = 1e-6)
})
