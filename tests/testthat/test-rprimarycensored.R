test_that("rprimarycensored generates samples within the correct range", {
  n <- 1000
  pwindow <- 5
  D <- 10
  samples <- rpcens(
    n, rlnorm, pwindow,
    D = D, meanlog = 0, sdlog = 1
  )

  expect_true(all(samples >= 0 & samples < D))
})

test_that("rprimarycensored handles different primary distributions", {
  n <- 1000
  pwindow <- 5
  D <- 10
  r <- 0.5
  samples <- rpcens(
    n, rlnorm, pwindow,
    D = D, rprimary = rexpgrowth, rprimary_args = list(r = r),
    meanlog = 0, sdlog = 1
  )

  expect_true(all(samples >= 0 & samples < D))
})

test_that("rprimarycensored handles very truncated distributions", {
  n <- 1000
  pwindow <- 0.1
  D <- 1

  samples <- rpcens(
    n, rnorm, pwindow,
    D = D, mean = 0.5, sd = 1
  )

  expect_true(all(samples >= 0 & samples < D))
})

test_that(
  "rprimarycensored supports non-secondary event censored distributions",
  {
    n <- 1000
    pwindow <- 5
    D <- 10
    swindow <- 0

    samples <- rpcens(
      n, rlnorm, pwindow,
      swindow = swindow,
      D = D, meanlog = 0, sdlog = 1
    )

    expect_true(all(samples >= 0 & samples < D))

    # Check that samples are not rounded to discrete intervals
    expect_false(all(samples == floor(samples)))

    # Check that at least some samples have fractional parts
    expect_true(any(samples != floor(samples)))
  }
)

# Left truncation (L parameter) tests

test_that("rprimarycensored generates samples within [L, D) with L > 0", {
  n <- 1000
  pwindow <- 5
  D <- 10
  L <- 2
  samples <- rpcens(
    n, rlnorm, pwindow,
    D = D, L = L, meanlog = 1, sdlog = 0.5
  )

  expect_true(all(samples >= L & samples < D))
})

test_that("rprimarycensored errors when L >= D", {
  expect_error(
    rpcens(10, rlnorm, pwindow = 1, D = 10, L = 10, meanlog = 1, sdlog = 1),
    "L must be less than D"
  )
  expect_error(
    rpcens(10, rlnorm, pwindow = 1, D = 10, L = 15, meanlog = 1, sdlog = 1),
    "L must be less than D"
  )
})

test_that("rprimarycensored errors when L < 0", {
  expect_error(
    rpcens(10, rlnorm, pwindow = 1, D = 10, L = -1, meanlog = 1, sdlog = 1),
    "L must be non-negative"
  )
})

test_that("rprimarycensored with L = 0 matches default behaviour", {
  set.seed(123)
  samples_default <- rpcens(
    100, rlnorm,
    pwindow = 1, D = 10,
    meanlog = 1, sdlog = 1
  )
  set.seed(123)
  samples_explicit <- rpcens(
    100, rlnorm,
    pwindow = 1, D = 10, L = 0,
    meanlog = 1, sdlog = 1
  )
  expect_identical(samples_default, samples_explicit)
})

test_that("rprimarycensored with L > 0 rounds correctly to swindow", {
  n <- 100
  pwindow <- 1
  D <- 10
  L <- 2
  swindow <- 1

  samples <- rpcens(
    n, rlnorm, pwindow,
    swindow = swindow, D = D, L = L, meanlog = 1, sdlog = 0.5
  )

  # All samples should be integers (multiples of swindow)
  expect_true(all(samples == floor(samples)))
  # All samples should be >= L
  expect_true(all(samples >= L))
})
