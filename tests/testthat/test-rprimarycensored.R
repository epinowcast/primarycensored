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
