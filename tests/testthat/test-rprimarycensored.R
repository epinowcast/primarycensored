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
    L = 0, D = D, mean = 0.5, sd = 1
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

test_that("rprimarycensored default matches L = 0 for positive support", {
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

test_that("rprimarycensored signed-support samples match analytical CDF", {
  set.seed(42)
  L <- -0.5
  D <- 1.5
  n <- 20000
  samples <- rpcens(
    n = n, rnorm,
    pwindow = 1, swindow = 0, L = L, D = D, mean = 0, sd = 1
  )
  expect_length(samples, n)
  expect_true(all(samples >= L & samples < D))

  qs <- seq(L + 0.1, D - 0.1, length.out = 5)
  empirical <- vapply(qs, function(q) mean(samples <= q), numeric(1))
  analytical <- ppcens(
    qs, pnorm,
    pwindow = 1, L = L, D = D, mean = 0, sd = 1
  )
  expect_equal(empirical, analytical, tolerance = 0.02)
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
