test_that("rprimarycensoreddist generates samples within the correct range", {
  n <- 1000
  pwindow <- 5
  D <- 10
  samples <- rpcens(
    n, rlnorm, pwindow,
    D = D, meanlog = 0, sdlog = 1
  )

  expect_true(all(samples >= 0 & samples < D))
})

test_that("rprimarycensoreddist handles different primary distributions", {
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
