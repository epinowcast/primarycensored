test_that("rprimarycensoreddist generates samples within the correct range", {
  n <- 1000
  pwindow <- 5
  D <- 10
  samples <- rpcens(
    n, rlnorm, pwindow,
    D = D, meanlog = 0, sdlog = 1
  )

  expect_true(all(samples > 0 & samples <= D))
})
