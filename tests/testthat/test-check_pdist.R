test_that("check_pdist works correctly for valid CDF", {
  expect_silent(check_pdist(pnorm, D = 10))
  expect_silent(check_pdist(punif, D = 1))
})

test_that("check_pdist throws error for invalid CDF", {
  invalid_cdf <- function(x, ...) x^2 - 1 # Not bounded between 0 and 1
  expect_error(
    check_pdist(invalid_cdf, D = 10),
    "pdist is not a valid cumulative distribution function"
  )
})

test_that("check_pdist handles infinite D correctly", {
  expect_silent(check_pdist(pnorm, D = Inf))
})
