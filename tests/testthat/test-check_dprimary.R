test_that("check_dprimary works correctly for valid PDF", {
  expect_silent(check_dprimary(dunif, pwindow = 1))
})

test_that("check_dprimary throws error if PDF doesn't have min and max args", {
  expect_error(
    check_dprimary(dnorm, pwindow = 1),
    "dprimary must take min and max as arguments"
  )
})

test_that("check_dprimary throws error for invalid PDF", {
  invalid_pdf <- function(x, min, max, ...) rep(0.5, length(x))
  expect_error(
    check_dprimary(invalid_pdf, pwindow = 1),
    "dprimary is not a valid probability density function"
  )
})

test_that("check_dprimary respects tolerance", {
  almost_valid_pdf <- function(x, min, max, ...) dunif(x, min, max) * 1.01
  expect_silent(check_dprimary(almost_valid_pdf, pwindow = 1, tolerance = 0.02))
  expect_error(
    check_dprimary(almost_valid_pdf, pwindow = 1, tolerance = 0.005),
    "dprimary is not a valid probability density function"
  )
})
