test_that("check_truncation is silent when D is appropriate", {
  expect_silent(
    check_truncation(delays = c(1, 2, 3, 4), D = 5, multiplier = 2)
  )
})

test_that("check_truncation gives a message when D is too large", {
  expect_message(
    check_truncation(delays = c(1, 2, 3, 4), D = 20, multiplier = 2),
    paste0(
      "The truncation time D \\(20\\) is larger than 2 times the maximum ", # nolint
      "observed delay \\(4\\). Consider setting D to Inf for better ", # nolint
      "efficiency with minimal accuracy cost for this case."
    )
  )
})

test_that("check_truncation warns for empty delays vector", {
  expect_warning(
    check_truncation(delays = numeric(0), D = 10, multiplier = 2),
    "No finite observed delays to check."
  )
})

test_that("check_truncation handles delays with NA and Inf values", {
  expect_silent(
    check_truncation(
      delays = c(1, 2, NA, 4, Inf), D = 10, multiplier = 2
    )
  )
})

test_that("check_truncation errors on non-numeric delays", {
  expect_error(
    check_truncation(delays = "not numeric", D = 10, multiplier = 2),
    "All arguments must be numeric."
  )
})

test_that("check_truncation errors on negative D", {
  expect_error(
    check_truncation(delays = c(1, 2, 3, 4), D = -1, multiplier = 2),
    paste0(
      "Invalid argument values. D must be positive and multiplier ",
      "must be greater than 1."
    )
  )
})

test_that("check_truncation errors on multiplier <= 1", {
  expect_error(
    check_truncation(delays = c(1, 2, 3, 4), D = 10, multiplier = 0.5),
    paste0(
      "Invalid argument values. D must be positive and multiplier ",
      "must be greater than 1."
    )
  )
})

test_that("check_truncation is silent when D is infinite", {
  expect_silent(
    check_truncation(delays = c(1, 2, 3, 4), D = Inf, multiplier = 2)
  )
})
