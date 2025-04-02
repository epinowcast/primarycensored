test_that("add_name_attribute adds name attribute correctly", {
  # Test with base R distribution function
  dist <- add_name_attribute(pnorm, "pnorm")
  expect_identical(attr(dist, "name"), "pnorm")

  # Test with custom function
  custom_func <- function(x) x^2
  named_func <- add_name_attribute(custom_func, "square")
  expect_identical(attr(named_func, "name"), "square")
  expect_identical(named_func(2), custom_func(2))
})

test_that("add_name_attribute preserves function behaviour", {
  # Test that the function still works as expected after adding name
  dist <- add_name_attribute(pnorm, "test")
  x <- 1.96
  expect_identical(dist(x), pnorm(x))

  # Test with additional arguments
  expect_identical(
    dist(x, mean = 0, sd = 2),
    pnorm(x, mean = 0, sd = 2)
  )
})
