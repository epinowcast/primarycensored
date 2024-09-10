# Skip tests if fitdistrplus is not installed
if (!requireNamespace("fitdistrplus", quietly = TRUE)) {
  skip("Package 'fitdistrplus' is required for these tests")
}

test_that("fitdistdoublecens works correctly", {
  # Set seed for reproducibility
  set.seed(123)

  # Define true distribution parameters
  n <- 1000
  shape <- 1.77
  rate <- 0.44

  # Generate samples
  samples <- rprimarycensoreddist(
    n, rgamma,
    shape = shape, rate = rate,
    pwindow = 1, swindow = 1, D = 8
  )

  # Create data frame
  delay_data <- data.frame(
    left = samples,
    right = samples + swindow
  )

  # Fit the model using fitdistdoublecens
  fit <- fitdistdoublecens(
    delay_data,
    distr = "gamma",
    start = list(shape = 1, rate = 1),
    D = 8, pwindow = 1
  )

  # Check that the function returns a fitdist object
  expect_s3_class(fit, "fitdist")

  # Check that the estimated parameters are close to the true values
  expect_equal(unname(fit$estimate["shape"]), shape, tolerance = 0.2)
  expect_equal(unname(fit$estimate["rate"]), rate, tolerance = 0.2)

  # Check that the log-likelihood is not NA or -Inf
  expect_false(is.na(fit$loglik))
  expect_false(is.infinite(fit$loglik))

  # Check that the AIC and BIC are calculated
  expect_false(is.na(fit$aic))
  expect_false(is.na(fit$bic))
})

test_that("fitdistdoublecens handles errors correctly", {
  # Test with invalid input
  expect_error(
    fitdistdoublecens(
      data.frame(x = 1:10), # Missing 'left' and 'right' columns
      distr = "gamma"
    ),
    "censdata must contain 'left' and 'right' columns"
  )

  # Test with non-existent distribution
  expect_error(
    fitdistdoublecens(
      data.frame(left = 1:10, right = 2:11),
      distr = "nonexistent_dist"
    )
  )
})

test_that("fitdistdoublecens works with different distributions", {
  set.seed(123)
  n <- 1000

  # Test with normal distribution
  true_mean <- 5
  true_sd <- 2
  samples <- rprimarycensoreddist(
    n, rnorm,
    mean = true_mean, sd = true_sd,
    pwindow = 2, swindow = 2, D = 10
  )

  delay_data <- data.frame(
    left = samples,
    right = samples + 2
  )

  fit_norm <- fitdistdoublecens(
    delay_data,
    distr = "norm",
    start = list(mean = 0, sd = 1),
    D = 10, pwindow = 2
  )

  expect_s3_class(fit_norm, "fitdist")
  expect_equal(unname(fit_norm$estimate["mean"]), true_mean, tolerance = 0.2)
  expect_equal(unname(fit_norm$estimate["sd"]), true_sd, tolerance = 0.2)
})
