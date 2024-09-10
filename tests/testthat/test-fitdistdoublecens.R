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
    left = samples
  )
  delay_data$right <- delay_data$left + 1

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
    left = samples
  )
  delay_data$right <- delay_data$left + 2

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

test_that("fitdistdoublecens works with mixed secondary windows", {
  set.seed(456)
  n <- 1000

  # True parameters for gamma distribution
  true_shape <- 3
  true_rate <- 0.5

  # Generate samples with mixed secondary windows
  # Generate samples with mixed secondary windows
  generate_sample <- function(pwindow, swindow, obs_time) {
    rpcens(
      1, rgamma,
      shape = true_shape, rate = true_rate,
      pwindow = pwindow, swindow = swindow, D = obs_time
    )
  }

  pwindows <- rep(1, n)
  swindows <- sample(c(1, 2), n, replace = TRUE)
  obs_times <- rep(10, n)

  samples <- mapply(generate_sample, pwindows, swindows, obs_times)

  delay_data <- data.frame(
    left = samples
  )
  delay_data$right <- delay_data$left + swindows

  fit_gamma <- fitdistdoublecens(
    delay_data,
    distr = "gamma",
    start = list(shape = 2, rate = 1),
    D = 10, pwindow = 1
  )

  expect_s3_class(fit_gamma, "fitdist")
  expect_equal(unname(fit_gamma$estimate["shape"]), true_shape, tolerance = 0.3)
  expect_equal(unname(fit_gamma$estimate["rate"]), true_rate, tolerance = 0.2)
})

test_that("fitdistdoublecens throws error when fitdistrplus is not installed", {
  with_mocked_bindings(
    {
      # Create dummy data
      dummy_data <- data.frame(left = 1:5, right = 2:6)

      # Expect an error when trying to use fitdistdoublecens
      expect_error(
        fitdistdoublecens(dummy_data, "norm"),
        "Package 'fitdistrplus' is required but not installed for this"
      )
    },
    requireNamespace = function(...) FALSE,
    .package = "base"
  )
})
