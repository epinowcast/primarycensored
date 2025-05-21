# fmt: skip file
# Skip tests if fitdistrplus is not installed
if (!requireNamespace("fitdistrplus", quietly = TRUE)) {
  skip("Package 'fitdistrplus' is required for these tests")
}

test_that("fitdistdoublecens works correctly with column names", {
  # Set seed for reproducibility
  set.seed(123)

  # Define true distribution parameters
  n <- 1000
  shape <- 1.77
  rate <- 0.44

  # Generate samples
  samples <- rprimarycensored(
    n,
    rgamma,
    shape = shape,
    rate = rate,
    pwindow = 1,
    swindow = 1,
    D = 8
  )

  # Create data frame with column names
  delay_data <- data.frame(
    left = samples,
    right = samples + 1,
    pwindow = rep(1, n),
    D = rep(8, n)
  )

  # Fit the model using fitdistdoublecens
  fit <- fitdistdoublecens(
    delay_data,
    distr = "gamma",
    start = list(shape = 1, rate = 1)
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

test_that("fitdistdoublecens works with deprecated numeric inputs", {
  # Set seed for reproducibility
  set.seed(123)

  # Define true distribution parameters
  n <- 1000
  shape <- 1.77
  rate <- 0.44

  # Generate samples
  samples <- rprimarycensored(
    n,
    rgamma,
    shape = shape,
    rate = rate,
    pwindow = 1,
    swindow = 1,
    D = 8
  )

  # Create data frame without pwindow and D columns
  delay_data <- data.frame(
    left = samples,
    right = samples + 1
  )

  # Test with deprecated numeric inputs for pwindow and D
  suppressWarnings(expect_warning(
    fit <- fitdistdoublecens( # nolint
      # nolint
      delay_data,
      distr = "gamma",
      start = list(shape = 1, rate = 1),
      pwindow = 1,
      D = 8
    )
  ))

  # Check that the function returns a fitdist object
  expect_s3_class(fit, "fitdist")

  # Check that the estimated parameters are close to the true values
  expect_equal(unname(fit$estimate["shape"]), shape, tolerance = 0.2)
  expect_equal(unname(fit$estimate["rate"]), rate, tolerance = 0.2)
})

test_that("fitdistdoublecens handles errors correctly", {
  # Test with missing columns
  expect_error(
    fitdistdoublecens(
      data.frame(x = 1:10), # Missing required columns
      distr = "gamma"
    ),
    "Missing required columns"
  )

  # Test with non-existent distribution
  expect_error(
    fitdistdoublecens(
      data.frame(
        left = 1:10,
        right = 2:11,
        pwindow = rep(1, 10),
        D = rep(10, 10)
      ),
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
  samples <- rprimarycensored(
    n,
    rnorm,
    mean = true_mean,
    sd = true_sd,
    pwindow = 2,
    swindow = 2,
    D = 10
  )

  delay_data <- data.frame(
    left = samples,
    right = samples + 2,
    pwindow = rep(2, n),
    D = rep(10, n)
  )

  fit_norm <- fitdistdoublecens(
    delay_data,
    distr = "norm",
    start = list(mean = 0, sd = 1)
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
  generate_sample <- function(pwindow, swindow, obs_time) {
    rpcens(
      1,
      rgamma,
      shape = true_shape,
      rate = true_rate,
      pwindow = pwindow,
      swindow = swindow,
      D = obs_time
    )
  }

  pwindows <- rep(1, n)
  swindows <- sample(c(1, 2), n, replace = TRUE)
  obs_times <- rep(10, n)

  samples <- mapply(generate_sample, pwindows, swindows, obs_times)

  delay_data <- data.frame(
    left = samples,
    right = samples + swindows,
    pwindow = pwindows,
    D = obs_times
  )

  fit_gamma <- fitdistdoublecens(
    delay_data,
    distr = "gamma",
    start = list(shape = 2, rate = 1)
  )

  expect_s3_class(fit_gamma, "fitdist")
  expect_equal(unname(fit_gamma$estimate["shape"]), true_shape, tolerance = 0.3)
  expect_equal(unname(fit_gamma$estimate["rate"]), true_rate, tolerance = 0.2)
})

test_that("fitdistdoublecens works with mixed D and primary windows", {
  set.seed(789)
  n <- 1000

  # True parameters for gamma distribution
  true_shape <- 2.5
  true_rate <- 0.6

  # Generate samples with mixed D and primary windows
  generate_sample <- function(pwindow, swindow, obs_time) {
    rpcens(
      1,
      rgamma,
      shape = true_shape,
      rate = true_rate,
      pwindow = pwindow,
      swindow = swindow,
      D = obs_time
    )
  }

  # Create mixed pwindows and D values
  pwindows <- sample(c(1, 2), n, replace = TRUE)
  swindows <- rep(1, n)
  obs_times <- sample(c(8, 12), n, replace = TRUE)

  samples <- mapply(generate_sample, pwindows, swindows, obs_times)

  delay_data <- data.frame(
    left = samples,
    right = samples + swindows,
    pwindow = pwindows,
    D = obs_times
  )

  fit_gamma <- fitdistdoublecens(
    delay_data,
    distr = "gamma",
    start = list(shape = 2, rate = 1)
  )

  expect_s3_class(fit_gamma, "fitdist")
  expect_equal(unname(fit_gamma$estimate["shape"]), true_shape, tolerance = 0.3)
  expect_equal(unname(fit_gamma$estimate["rate"]), true_rate, tolerance = 0.3)
})

test_that("fitdistdoublecens works with custom column names", {
  set.seed(123)
  n <- 1000
  shape <- 1.77
  rate <- 0.44

  samples <- rprimarycensored(
    n,
    rgamma,
    shape = shape,
    rate = rate,
    pwindow = 1,
    swindow = 1,
    D = 8
  )

  # Create data frame with custom column names
  delay_data <- data.frame(
    lower_bound = samples,
    upper_bound = samples + 1,
    primary_window = rep(1, n),
    truncation_time = rep(8, n)
  )

  fit <- fitdistdoublecens(
    delay_data,
    distr = "gamma",
    start = list(shape = 1, rate = 1),
    left = "lower_bound",
    right = "upper_bound",
    pwindow = "primary_window",
    D = "truncation_time"
  )

  expect_s3_class(fit, "fitdist")
  expect_equal(unname(fit$estimate["shape"]), shape, tolerance = 0.2)
  expect_equal(unname(fit$estimate["rate"]), rate, tolerance = 0.2)
})

test_that("fitdistdoublecens handles truncation_check_multiplier correctly", {
  set.seed(123)
  n <- 100
  shape <- 1.77
  rate <- 0.44

  samples <- rprimarycensored(
    n,
    rgamma,
    shape = shape,
    rate = rate,
    pwindow = 1,
    swindow = 1,
    D = 100 # Very large D
  )

  delay_data <- data.frame(
    left = samples,
    right = samples + 1,
    pwindow = rep(1, n),
    D = rep(100, n)
  )

  # Should show a message about large D
  expect_message(
    fitdistdoublecens(
      delay_data,
      distr = "gamma",
      start = list(shape = 1, rate = 1),
      truncation_check_multiplier = 2
    ),
    "truncation time D"
  )

  # Should not show a message when check is disabled
  expect_no_message(
    fitdistdoublecens(
      delay_data,
      distr = "gamma",
      start = list(shape = 1, rate = 1),
      truncation_check_multiplier = NULL
    )
  )
})

test_that(
  "fitdistdoublecens throws error when required packages are not installed", {
  # Create dummy data
  dummy_data <- data.frame(
    left = 1:5,
    right = 2:6,
    pwindow = rep(1, 5),
    D = rep(10, 5)
  )

  # Test for fitdistrplus
  with_mocked_bindings(
    expect_error(
      fitdistdoublecens(dummy_data, "norm"),
      "Package 'fitdistrplus' is required but not installed for this",
      fixed = TRUE
    ),
    requireNamespace = function(pkg, ...) {
      if (pkg == "fitdistrplus") {
        return(FALSE)
      }
      TRUE
    },
    .package = "base"
  )

  # Test for withr
  with_mocked_bindings(
    expect_error(
      fitdistdoublecens(dummy_data, "norm"),
      "Package 'withr' is required but not installed for this function.",
      fixed = TRUE
    ),
    requireNamespace = function(pkg, ...) {
      if (pkg == "withr") {
        return(FALSE)
      }
      TRUE
    },
    .package = "base"
  )

  # Test when both packages are missing
  with_mocked_bindings(
    expect_error(
      fitdistdoublecens(dummy_data, "norm"),
      "Package 'fitdistrplus' is required but not installed",
      fixed = TRUE
    ),
    requireNamespace = function(...) FALSE,
    .package = "base"
  )
})
