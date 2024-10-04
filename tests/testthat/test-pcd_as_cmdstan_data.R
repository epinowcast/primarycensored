test_that("pcd_as_stan_data correctly formats data", {
  data <- data.frame(
    delay = c(1, 2, 3),
    delay_upper = c(2, 3, 4),
    n = c(10, 20, 15),
    pwindow = c(1, 1, 2),
    relative_obs_time = c(10, 10, 10)
  )

  dist_id <- 1
  primary_dist_id <- 1
  param_bounds <- list(lower = c(0, 0), upper = c(10, 10))
  primary_param_bounds <- list(lower = numeric(0), upper = numeric(0))
  priors <- list(location = c(1, 1), scale = c(1, 1))
  primary_priors <- list(location = numeric(0), scale = numeric(0))

  result <- pcd_as_stan_data(
    data,
    dist_id = dist_id,
    primary_dist_id = primary_dist_id,
    param_bounds = param_bounds,
    primary_param_bounds = primary_param_bounds,
    priors = priors,
    primary_priors = primary_priors
  )

  expect_type(result, "list")
  expect_identical(result$N, nrow(data))
  expect_identical(result$d, data$delay)
  expect_identical(result$d_upper, data$delay_upper)
  expect_identical(result$n, data$n)
  expect_identical(result$pwindow, data$pwindow)
  expect_identical(result$D, data$relative_obs_time)
  expect_identical(result$dist_id, dist_id)
  expect_identical(result$primary_dist_id, primary_dist_id)
  expect_identical(result$n_params, length(param_bounds$lower))
  expect_identical(result$n_primary_params, length(primary_param_bounds$lower))
  expect_identical(result$compute_log_lik, 0L)
  expect_identical(result$use_reduce_sum, 0L)
  expect_identical(result$param_lower_bounds, param_bounds$lower)
  expect_identical(result$param_upper_bounds, param_bounds$upper)
  expect_identical(
    result$primary_param_lower_bounds, primary_param_bounds$lower
  )
  expect_identical(
    result$primary_param_upper_bounds, primary_param_bounds$upper
  )
  expect_identical(result$prior_location, priors$location)
  expect_identical(result$prior_scale, priors$scale)
  expect_identical(result$primary_prior_location, primary_priors$location)
  expect_identical(result$primary_prior_scale, primary_priors$scale)
})

test_that("pcd_as_stan_data handles missing columns correctly", {
  data <- data.frame(
    delay = c(1, 2, 3),
    n = c(10, 20, 15),
    pwindow = c(1, 1, 2)
  )

  expect_error(
    pcd_as_stan_data(
      data,
      dist_id = 1,
      primary_dist_id = 1,
      param_bounds = list(lower = c(0, 0), upper = c(10, 10)),
      primary_param_bounds = list(lower = numeric(0), upper = numeric(0)),
      priors = list(location = c(1, 1), scale = c(1, 1)),
      primary_priors = list(location = numeric(0), scale = numeric(0))
    ),
    "Missing required columns: delay_upper, relative_obs_time"
  )
})

test_that("pcd_as_stan_data handles optional parameters correctly", {
  data <- data.frame(
    delay = c(1, 2, 3),
    delay_upper = c(2, 3, 4),
    n = c(10, 20, 15),
    pwindow = c(1, 1, 2),
    relative_obs_time = c(10, 10, 10)
  )

  result <- pcd_as_stan_data(
    data,
    dist_id = 1,
    primary_dist_id = 1,
    param_bounds = list(lower = c(0, 0), upper = c(10, 10)),
    primary_param_bounds = list(lower = numeric(0), upper = numeric(0)),
    priors = list(location = c(1, 1), scale = c(1, 1)),
    primary_priors = list(location = numeric(0), scale = numeric(0)),
    compute_log_lik = TRUE,
    use_reduce_sum = TRUE
  )

  expect_identical(result$compute_log_lik, 1L)
  expect_identical(result$use_reduce_sum, 1L)
})

test_that("pcd_as_stan_data handles custom column names correctly", {
  data <- data.frame(
    obs_delay = c(1, 2, 3),
    obs_delay_upper = c(2, 3, 4),
    count = c(10, 20, 15),
    primary_window = c(1, 1, 2),
    obs_time = c(10, 10, 10)
  )

  result <- pcd_as_stan_data(
    data,
    delay = "obs_delay",
    delay_upper = "obs_delay_upper",
    n = "count",
    pwindow = "primary_window",
    relative_obs_time = "obs_time",
    dist_id = 1,
    primary_dist_id = 1,
    param_bounds = list(lower = c(0, 0), upper = c(10, 10)),
    primary_param_bounds = list(lower = numeric(0), upper = numeric(0)),
    priors = list(location = c(1, 1), scale = c(1, 1)),
    primary_priors = list(location = numeric(0), scale = numeric(0))
  )

  expect_identical(result$d, data$obs_delay)
  expect_identical(result$d_upper, data$obs_delay_upper)
  expect_identical(result$n, data$count)
  expect_identical(result$pwindow, data$primary_window)
  expect_identical(result$D, data$obs_time)
})
