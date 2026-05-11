test_that("pcd_as_stan_data correctly formats data", {
  data <- data.frame(
    delay = c(1, 2, 3),
    delay_upper = c(2, 3, 4),
    n = c(10, 20, 15),
    pwindow = c(1, 1, 2),
    relative_obs_time = c(10, 10, 10)
  )

  dist_id <- 1
  primary_id <- 1
  param_bounds <- list(lower = c(0, 0), upper = c(10, 10))
  primary_param_bounds <- list(lower = numeric(0), upper = numeric(0))
  priors <- list(location = c(1, 1), scale = c(1, 1))
  primary_priors <- list(location = numeric(0), scale = numeric(0))

  expect_message(
    result <- pcd_as_stan_data( # nolint
      data,
      dist_id = dist_id,
      primary_id = primary_id,
      param_bounds = param_bounds,
      primary_param_bounds = primary_param_bounds,
      priors = priors,
      primary_priors = primary_priors
    )
  )
  # Default (parametric) path: nonparametric off, np_* fields sized 0.
  expect_identical(result$nonparametric, 0L)
  expect_identical(result$K_np, 0L)
  # np_boundaries is declared `vector[K_np + 1]` so length 1 when K_np = 0.
  expect_length(result$np_boundaries, 1)
  expect_length(result$np_dirichlet_alpha, 0)

  expect_type(result, "list")
  expect_identical(result$N, nrow(data))
  expect_identical(result$d, data$delay)
  expect_identical(result$d_upper, data$delay_upper)
  expect_identical(result$n, data$n)
  expect_identical(result$pwindow, data$pwindow)
  expect_identical(result$D, data$relative_obs_time)
  expect_identical(result$dist_id, dist_id)
  expect_identical(result$primary_id, primary_id)
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
      primary_id = 1,
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

  expect_message(
    result <- pcd_as_stan_data( # nolint
      data,
      dist_id = 1,
      primary_id = 1,
      param_bounds = list(lower = c(0, 0), upper = c(10, 10)),
      primary_param_bounds = list(lower = numeric(0), upper = numeric(0)),
      priors = list(location = c(1, 1), scale = c(1, 1)),
      primary_priors = list(location = numeric(0), scale = numeric(0)),
      compute_log_lik = TRUE,
      use_reduce_sum = TRUE
    )
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

  expect_message(
    result <- pcd_as_stan_data( # nolint
      data,
      delay = "obs_delay",
      delay_upper = "obs_delay_upper",
      n = "count",
      pwindow = "primary_window",
      relative_obs_time = "obs_time",
      dist_id = 1,
      primary_id = 1,
      param_bounds = list(lower = c(0, 0), upper = c(10, 10)),
      primary_param_bounds = list(lower = numeric(0), upper = numeric(0)),
      priors = list(location = c(1, 1), scale = c(1, 1)),
      primary_priors = list(location = numeric(0), scale = numeric(0))
    )
  )

  expect_identical(result$d, data$obs_delay)
  expect_identical(result$d_upper, data$obs_delay_upper)
  expect_identical(result$n, data$count)
  expect_identical(result$pwindow, data$primary_window)
  expect_identical(result$D, data$obs_time)
})

test_that("pcd_as_stan_data handles additional columns correctly", {
  # Data with extra columns that should be ignored

  data <- data.frame(
    delay = c(1, 2, 3),
    delay_upper = c(2, 3, 4),
    n = c(10, 20, 15),
    pwindow = c(1, 1, 2),
    relative_obs_time = c(10, 10, 10),
    extra_column = c("a", "b", "c"),
    another_extra = c(100, 200, 300),
    yet_another = c(TRUE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )

  expect_message(
    result <- pcd_as_stan_data( # nolint
      data,
      dist_id = 1,
      primary_id = 1,
      param_bounds = list(lower = c(0, 0), upper = c(10, 10)),
      primary_param_bounds = list(lower = numeric(0), upper = numeric(0)),
      priors = list(location = c(1, 1), scale = c(1, 1)),
      primary_priors = list(location = numeric(0), scale = numeric(0))
    )
  )

  # Verify the result contains correct data from required columns
  expect_type(result, "list")
  expect_identical(result$N, nrow(data))
  expect_identical(result$d, data$delay)
  expect_identical(result$d_upper, data$delay_upper)
  expect_identical(result$n, data$n)
  expect_identical(result$pwindow, data$pwindow)
  expect_identical(result$D, data$relative_obs_time)

  # Verify extra columns are not included in result

  expect_null(result$extra_column)
  expect_null(result$another_extra)
  expect_null(result$yet_another)
})

test_that("pcd_as_stan_data handles additional columns with custom names", {
  # Data with custom column names AND extra columns
  data <- data.frame(
    obs_delay = c(1, 2, 3),
    obs_delay_upper = c(2, 3, 4),
    count = c(10, 20, 15),
    primary_window = c(1, 1, 2),
    obs_time = c(10, 10, 10),
    spurious_data = c(1.5, 2.5, 3.5),
    metadata = c("x", "y", "z"),
    stringsAsFactors = FALSE
  )

  expect_message(
    result <- pcd_as_stan_data( # nolint
      data,
      delay = "obs_delay",
      delay_upper = "obs_delay_upper",
      n = "count",
      pwindow = "primary_window",
      relative_obs_time = "obs_time",
      dist_id = 1,
      primary_id = 1,
      param_bounds = list(lower = c(0, 0), upper = c(10, 10)),
      primary_param_bounds = list(lower = numeric(0), upper = numeric(0)),
      priors = list(location = c(1, 1), scale = c(1, 1)),
      primary_priors = list(location = numeric(0), scale = numeric(0))
    )
  )

  # Verify the result contains correct data from required columns
  expect_identical(result$d, data$obs_delay)
  expect_identical(result$d_upper, data$obs_delay_upper)
  expect_identical(result$n, data$count)
  expect_identical(result$pwindow, data$primary_window)
  expect_identical(result$D, data$obs_time)

  # Verify extra columns are not included in result
  expect_null(result$spurious_data)
  expect_null(result$metadata)
})

test_that("pcd_as_stan_data populates np_* fields for simplex paramtype", {
  data <- data.frame(
    delay = c(1, 2, 3),
    delay_upper = c(2, 3, 4),
    n = c(10, 20, 15),
    pwindow = c(1, 1, 2),
    relative_obs_time = c(10, 10, 10)
  )

  expect_message(
    result <- pcd_as_stan_data( # nolint
      data,
      dist_id = 1, # ignored when nonparametric is supplied
      primary_id = 1,
      param_bounds = list(lower = c(0, 0), upper = c(10, 10)),
      primary_param_bounds = list(lower = numeric(0), upper = numeric(0)),
      priors = list(location = c(1, 1), scale = c(1, 1)),
      primary_priors = list(location = numeric(0), scale = numeric(0)),
      nonparametric = list(
        K = 5,
        boundaries = 0:5,
        paramtype = "simplex"
      )
    )
  )

  expect_identical(result$nonparametric, 1L)
  expect_identical(result$dist_id, 26L)
  expect_identical(result$K_np, 5L)
  expect_identical(result$np_paramtype, 1L)
  expect_identical(result$np_boundaries, as.numeric(0:5))
  expect_identical(result$np_dirichlet_alpha, rep(1, 5))
  expect_identical(result$n_params, 0L)
})

test_that("pcd_as_stan_data populates np_* fields for hazard paramtype", {
  data <- data.frame(
    delay = c(1, 2, 3),
    delay_upper = c(2, 3, 4),
    n = c(10, 20, 15),
    pwindow = c(1, 1, 2),
    relative_obs_time = c(10, 10, 10)
  )

  expect_message(
    result <- pcd_as_stan_data( # nolint
      data,
      dist_id = 1,
      primary_id = 1,
      param_bounds = list(lower = c(0, 0), upper = c(10, 10)),
      primary_param_bounds = list(lower = numeric(0), upper = numeric(0)),
      priors = list(location = c(1, 1), scale = c(1, 1)),
      primary_priors = list(location = numeric(0), scale = numeric(0)),
      nonparametric = list(
        K = 4,
        boundaries = 0:4,
        paramtype = "hazard",
        hazard_priors = list(alpha_mean = -1, alpha_sd = 2)
      )
    )
  )

  expect_identical(result$nonparametric, 1L)
  expect_identical(result$dist_id, 27L)
  expect_identical(result$K_np, 4L)
  expect_identical(result$np_paramtype, 2L)
  expect_identical(result$np_alpha_mean, -1)
  expect_identical(result$np_alpha_sd, 2)
  # Defaults for unspecified hazard_priors entries.
  expect_identical(result$np_log_sigma_mean, 0)
  expect_identical(result$np_log_sigma_sd, 1)
})

test_that("pcd_as_stan_data populates np_* fields for hazardre paramtype", {
  data <- data.frame(
    delay = c(1, 2, 3),
    delay_upper = c(2, 3, 4),
    n = c(10, 20, 15),
    pwindow = c(1, 1, 2),
    relative_obs_time = c(10, 10, 10)
  )

  expect_message(
    result <- pcd_as_stan_data( # nolint
      data,
      dist_id = 1,
      primary_id = 1,
      param_bounds = list(lower = c(0, 0), upper = c(10, 10)),
      primary_param_bounds = list(lower = numeric(0), upper = numeric(0)),
      priors = list(location = c(1, 1), scale = c(1, 1)),
      primary_priors = list(location = numeric(0), scale = numeric(0)),
      nonparametric = list(
        K = 4,
        boundaries = 0:4,
        paramtype = "hazardre"
      )
    )
  )

  # Same dist_id as the RW hazard path; the data flag picks the variant.
  expect_identical(result$dist_id, 27L)
  expect_identical(result$np_paramtype, 3L)
})

test_that("pcd_as_stan_data validates nonparametric input", {
  data <- data.frame(
    delay = c(1, 2, 3),
    delay_upper = c(2, 3, 4),
    n = c(10, 20, 15),
    pwindow = c(1, 1, 2),
    relative_obs_time = c(10, 10, 10)
  )
  base_args <- list(
    data,
    dist_id = 1,
    primary_id = 1,
    param_bounds = list(lower = numeric(0), upper = numeric(0)),
    primary_param_bounds = list(lower = numeric(0), upper = numeric(0)),
    priors = list(location = numeric(0), scale = numeric(0)),
    primary_priors = list(location = numeric(0), scale = numeric(0))
  )
  # Wrong boundaries length.
  expect_error(
    suppressMessages(do.call(pcd_as_stan_data, c(base_args, list(
      nonparametric = list(K = 3, boundaries = 0:2, paramtype = "simplex")
    )))),
    "boundaries.*length K \\+ 1"
  )
  # Unknown paramtype.
  expect_error(
    suppressMessages(do.call(pcd_as_stan_data, c(base_args, list(
      nonparametric = list(K = 3, boundaries = 0:3, paramtype = "wibble")
    ))))
  )
  # Missing K.
  expect_error(
    suppressMessages(do.call(pcd_as_stan_data, c(base_args, list(
      nonparametric = list(boundaries = 0:3, paramtype = "simplex")
    )))),
    "missing required elements"
  )
})
