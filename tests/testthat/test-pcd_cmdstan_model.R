skip_on_cran()
if (on_ci()) {
  skip_on_os("windows")
  skip_on_os("mac")
}
skip_if_not_installed("cmdstanr")

test_that("pcd_cmdstan_model creates a valid CmdStanModel object", {
  model <- pcd_cmdstan_model()
  expect_s3_class(model, "CmdStanModel")
})

test_that("pcd_cmdstan_model handles custom include paths", {
  custom_path <- tempdir()
  model <- pcd_cmdstan_model(include_paths = custom_path)
  expect_true(custom_path %in% model$include_paths())
})

test_that("pcd_cmdstan_model recovers true values for simple lognormal data", {
  # Simulate data
  set.seed(123)
  n <- 2000
  true_meanlog <- 1.5
  true_sdlog <- 0.5

  simulated_delays <- rprimarycensoreddist(
    n = n,
    rdist = rlnorm,
    meanlog = true_meanlog,
    sdlog = true_sdlog,
    pwindow = 1,
    D = 10
  )

  simulated_data <- data.frame(
    delay = simulated_delays,
    delay_upper = simulated_delays + 1,
    n = 1,
    pwindow = 1,
    relative_obs_time = 10
  )

  # Prepare data for Stan
  stan_data <- pcd_as_cmdstan_data(
    simulated_data,
    dist_id = 1, # Lognormal
    primary_dist_id = 1, # Uniform
    param_bounds = list(lower = c(-Inf, 0), upper = c(Inf, Inf)),
    primary_param_bounds = list(lower = numeric(0), upper = numeric(0)),
    priors = list(location = c(0, 1), scale = c(1, 1)),
    primary_priors = list(location = numeric(0), scale = numeric(0))
  )

  # Fit model
  model <- pcd_cmdstan_model()
  fit <- model$sample(
    data = stan_data,
    seed = 123,
    chains = 4,
    parallel_chains = 4,
    refresh = 0
  )

  # Extract posterior
  posterior <- fit$draws(c("params[1]", "params[2]"), format = "df")

  # Check mean estimates
  expect_equal(mean(posterior$`params[1]`), true_meanlog, tolerance = 0.1)
  expect_equal(mean(posterior$`params[2]`), true_sdlog, tolerance = 0.1)

  # Check credible intervals
  ci_meanlog <- quantile(posterior$`params[1]`, c(0.05, 0.95))
  ci_sdlog <- quantile(posterior$`params[2]`, c(0.05, 0.95))

  expect_gt(true_meanlog, ci_meanlog[1])
  expect_lt(true_meanlog, ci_meanlog[2])
  expect_gt(true_sdlog, ci_sdlog[1])
  expect_lt(true_sdlog, ci_sdlog[2])
})

test_that(
  "pcd_cmdstan_model recovers true values for complex gamma data with
   censoring",
  {
    # Simulate data
    set.seed(456)
    n <- 2000
    true_shape <- 2
    true_rate <- 0.5

    simulated_delays <- rprimarycensoreddist(
      n = n,
      rdist = rgamma,
      shape = true_shape,
      rate = true_rate,
      pwindow = 3,
      D = 8,
      rprimary = rexpgrowth,
      rprimary_args = list(r = 0.1)
    )

    simulated_data <- data.frame(
      delay = simulated_delays,
      delay_upper = simulated_delays + 1,
      n = 1,
      pwindow = 3,
      relative_obs_time = 8
    )

    # Prepare data for Stan
    stan_data <- pcd_as_cmdstan_data(
      simulated_data,
      dist_id = 2, # Gamma
      primary_dist_id = 2, # Exponential growth
      param_bounds = list(lower = c(0, 0), upper = c(Inf, Inf)),
      primary_param_bounds = list(lower = 0, upper = Inf),
      priors = list(location = c(1, 1), scale = c(1, 1)),
      primary_priors = list(location = 0.1, scale = 0.1)
    )

    # Fit model
    model <- pcd_cmdstan_model()
    fit <- model$sample(
      data = stan_data,
      seed = 456,
      chains = 4,
      parallel_chains = 4,
      refresh = 0
    )

    # Extract posterior
    posterior <- fit$draws(
      c("params[1]", "params[2]", "primary_params[1]"),
      format = "df"
    )

    # Check mean estimates
    expect_equal(mean(posterior$`params[1]`), true_shape, tolerance = 0.2)
    expect_equal(mean(posterior$`params[2]`), true_rate, tolerance = 0.1)
    expect_equal(mean(posterior$`primary_params[1]`), 0.1, tolerance = 0.05)

    # Check credible intervals
    ci_shape <- quantile(posterior$`params[1]`, c(0.05, 0.95))
    ci_rate <- quantile(posterior$`params[2]`, c(0.05, 0.95))
    ci_growth <- quantile(posterior$`primary_params[1]`, c(0.05, 0.95))

    expect_gt(true_shape, ci_shape[1])
    expect_lt(true_shape, ci_shape[2])
    expect_gt(true_rate, ci_rate[1])
    expect_lt(true_rate, ci_rate[2])
    expect_gt(0.1, ci_growth[1])
    expect_lt(0.1, ci_growth[2])
  }
)

test_that("pcd_cmdstan_model works with and without within-chain parallelization", {
  # Simulate simple data
  set.seed(789)
  n <- 2000
  simulated_delays <- rprimarycensoreddist(
    n = n,
    rdist = rlnorm,
    meanlog = 1,
    sdlog = 0.5,
    pwindow = 1,
    D = 8
  )

  simulated_data <- data.frame(
    delay = simulated_delays,
    delay_upper = simulated_delays + 1,
    n = 1,
    pwindow = 1,
    relative_obs_time = 8
  )

  # Prepare data for Stan
  stan_data <- pcd_as_cmdstan_data(
    simulated_data,
    dist_id = 1, # Lognormal
    primary_dist_id = 1, # Uniform
    param_bounds = list(lower = c(-Inf, 0), upper = c(Inf, Inf)),
    primary_param_bounds = list(lower = numeric(0), upper = numeric(0)),
    priors = list(location = c(0, 1), scale = c(1, 1)),
    primary_priors = list(location = numeric(0), scale = numeric(0)),
    use_reduce_sum = FALSE
  )

  # Fit model without within-chain parallelization
  model_no_parallel <- pcd_cmdstan_model()
  fit_no_parallel <- model_no_parallel$sample(
    data = stan_data,
    seed = 789,
    chains = 2,
    parallel_chains = 2,
    refresh = 0
  )

  # Prepare data for Stan with within-chain parallelization
  stan_data_parallel <- pcd_as_cmdstan_data(
    simulated_data,
    dist_id = 1, # Lognormal
    primary_dist_id = 1, # Uniform
    param_bounds = list(lower = c(-Inf, 0), upper = c(Inf, Inf)),
    primary_param_bounds = list(lower = numeric(0), upper = numeric(0)),
    priors = list(location = c(0, 1), scale = c(1, 1)),
    primary_priors = list(location = numeric(0), scale = numeric(0)),
    use_reduce_sum = TRUE
  )

  # Fit model with within-chain parallelization
  model_parallel <- pcd_cmdstan_model(cpp_options = list(stan_threads = TRUE))
  fit_parallel <- model_parallel$sample(
    data = stan_data_parallel,
    seed = 789,
    chains = 2,
    parallel_chains = 2,
    threads_per_chain = 2,
    refresh = 0
  )

  # Check that both fits produce similar results
  summary_no_parallel <- fit_no_parallel$summary()
  summary_parallel <- fit_parallel$summary()

  expect_equal(summary_no_parallel$mean, summary_parallel$mean, tolerance = 0.1)
  expect_equal(summary_no_parallel$sd, summary_parallel$sd, tolerance = 0.1)
})


test_that("pcd_cmdstan_model recovers parameters with swindow and pwindow of 2", {
  # Set seed for reproducibility
  set.seed(123)

  # Generate simulated data with swindow and pwindow of 2
  n_obs <- 1000
  true_meanlog <- 1.5
  true_sdlog <- 0.5
  pwindow <- 2
  swindow <- 2
  D <- 30

  simulated_data <- data.frame(
    delay = rprimarycensoreddist(
      n = n_obs,
      rdist = rlnorm,
      meanlog = true_meanlog,
      sdlog = true_sdlog,
      pwindow = pwindow,
      swindow = swindow,
      D = D
    ),
    n = 1
  )
  simulated_data$delay_upper <- simulated_data$delay + swindow
  simulated_data$pwindow <- pwindow
  simulated_data$relative_obs_time <- D

  # Prepare data for Stan
  stan_data <- pcd_as_cmdstan_data(
    simulated_data,
    dist_id = 1, # Lognormal
    primary_dist_id = 1, # Uniform
    param_bounds = list(lower = c(-Inf, 0), upper = c(Inf, Inf)),
    primary_param_bounds = list(lower = numeric(0), upper = numeric(0)),
    priors = list(location = c(0, 1), scale = c(1, 1)),
    primary_priors = list(location = numeric(0), scale = numeric(0))
  )

  # Fit model
  model <- pcd_cmdstan_model()
  fit <- model$sample(
    data = stan_data,
    seed = 456,
    chains = 4,
    parallel_chains = 4,
    refresh = 0
  )

  # Extract posterior summaries
  summary <- fit$summary()

  # Check if true parameters are within 95% credible intervals
  expect_lt(summary$q2.5[summary$variable == "params[1]"], true_meanlog)
  expect_gt(summary$q97.5[summary$variable == "params[1]"], true_meanlog)
  expect_lt(summary$q2.5[summary$variable == "params[2]"], true_sdlog)
  expect_gt(summary$q97.5[summary$variable == "params[2]"], true_sdlog)

  # Check if posterior means are close to true values
  expect_equal(
    summary$mean[summary$variable == "params[1]"],
    true_meanlog,
    tolerance = 0.1
  )
  expect_equal(
    summary$mean[summary$variable == "params[2]"],
    true_sdlog,
    tolerance = 0.1
  )
})
