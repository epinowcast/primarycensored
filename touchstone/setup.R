library(primarycensoreddist)
library(dplyr)
library(cmdstanr)

# Set cmdstan path
cmdstanr::set_cmdstan_path()

# Use 2 cores
options(mc.cores = 2)

# Example 1: Lognormal distribution with uniform primary distribution
set.seed(123)
n1 <- 2000
true_meanlog1 <- 1.5
true_sdlog1 <- 0.5

simulated_delays1 <- rprimarycensoreddist(
  n = n1,
  rdist = rlnorm,
  meanlog = true_meanlog1,
  sdlog = true_sdlog1,
  pwindow = 1,
  D = 10
)

example_data1 <- data.frame(
  delay = simulated_delays1,
  delay_upper = simulated_delays1 + 1,
  pwindow = 1,
  relative_obs_time = 10
)

delay_counts1 <- example_data1 |>
  summarise(
    n = n(),
    .by = c(pwindow, relative_obs_time, delay, delay_upper)
  )

# Convert to Stan data format for Example 1
stan_data1 <- pcd_as_cmdstan_data(
  delay_counts1,
  dist_id = 1, # Lognormal
  primary_dist_id = 1, # Uniform
  param_bounds = list(lower = c(-Inf, 0), upper = c(Inf, Inf)),
  primary_param_bounds = list(lower = numeric(0), upper = numeric(0)),
  priors = list(location = c(0, 1), scale = c(1, 1)),
  primary_priors = list(location = numeric(0), scale = numeric(0))
)

# Example 2: Gamma distribution with exponential growth primary distribution
set.seed(456)
n2 <- 2000
true_shape2 <- 2
true_rate2 <- 0.5

simulated_delays2 <- rprimarycensoreddist(
  n = n2,
  rdist = rgamma,
  shape = true_shape2,
  rate = true_rate2,
  pwindow = 2,
  swindow = 2,
  D = 8,
  rprimary = rexpgrowth,
  rprimary_args = list(r = 0.1)
)

example_data2 <- data.frame(
  delay = simulated_delays2,
  delay_upper = simulated_delays2 + 2,
  pwindow = 2,
  relative_obs_time = 8
)

delay_counts2 <- example_data2 |>
  summarise(
    n = n(),
    .by = c(pwindow, relative_obs_time, delay, delay_upper)
  )

# Convert to Stan data format for Example 2
stan_data2 <- pcd_as_cmdstan_data(
  delay_counts2,
  dist_id = 2, # Gamma
  primary_dist_id = 2, # Exponential growth
  param_bounds = list(lower = c(0, 0), upper = c(Inf, Inf)),
  primary_param_bounds = list(lower = 0, upper = Inf),
  priors = list(location = c(2, 1), scale = c(0.5, 0.5)),
  primary_priors = list(location = 0.1, scale = 0.1)
)

# Precompile the model
model <- suppressMessages(suppressWarnings(pcd_cmdstan_model()))
