# see `help(run_script, package = 'touchstone')` on how to run this
# interactively

# installs branches to benchmark
touchstone::branch_install()

# Benchmark for pprimarycensoreddist with lognormal distribution
touchstone::benchmark_run(
  expr_before_benchmark = {
    library(primarycensoreddist)
    q <- seq(0, 10, by = 0.01)
  },
  pprimarycensoreddist_lnorm = {
    pprimarycensoreddist(q, plnorm, meanlog = 0, sdlog = 1, D = 12)
  },
  n = 20
)

# Benchmark for pprimarycensoreddist with exponential growth
touchstone::benchmark_run(
  expr_before_benchmark = {
    library(primarycensoreddist)
    q <- seq(0, 10, by = 0.01)
  },
  pprimarycensoreddist_plnrom_expgrowth = {
    pprimarycensoreddist(
      q, plnorm,
      dprimary = dexpgrowth,
      dprimary_args = list(r = 0.2),
      meanlog = 0, sdlog = 1, D = 12
    )
  },
  n = 20
)

# Benchmark for dprimarycensoreddist with Weibull distribution
touchstone::benchmark_run(
  expr_before_benchmark = {
    library(primarycensoreddist)
    x <- seq(0, 10, by = 1)
  },
  dprimarycensoreddist_weibull = {
    dprimarycensoreddist(x, pweibull, shape = 1.5, scale = 2.0, D = 12)
  },
  n = 20
)

# Benchmark for dprimarycensoreddist with exponential growth
touchstone::benchmark_run(
  expr_before_benchmark = {
    library(primarycensoreddist)
    x <- seq(0, 10, by = 1)
  },
  dprimarycensoreddist_pweibull_expgrowth = {
    dprimarycensoreddist(
      x, pweibull,
      dprimary = dexpgrowth,
      dprimary_args = list(r = 0.2),
      shape = 1.5, scale = 2.0, D = 12
    )
  },
  n = 20
)

# Benchmark for fitdistdoublecens with normal distribution
touchstone::benchmark_run(
  expr_before_benchmark = {
    library(primarycensoreddist)
    library(fitdistrplus)
    set.seed(123)
    n <- 1000
    true_mean <- 5
    true_sd <- 2
    pwindow <- 1
    swindow <- 1
    D <- 10
    samples <- rprimarycensoreddist(
      n, rnorm,
      mean = true_mean, sd = true_sd,
      pwindow = pwindow, swindow = swindow, D = D
    )
    delay_data <- data.frame(
      left = samples,
      right = samples + swindow
    )
  },
  fitdistdoublecens_normal = {
    fitdistdoublecens(
      delay_data,
      distr = "norm",
      start = list(mean = 0, sd = 1),
      D = D,
      pwindow = pwindow
    )
  },
  n = 10
)

# Benchmark for fitdistdoublecens with exponential growth
touchstone::benchmark_run(
  expr_before_benchmark = {
    library(primarycensoreddist)
    library(fitdistrplus)
    set.seed(456)
    n <- 1000
    true_shape <- 2
    true_rate <- 0.5
    pwindow <- 1
    swindow <- 1
    D <- 8
    samples <- rprimarycensoreddist(
      n, rgamma,
      shape = true_shape, rate = true_rate,
      pwindow = pwindow, swindow = swindow, D = D,
    )
    delay_data <- data.frame(
      left = samples,
      right = samples + swindow
    )
  },
  fitdistdoublecens_gamma_expgrowth = {
    fitdistdoublecens(
      delay_data,
      distr = "gamma",
      start = list(shape = 1, rate = 1),
      D = D,
      pwindow = pwindow
    )
  },
  n = 10
)

# Benchmark for fitting lognormal distribution
touchstone::benchmark_run(
  expr_before_benchmark = {
    library(primarycensoreddist)
    library(dplyr)
    library(cmdstanr)

    cmdstanr::set_cmdstan_path()
    options(mc.cores = 2)

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

    stan_data1 <- pcd_as_cmdstan_data(
      delay_counts1,
      dist_id = 1,
      primary_dist_id = 1,
      param_bounds = list(lower = c(-Inf, 0), upper = c(Inf, Inf)),
      primary_param_bounds = list(lower = numeric(0), upper = numeric(0)),
      priors = list(location = c(0, 1), scale = c(1, 1)),
      primary_priors = list(location = numeric(0), scale = numeric(0))
    )

    model <- suppressMessages(suppressWarnings(pcd_cmdstan_model()))
  },
  cmdstan_fit_lognormal_unif = {
    fit1 <- model$sample(
      data = stan_data1,
      seed = 123,
      chains = 2,
      parallel_chains = 2,
      refresh = 0,
      show_messages = FALSE,
      iter_warmup = 500,
      iter_sampling = 500
    )
  },
  n = 10
)

# Benchmark for fitting gamma distribution
touchstone::benchmark_run(
  expr_before_benchmark = {
    library(primarycensoreddist)
    library(dplyr)
    library(cmdstanr)

    cmdstanr::set_cmdstan_path()
    options(mc.cores = 2)

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

    stan_data2 <- pcd_as_cmdstan_data(
      delay_counts2,
      dist_id = 2,
      primary_dist_id = 2,
      param_bounds = list(lower = c(0, 0), upper = c(Inf, Inf)),
      primary_param_bounds = list(lower = 0, upper = Inf),
      priors = list(location = c(2, 1), scale = c(0.5, 0.5)),
      primary_priors = list(location = 0.1, scale = 0.1)
    )

    model <- suppressMessages(suppressWarnings(pcd_cmdstan_model()))
  },
  cmdstan_fit_gamma_expgrowth = {
    fit2 <- model$sample(
      data = stan_data2,
      seed = 456,
      chains = 2,
      parallel_chains = 2,
      refresh = 0,
      show_messages = FALSE,
      iter_warmup = 500,
      iter_sampling = 500
    )
  },
  n = 10
)

# create artifacts used downstream in the GitHub Action.
touchstone::benchmark_analyze()
