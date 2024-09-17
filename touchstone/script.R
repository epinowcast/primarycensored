# see `help(run_script, package = 'touchstone')` on how to run this
# interactively

# installs branches to benchmark
touchstone::branch_install()

# Benchmark for pprimarycensoreddist with lognormal distribution
touchstone::benchmark_run(
  expr_before_benchmark = {
    library(primarycensoreddist)
    q <- seq(0, 10, by = 0.1)
  },
  pprimarycensoreddist_lnorm = {
    pprimarycensoreddist(q, plnorm, meanlog = 0, sdlog = 1)
  },
  n = 100
)

# Benchmark for pprimarycensoreddist with exponential growth
touchstone::benchmark_run(
  expr_before_benchmark = {
    library(primarycensoreddist)
    q <- seq(0, 10, by = 0.1)
  },
  pprimarycensoreddist_expgrowth = {
    pprimarycensoreddist(
      q, plnorm,
      dprimary = dexpgrowth,
      dprimary_args = list(r = 0.2),
      meanlog = 0, sdlog = 1
    )
  },
  n = 100
)

# Benchmark for dprimarycensoreddist with Weibull distribution
touchstone::benchmark_run(
  expr_before_benchmark = {
    library(primarycensoreddist)
    x <- seq(0, 10, by = 0.1)
  },
  dprimarycensoreddist_weibull = {
    dprimarycensoreddist(x, pweibull, shape = 1.5, scale = 2.0)
  },
  n = 100
)

# Benchmark for dprimarycensoreddist with exponential growth
touchstone::benchmark_run(
  expr_before_benchmark = {
    library(primarycensoreddist)
    x <- seq(0, 10, by = 0.1)
  },
  dprimarycensoreddist_expgrowth = {
    dprimarycensoreddist(
      x, pweibull,
      dprimary = dexpgrowth,
      dprimary_args = list(r = 0.2),
      shape = 1.5, scale = 2.0
    )
  },
  n = 100
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
      rprimary = rexpgrowth,
      rprimary_args = list(r = 0.1)
    )
    delay_data <- data.frame(
      left = samples,
      right = samples + swindow
    )
  },
  fitdistdoublecens_expgrowth = {
    fitdistdoublecens(
      delay_data,
      distr = "gamma",
      start = list(shape = 1, rate = 1),
      D = D,
      pwindow = pwindow,
      dprimary = dexpgrowth,
      dprimary_args = list(r = 0.1)
    )
  },
  n = 10
)


# Benchmark for fitting lognormal distribution
touchstone::benchmark_run(
  expr_before_benchmark = {
    source("touchstone/setup.R")
  },
  cmdstan_fit_lognormal = {
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
    source("touchstone/setup.R")
  },
  cmdstan_fit_gamma = {
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
