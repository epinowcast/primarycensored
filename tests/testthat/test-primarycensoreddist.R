# Test the interactions between dprimarycensoreddist, pprimarycensoreddist,
# and the random number generators for primary events

test_that(
  "rprimarycensoreddist is consistent with dprimarycensoreddist and pprimarycensoreddist", { # nolint
  n <- 10000
  pwindow <- 4
  D <- 10
  samples <- rpcens(
    n, rlnorm, pwindow,
    D = D, meanlog = 0, sdlog = 1
  )

  # Check empirical mean and pmf
  empirical_pmf <- as.vector(table(samples) / n)
  empirical_mean <- mean(samples)
  expect_equal(empirical_mean, 2.9, tolerance = 0.05)
  empirical_sd <- sd(samples)
  expect_equal(empirical_sd, 1.7, tolerance = 0.05)

  # Check empirical cdf against theoretical cdf
  x_values <- 0:(D - 1)
  pmf <- dpcens(x_values, plnorm, pwindow, D = D, meanlog = 0, sdlog = 1)
  theoretical_mean <- sum(x_values * pmf)

  expect_equal(empirical_mean, theoretical_mean, tolerance = 0.05)
  expect_equal(empirical_pmf, pmf, tolerance = 0.05)

  # Check empirical cdf against theoretical cdf
  empirical_cdf <- ecdf(samples)(x_values)
  theoretical_cdf <- ppcens(
    c(x_values[-1], D), plnorm, pwindow, D,
    meanlog = 0, sdlog = 1
  )
  expect_equal(empirical_cdf, theoretical_cdf, tolerance = 0.05)
})


test_that(
  "rprimarycensoreddist is consistent with dprimarycensoreddist and pprimarycensoreddist for exponential growth primary distribution", { # nolint
  n <- 10000
  pwindow <- 3
  D <- 10
  r <- 0.5
  samples <- rpcens(
    n, rlnorm, pwindow,
    D = D,
    rprimary = rexpgrowth,
    rprimary_args = list(r = r),
    meanlog = 1, sdlog = 0.5
  )

  # Check empirical mean and pmf
  empirical_pmf <- as.vector(table(samples) / n)
  empirical_mean <- mean(samples)
  empirical_sd <- sd(samples)

  expect_equal(empirical_mean, 4.3, tolerance = 0.05)
  expect_equal(empirical_sd, 1.6, tolerance = 0.05)

  # Check empirical cdf against theoretical cdf
  x_values <- 0:(D - 1)
  pmf <- dpcens(
    x_values, plnorm, pwindow, D = D,
    dprimary = dexpgrowth,
    dprimary_args = list(r = r),
    meanlog = 1, sdlog = 0.5
  )
  theoretical_mean <- sum(x_values * pmf)

  expect_equal(empirical_mean, theoretical_mean, tolerance = 0.05)
  expect_equal(empirical_pmf, pmf, tolerance = 0.05)

  # Check empirical cdf against theoretical cdf
  empirical_cdf <- ecdf(samples)(x_values)
  theoretical_cdf <- ppcens(
    c(x_values[-1], D), plnorm, pwindow, D,
    dprimary = dexpgrowth,
    dprimary_args = list(r = r),
    meanlog = 1, sdlog = 0.5
  )
  expect_equal(empirical_cdf, theoretical_cdf, tolerance = 0.05)
})

test_that(
  "Exponential distribution with daily censoring matches analytical solution", {
    # Parameters
    rate <- 1 # Exponential rate parameter
    pwindow <- 1 # Primary event window (1 day)
    swindow <- 1 # Secondary event window (1 day)
    D <- 20 # Maximum delay to consider

    # Analytical solution
    analytical_pmf <- function(s) {
      (1 - exp(-1)) * (exp(1) - 1) * exp(-s)
    }
    analytical_pmf <- function(t, D) {
      (exp(-(t - 1)) - exp(-t)) / (1 - exp(-D))
    }

            expected_pmf = [(exp(-(t - 1)) - exp(-t)) / (1 - exp(-5)) for t in 1:5]
        pmf = censored_pmf(dist,
            Val(:single_censored);
            primary_approximation_point = 0.0,
            Δd = 1.0,
            D = 5.0)

    # Numerical solution using our functions
    numerical_pmf <- dpcens(
      0:(D - 1), dexp, pwindow, swindow, 20,
      rate = rate
    )
    numerical_cdf <- ppcens(
      0:(D - 1), dexp, pwindow, 20,
      rate = rate
    )

    # Compare PMF
    analytical_values <- sapply(0:(D - 1), analytical_pmf, D = D)
    expect_equal(numerical_pmf, analytical_values, tolerance = 1e-6)

    # Compare CDF
    analytical_cdf <- cumsum(analytical_values)
    expect_equal(numerical_cdf, analytical_cdf, tolerance = 1e-6)

    # Check that PMF sums to 1 (approximately)
    expect_equal(sum(numerical_pmf), 1, tolerance = 1e-6)

    # Check that CDF approaches 1
    expect_equal(numerical_cdf[length(numerical_cdf)], 1, tolerance = 1e-6)
  }
)

test_that(
  "rprimarycensoreddist matches analytical solution for exponential
   distribution",
  {
    # Parameters
    n <- 100000
    rate <- 1 # Exponential rate parameter
    pwindow <- 1 # Primary event window (1 day)
    swindow <- 1 # Secondary event window (1 day)
    D <- 20 # Maximum delay to consider

    # Generate samples
    samples <- rpcens(n, rexp, pwindow, swindow, D, rate = rate)

    # Analytical solution
    analytical_pmf <- function(s) {
      if (s == 0) {
        exp(-1)
      } else {
        (1 - exp(-1)) * (exp(1) - 1) * exp(-s)
      }
    }

    # Compare empirical distribution to analytical solution
    empirical_pmf <- table(samples) / n
    analytical_values <- sapply(
      as.numeric(names(empirical_pmf)), analytical_pmf
    )

    expect_equal(as.numeric(empirical_pmf), analytical_values, tolerance = 0.01)
  }
)
