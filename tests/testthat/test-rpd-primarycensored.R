# Test the interactions between dprimarycensored, pprimarycensored,
# and the random number generators for primary events

test_that(
  "rprimarycensored is consistent with dprimarycensored and
   pprimarycensored",
  { # nolint
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
    expect_equal(empirical_sd, 1.8, tolerance = 0.05)

    # Check empirical cdf against theoretical cdf
    x_values <- 0:(D - 1)
    pmf <- dpcens(x_values, plnorm, pwindow, D = D, meanlog = 0, sdlog = 1)
    theoretical_mean <- sum(x_values * pmf)

    expect_equal(empirical_mean, theoretical_mean, tolerance = 0.05)
    expect_equal(empirical_pmf, pmf, tolerance = 0.05)

    # Check empirical cdf against theoretical cdf
    empirical_cdf <- ecdf(samples)(x_values)
    theoretical_cdf <- ppcens(
      c(x_values[-1], D), plnorm,
      pwindow = pwindow, D = D,
      meanlog = 0, sdlog = 1
    )
    expect_equal(empirical_cdf, theoretical_cdf, tolerance = 0.05)
  }
)


test_that(
  "rprimarycensored is consistent with dprimarycensored and
   pprimarycensored for exponential growth primary distribution",
  { # nolint
    n <- 1e6
    pwindow <- 3
    D <- 6
    r <- 0.5
    samples <- rpcens(
      n, rlnorm, pwindow,
      D = D,
      rprimary = rexpgrowth,
      rprimary_args = list(r = r),
      meanlog = 1.5, sdlog = 0.5
    )

    # Check empirical mean and pmf
    empirical_pmf <- as.vector(table(samples) / n)
    empirical_mean <- mean(samples)
    empirical_sd <- sd(samples)

    expect_equal(empirical_mean, 4.12, tolerance = 0.01)
    expect_equal(empirical_sd, 0.94, tolerance = 0.01)

    # Check empirical cdf against theoretical cdf
    x_values <- 0:(D - 1)
    pmf <- dpcens(
      x_values, plnorm, pwindow,
      D = D,
      dprimary = dexpgrowth,
      dprimary_args = list(r = r),
      meanlog = 1.5, sdlog = 0.5
    )
    theoretical_mean <- sum(x_values * pmf)

    expect_equal(empirical_mean, theoretical_mean, tolerance = 0.01)
    expect_equal(empirical_pmf, pmf, tolerance = 0.01)

    # Check empirical cdf against theoretical cdf
    empirical_cdf <- ecdf(samples)(x_values)
    theoretical_cdf <- ppcens(
      c(x_values[-1], D), plnorm,
      pwindow = pwindow, D = D,
      dprimary = dexpgrowth,
      dprimary_args = list(r = r),
      meanlog = 1.5, sdlog = 0.5
    )
    expect_equal(empirical_cdf, theoretical_cdf, tolerance = 0.01)
  }
)

test_that(
  "rprimarycensored with wider windows and different delay distribution
   mathches p and d numerically",
  {
    n <- 1e6
    pwindow <- 7 # One week primary window
    swindow <- 3 # Three-day secondary window
    D <- 30 # Maximum delay of 30 days

    # Using Weibull distribution for delay
    shape <- 2
    scale <- 10

    samples <- rpcens(
      n, rweibull, pwindow, swindow,
      rprimary = rexpgrowth,
      rprimary_args = list(r = 0.5),
      D = D, shape = shape, scale = scale
    )

    # Check empirical mean and standard deviation
    empirical_mean <- mean(samples)
    empirical_sd <- sd(samples)

    # Calculate theoretical mean and standard deviation
    x_values <- seq(0, D - swindow, by = swindow)
    pmf <- dpcens(
      x_values, pweibull, pwindow, swindow,
      dprimary = dexpgrowth,
      dprimary_args = list(r = 0.5),
      D = D, shape = shape, scale = scale
    )
    theoretical_mean <- sum(x_values * pmf)
    theoretical_sd <- sqrt(sum((x_values - theoretical_mean)^2 * pmf))

    # Compare empirical and theoretical statistics
    expect_equal(empirical_mean, theoretical_mean, tolerance = 0.1)
    expect_equal(empirical_sd, theoretical_sd, tolerance = 0.1)

    # Check empirical PMF against theoretical PMF
    empirical_pmf <- as.vector(table(samples) / n)
    expect_equal(empirical_pmf, pmf, tolerance = 0.01)

    # Check empirical CDF against theoretical CDF
    empirical_cdf <- ecdf(samples)(x_values)
    theoretical_cdf <- ppcens(
      c(x_values[-1], D), pweibull,
      pwindow = pwindow, D = D,
      dprimary = dexpgrowth,
      dprimary_args = list(r = 0.5),
      shape = shape, scale = scale
    )
    expect_equal(cumsum(pmf), theoretical_cdf, tolerance = 0.01)
    expect_equal(empirical_cdf, theoretical_cdf, tolerance = 0.01)
  }
)

test_that(
  "Exponential distribution with daily censoring matches analytical solution",
  {
    # Parameters
    rate <- 1 # Exponential rate parameter
    pwindow <- 1 # Primary event window (1 day)
    swindow <- 1 # Secondary event window (1 day)
    D <- 10 # Maximum delay to consider

    # Analytical solution
    analytical_pmf <- function(t) {
      if (t == 0) {
        exp(-1)
      } else {
        (1 - exp(-1)) * (exp(1) - 1) * exp(-t)
      }
    }
    # Numerical solution using our functions
    numerical_pmf <- dpcens(
      0:(D - 1), pexp,
      pwindow = pwindow, swindow = swindow, D = 20,
      rate = rate
    )
    numerical_cdf <- ppcens(
      1:D, pexp,
      pwindow = pwindow, D = 20,
      rate = rate
    )

    # Compare PMF
    analytical_values <- sapply(0:D, analytical_pmf)
    analytical_values <- analytical_values[-length(analytical_values)] /
      sum(analytical_values)
    expect_equal(numerical_pmf, analytical_values, tolerance = 1e-3)

    # Compare CDF
    analytical_cdf <- cumsum(analytical_values)
    expect_equal(numerical_cdf, analytical_cdf, tolerance = 1e-3)

    # Check that PMF sums to 1 (approximately)
    expect_equal(sum(numerical_pmf), 1, tolerance = 1e-3)

    # Check that CDF approaches 1
    expect_equal(numerical_cdf[length(numerical_cdf)], 1, tolerance = 1e-3)

    # Compare random number generator to analytical solution
    samples <- rpcens(
      10000, rexp,
      pwindow = pwindow, swindow = swindow, D = D,
      rate = rate
    )
    empirical_pmf <- as.vector(table(samples) / 10000)
    expect_equal(
      empirical_pmf, analytical_values[seq_along(empirical_pmf)],
      tolerance = 0.05
    )
  }
)

# Left truncation (L parameter) tests

test_that(
  "rprimarycensored is consistent with dprimarycensored and
   pprimarycensored with L > 0",
  { # nolint
    n <- 1e5
    pwindow <- 1
    D <- 10
    L <- 2
    samples <- rpcens(
      n, rlnorm, pwindow,
      D = D, L = L, meanlog = 1, sdlog = 0.5
    )

    # Check samples are in valid range
    expect_true(all(samples >= L & samples < D))

    # Check empirical PMF against theoretical PMF
    x_values <- L:(D - 1)
    empirical_pmf <- as.vector(table(factor(samples, levels = x_values)) / n)
    theoretical_pmf <- dpcens(
      x_values, plnorm, pwindow,
      D = D, L = L, meanlog = 1, sdlog = 0.5
    )
    expect_equal(empirical_pmf, theoretical_pmf, tolerance = 0.02)

    # Check empirical CDF against theoretical CDF
    empirical_cdf <- ecdf(samples)(x_values)
    theoretical_cdf <- ppcens(
      x_values + 1, plnorm, pwindow,
      D = D, L = L, meanlog = 1, sdlog = 0.5
    )
    # Adjust last value since CDF at D should be 1
    theoretical_cdf[length(theoretical_cdf)] <- 1
    expect_equal(empirical_cdf, theoretical_cdf, tolerance = 0.02)

    # Check empirical mean matches theoretical mean
    empirical_mean <- mean(samples)
    theoretical_mean <- sum(x_values * theoretical_pmf)
    expect_equal(empirical_mean, theoretical_mean, tolerance = 0.05)
  }
)

test_that(
  "rprimarycensored is consistent with dprimarycensored and
   pprimarycensored with L > 0 and exponential growth primary",
  { # nolint
    n <- 1e5
    pwindow <- 2
    D <- 12
    L <- 3
    r <- 0.3
    samples <- rpcens(
      n, rlnorm, pwindow,
      D = D, L = L,
      rprimary = rexpgrowth,
      rprimary_args = list(r = r),
      meanlog = 1.5, sdlog = 0.5
    )

    # Check samples are in valid range
    expect_true(all(samples >= L & samples < D))

    # Check empirical PMF against theoretical PMF
    x_values <- L:(D - 1)
    empirical_pmf <- as.vector(table(factor(samples, levels = x_values)) / n)
    theoretical_pmf <- dpcens(
      x_values, plnorm, pwindow,
      D = D, L = L,
      dprimary = dexpgrowth,
      dprimary_args = list(r = r),
      meanlog = 1.5, sdlog = 0.5
    )
    expect_equal(empirical_pmf, theoretical_pmf, tolerance = 0.02)

    # Check empirical mean matches theoretical mean
    empirical_mean <- mean(samples)
    theoretical_mean <- sum(x_values * theoretical_pmf)
    expect_equal(empirical_mean, theoretical_mean, tolerance = 0.1)
  }
)
