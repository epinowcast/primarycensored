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
      0:(D - 1), pexp, pwindow, swindow, 20,
      rate = rate
    )
    numerical_cdf <- ppcens(
      1:D, pexp, pwindow, 20,
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
      10000, rexp, pwindow, swindow, D, rate = rate
    )
    empirical_pmf <- as.vector(table(samples) / 10000)
    expect_equal(
      empirical_pmf, analytical_values[1:length(empirical_pmf)],
      tolerance = 0.05
    )
  }
)
