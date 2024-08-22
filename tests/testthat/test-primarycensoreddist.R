test_that("pprimarycensoreddist integrates to 1", {
  pwindow <- 5
  D <- 10
  integral <- integrate(
    function(x) {
      dprimarycensoreddist(x, plnorm, pwindow, D = D, meanlog = 0, sdlog = 1)
    },
    0, D
  )$value
  expect_equal(integral, 1, tolerance = 1e-6)
})

test_that("dprimarycensoreddist matches difference of pprimarycensoreddist", {
  x <- c(1, 2, 3)
  pwindow <- 5
  swindow <- 0.5
  D <- 10

  pmf <- dprimarycensoreddist(
    x, plnorm, pwindow, swindow, D,
    meanlog = 0, sdlog = 1
  )
  cdf_diff <- sapply(x, function(xi) {
    pprimarycensoreddist(
      xi + swindow, plnorm, pwindow, D,
      meanlog = 0, sdlog = 1
    ) -
      pprimarycensoreddist(
        xi, plnorm, pwindow, D,
        meanlog = 0, sdlog = 1
      )
  })

  expect_equal(pmf, cdf_diff, tolerance = 1e-6)
})

test_that("rprimarycensoreddist generates samples within the correct range", {
  n <- 1000
  pwindow <- 5
  D <- 10
  samples <- rprimarycensoreddist(
    n, rlnorm, pwindow,
    D = D, meanlog = 0, sdlog = 1
  )

  expect_true(all(samples > 0 & samples <= D))
})

test_that("rprimarycensoreddist mean approximates theoretical mean", {
  n <- 100000
  pwindow <- 5
  D <- 10
  samples <- rprimarycensoreddist(
    n, rlnorm, pwindow,
    D = D, meanlog = 0, sdlog = 1
  )

  theoretical_mean <- integrate(
    function(x) {
      x * dprimarycensoreddist(x, plnorm, pwindow, D = D, meanlog = 0, sdlog = 1)
    },
    0, D
  )$value
  sample_mean <- mean(samples)

  expect_equal(sample_mean, theoretical_mean, tolerance = 0.05)
})

test_that(
  "pprimarycensoreddist, dprimarycensoreddist, and rprimarycensoreddist are
   consistent",
  {
    n <- 10000
    pwindow <- 5
    D <- 10
    samples <- rprimarycensoreddist(
      n, rlnorm, pwindow,
      D = D, meanlog = 0, sdlog = 1
    )

    # Compare empirical CDF with theoretical CDF
    empirical_cdf <- ecdf(samples)
    x_values <- seq(0, D, length.out = 100)
    theoretical_cdf <- pprimarycensoreddist(
      x_values, plnorm, pwindow, D,
      meanlog = 0, sdlog = 1
    )

    expect_equal(empirical_cdf(x_values), theoretical_cdf, tolerance = 0.05)

    # Compare empirical PDF with theoretical PDF using histogram
    hist_data <- hist(samples, breaks = 50, plot = FALSE)
    midpoints <- (
      hist_data$breaks[-1] + hist_data$breaks[-length(hist_data$breaks)]
    ) / 2
    empirical_pdf <- hist_data$density
    theoretical_pdf <- dprimarycensoreddist(
      midpoints, plnorm, pwindow,
      D = D, meanlog = 0, sdlog = 1
    )

    expect_equal(empirical_pdf, theoretical_pdf, tolerance = 0.1)
  }
)

test_that(
  "primarycensoreddist functions work with exponential growth primary
   distribution",
  {
    n <- 10000
    pwindow <- 5
    D <- 10
    r <- 0.5

    samples <- rprimarycensoreddist(
      n, rlnorm, pwindow,
      D = D,
      rprimary = rexpgrowth,
      rprimary_args = list(min = 0, max = pwindow, r = r),
      meanlog = 0, sdlog = 1
    )

    # Check CDF
    empirical_cdf <- ecdf(samples)
    x_values <- seq(0, D, length.out = 100)
    theoretical_cdf <- pprimarycensoreddist(
      x_values, plnorm, pwindow, D,
      dprimary = dexpgrowth,
      dprimary_args = list(min = 0, max = pwindow, r = r),
      meanlog = 0, sdlog = 1
    )

    expect_equal(empirical_cdf(x_values), theoretical_cdf, tolerance = 0.05)

    # Check PDF
    hist_data <- hist(samples, breaks = 50, plot = FALSE)
    midpoints <- (
      hist_data$breaks[-1] + hist_data$breaks[-length(hist_data$breaks)]
    ) / 2
    empirical_pdf <- hist_data$density
    theoretical_pdf <- dprimarycensoreddist(
      midpoints, plnorm, pwindow,
      D = D,
      dprimary = dexpgrowth,
      dprimary_args = list(min = 0, max = pwindow, r = r),
      meanlog = 0, sdlog = 1
    )

    expect_equal(empirical_pdf, theoretical_pdf, tolerance = 0.1)
  }
)

test_that(
  "Exponential distribution with daily censoring matches analytical solution",
  {
    # Parameters
    rate <- 1 # Exponential rate parameter
    pwindow <- 1 # Primary event window (1 day)
    swindow <- 1 # Secondary event window (1 day)
    D <- 20 # Maximum delay to consider

    # Analytical solution
    analytical_pmf <- function(s) {
      if (s == 0) {
        exp(-1)
      } else {
        (1 - exp(-1)) * (exp(1) - 1) * exp(-s)
      }
    }

    # Numerical solution using our functions
    numerical_pmf <- dprimarycensoreddist(
      0:D, dexp, pwindow, swindow, D,
      rate = rate
    )
    numerical_cdf <- pprimarycensoreddist(
      0:D, dexp, pwindow, D,
      rate = rate
    )

    # Compare PMF
    analytical_values <- sapply(0:D, analytical_pmf)
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
    samples <- rprimarycensoreddist(n, rexp, pwindow, swindow, D, rate = rate)

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
    analytical_values <- sapply(as.numeric(names(empirical_pmf)), analytical_pmf)

    expect_equal(as.numeric(empirical_pmf), analytical_values, tolerance = 0.01)
  }
)
