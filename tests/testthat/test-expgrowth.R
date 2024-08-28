test_that("dexpgrowth integrates to 1", {
  min <- 0
  max <- 1
  r <- 0.5
  integral <- integrate(function(x) dexpgrowth(x, min, max, r), min, max)$value
  expect_equal(integral, 1, tolerance = 1e-6)
})

test_that("pexpgrowth matches integral of dexpgrowth", {
  min <- 0
  max <- 2
  r <- 1.5
  q <- c(0.5, 1, 1.5)

  for (x in q) {
    cdf_integral <- integrate(function(t) {
      dexpgrowth(t, min, max, r)
    }, min, x)$value
    cdf_pexpgrowth <- pexpgrowth(x, min, max, r)
    expect_equal(cdf_integral, cdf_pexpgrowth, tolerance = 1e-6)
  }
})

test_that("rexpgrowth generates samples within the correct range", {
  min <- 1
  max <- 5
  r <- 0.8
  n <- 1000
  samples <- rexpgrowth(n, min, max, r)

  expect_true(all(samples >= min & samples < max))
})

test_that("rexpgrowth mean approximates theoretical mean", {
  min <- 0
  max <- 1
  r <- 2
  n <- 100000
  samples <- rexpgrowth(n, min, max, r)

  theoretical_mean <- (
    r * max * exp(r * max) - r * min * exp(r * min) -
      exp(r * max) + exp(r * min)
  ) / (r * (exp(r * max) - exp(r * min)))
  sample_mean <- mean(samples)

  expect_equal(sample_mean, theoretical_mean, tolerance = 0.01)
})

test_that("dexpgrowth, pexpgrowth, and rexpgrowth are consistent", {
  min <- 0
  max <- 2
  r <- 1.2
  n <- 10000
  samples <- rexpgrowth(n, min, max, r)

  # Compare empirical CDF with theoretical CDF
  empirical_cdf <- ecdf(samples)
  x_values <- seq(min, max, length.out = 100)
  theoretical_cdf <- pexpgrowth(x_values, min, max, r)

  expect_equal(empirical_cdf(x_values), theoretical_cdf, tolerance = 0.05)

  # Compare empirical PDF with theoretical PDF using histogram
  hist_data <- hist(samples, breaks = 50, plot = FALSE)
  midpoints <- (
    hist_data$breaks[-1] + hist_data$breaks[-length(hist_data$breaks)]
  ) / 2
  empirical_pdf <- hist_data$density
  theoretical_pdf <- dexpgrowth(midpoints, min, max, r)

  expect_equal(empirical_pdf, theoretical_pdf, tolerance = 0.1)
})
