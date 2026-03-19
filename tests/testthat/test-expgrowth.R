# Parameter sets covering positive r, negative r, and non-zero min
test_cases <- list(
  list(min = 0, max = 1, r = 0.5, label = "positive r"),
  list(min = 0, max = 2, r = 1.5, label = "positive r wide"),
  list(min = 0, max = 1, r = -0.5, label = "negative r"),
  list(min = 0, max = 2, r = -1.5, label = "negative r wide"),
  list(min = 10, max = 20, r = 0.1, label = "non-zero min positive r"),
  list(min = 5, max = 15, r = -0.3, label = "non-zero min negative r")
)

for (tc in test_cases) {
  test_that(
    paste("dexpgrowth integrates to 1:", tc$label),
    {
      integral <- integrate(
        function(x) dexpgrowth(x, tc$min, tc$max, tc$r),
        tc$min, tc$max
      )$value
      expect_equal(integral, 1, tolerance = 1e-6)
    }
  )

  test_that(
    paste("pexpgrowth matches integral of dexpgrowth:", tc$label),
    {
      q <- seq(tc$min, tc$max, length.out = 5)[2:4]
      for (x in q) {
        cdf_integral <- integrate(
          function(t) dexpgrowth(t, tc$min, tc$max, tc$r),
          tc$min, x
        )$value
        expect_equal(
          pexpgrowth(x, tc$min, tc$max, tc$r),
          cdf_integral,
          tolerance = 1e-6
        )
      }
    }
  )
}

# Parameterised range and consistency tests
consistency_cases <- list(
  list(min = 0, max = 2, r = 1.2, label = "positive r"),
  list(min = 0, max = 2, r = -1.2, label = "negative r"),
  list(min = 10, max = 20, r = 0.2, label = "non-zero min")
)

for (tc in consistency_cases) {
  test_that(
    paste("rexpgrowth generates samples in correct range:", tc$label),
    {
      samples <- rexpgrowth(1000, tc$min, tc$max, tc$r)
      expect_true(all(samples >= tc$min & samples < tc$max))
    }
  )

  test_that(
    paste("d, p, r expgrowth are consistent:", tc$label),
    {
      set.seed(42)
      n <- 10000
      samples <- rexpgrowth(n, tc$min, tc$max, tc$r)

      empirical_cdf <- ecdf(samples)
      x_values <- seq(tc$min, tc$max, length.out = 50)
      theoretical_cdf <- pexpgrowth(x_values, tc$min, tc$max, tc$r)
      expect_equal(
        empirical_cdf(x_values), theoretical_cdf,
        tolerance = 0.05
      )
    }
  )
}

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

  expect_equal(mean(samples), theoretical_mean, tolerance = 0.01)
})

test_that("expgrowth functions handle very small r correctly", {
  min <- 0
  max <- 1
  r <- 1e-11
  n <- 100000

  samples <- rexpgrowth(n, min, max, r)
  expect_true(all(samples >= min & samples < max))

  # For very small r, the distribution should be close to uniform
  expect_equal(mean(samples), 0.5, tolerance = 0.01)
  expect_equal(var(samples), 1 / 12, tolerance = 0.02)

  x_values <- seq(min, max, length.out = 100)
  densities <- dexpgrowth(x_values, min, max, r)
  expect_true(all(densities >= 0.99 & densities <= 1.01))

  cdfs <- pexpgrowth(x_values, min, max, r)
  expect_equal(cdfs, x_values, tolerance = 0.01)
})

test_that("dexpgrowth with log = TRUE matches log of density", {
  min <- 0
  max <- 1
  r <- -0.5
  x <- seq(0, 1, by = 0.1)

  log_density <- dexpgrowth(x, min, max, r, log = TRUE)
  density <- dexpgrowth(x, min, max, r, log = FALSE)

  expect_equal(log_density, log(density), tolerance = 1e-10)
  expect_true(all(is.finite(log_density)))
})

test_that("pexpgrowth boundary values are correct with non-zero min", {
  min <- 10
  max <- 20
  r <- 0.1
  expect_equal(pexpgrowth(max, min, max, r), 1, tolerance = 1e-10)
  expect_equal(pexpgrowth(min, min, max, r), 0, tolerance = 1e-10)
})

test_that("pexpgrowth handles lower.tail argument correctly", {
  min <- 0
  max <- 10
  r <- 0.5
  x_values <- seq(min, max, length.out = 100)

  cdf_lower <- pexpgrowth(x_values, min, max, r, lower.tail = TRUE)
  cdf_upper <- pexpgrowth(x_values, min, max, r, lower.tail = FALSE)

  expect_equal(
    cdf_lower + cdf_upper, rep(1, length(x_values)),
    tolerance = 1e-10
  )

  expect_identical(pexpgrowth(min, min, max, r, lower.tail = TRUE), 0)
  expect_identical(pexpgrowth(min, min, max, r, lower.tail = FALSE), 1)
  expect_identical(pexpgrowth(max, min, max, r, lower.tail = TRUE), 1)
  expect_identical(pexpgrowth(max, min, max, r, lower.tail = FALSE), 0)
})
