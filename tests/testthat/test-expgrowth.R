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

  test_that(
    paste("rexpgrowth generates samples in correct range:", tc$label),
    {
      samples <- rexpgrowth(1000, tc$min, tc$max, tc$r)
      expect_true(all(samples >= tc$min & samples <= tc$max))
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

  test_that(
    paste("pexpgrowth boundary values:", tc$label),
    {
      expect_equal(
        pexpgrowth(tc$max, tc$min, tc$max, tc$r), 1,
        tolerance = 1e-10
      )
      expect_equal(
        pexpgrowth(tc$min, tc$min, tc$max, tc$r), 0,
        tolerance = 1e-10
      )
    }
  )

  test_that(
    paste("pexpgrowth lower.tail sums to 1:", tc$label),
    {
      x_values <- seq(tc$min, tc$max, length.out = 50)
      cdf_lower <- pexpgrowth(
        x_values, tc$min, tc$max, tc$r,
        lower.tail = TRUE
      )
      cdf_upper <- pexpgrowth(
        x_values, tc$min, tc$max, tc$r,
        lower.tail = FALSE
      )
      expect_equal(
        cdf_lower + cdf_upper, rep(1, length(x_values)),
        tolerance = 1e-10
      )
    }
  )

  test_that(
    paste("rexpgrowth mean approximates theoretical mean:", tc$label),
    {
      n <- 100000
      samples <- rexpgrowth(n, tc$min, tc$max, tc$r)

      theoretical_mean <- (
        tc$r * tc$max * exp(tc$r * tc$max) -
          tc$r * tc$min * exp(tc$r * tc$min) -
          exp(tc$r * tc$max) + exp(tc$r * tc$min)
      ) / (tc$r * (exp(tc$r * tc$max) - exp(tc$r * tc$min)))

      expect_equal(mean(samples), theoretical_mean, tolerance = 0.01)
    }
  )
}

test_that("expgrowth functions handle very small r correctly", {
  min <- 0
  max <- 1
  r <- 1e-11
  n <- 100000

  samples <- rexpgrowth(n, min, max, r)
  expect_true(all(samples >= min & samples <= max))

  # For very small r, the distribution should be close to uniform
  expect_equal(mean(samples), 0.5, tolerance = 0.01)
  expect_equal(var(samples), 1 / 12, tolerance = 0.02)

  x_values <- seq(min, max, length.out = 100)
  densities <- dexpgrowth(x_values, min, max, r)
  expect_true(all(densities >= 0.99 & densities <= 1.01))

  cdfs <- pexpgrowth(x_values, min, max, r)
  expect_equal(cdfs, x_values, tolerance = 0.01)

  # Consistency check
  empirical_cdf <- ecdf(samples)
  theoretical_cdf <- pexpgrowth(x_values, min, max, r)
  expect_equal(empirical_cdf(x_values), theoretical_cdf, tolerance = 0.01)
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
