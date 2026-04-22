skip_on_cran()

# Parameter sets covering positive r, negative r, and non-zero min
stan_test_cases <- list(
  list(
    x = seq(0, 1, by = 0.1), min = 0, max = 1,
    r = 0.5, label = "positive r"
  ),
  list(
    x = seq(0, 1, by = 0.1), min = 0, max = 1,
    r = -0.5, label = "negative r"
  ),
  list(
    x = seq(10, 20, by = 1), min = 10, max = 20,
    r = 0.1, label = "non-zero min"
  )
)

for (tc in stan_test_cases) {
  test_that(
    paste("Stan expgrowth_pdf matches R dexpgrowth:", tc$label),
    {
      stan_pdf <- sapply(tc$x, expgrowth_pdf, tc$min, tc$max, tc$r)
      r_pdf <- dexpgrowth(tc$x, tc$min, tc$max, tc$r)
      expect_equal(stan_pdf, r_pdf, tolerance = 1e-4)
    }
  )

  test_that(
    paste("Stan expgrowth_lpdf matches R dexpgrowth:", tc$label),
    {
      stan_lpdf <- sapply(
        tc$x, expgrowth_lpdf, tc$min, tc$max, tc$r
      )
      r_lpdf <- dexpgrowth(tc$x, tc$min, tc$max, tc$r, log = TRUE)
      expect_equal(stan_lpdf, r_lpdf, tolerance = 1e-4)
    }
  )

  test_that(
    paste("Stan expgrowth_cdf matches R pexpgrowth:", tc$label),
    {
      stan_cdf <- sapply(tc$x, expgrowth_cdf, tc$min, tc$max, tc$r)
      r_cdf <- pexpgrowth(tc$x, tc$min, tc$max, tc$r)
      expect_equal(stan_cdf, r_cdf, tolerance = 1e-4)
    }
  )

  test_that(
    paste("Stan expgrowth_lcdf matches R pexpgrowth:", tc$label),
    {
      # Avoid boundary points where log(0) = -Inf
      x_inner <- tc$x[tc$x > tc$min & tc$x < tc$max]
      stan_lcdf <- sapply(
        x_inner, expgrowth_lcdf, tc$min, tc$max, tc$r
      )
      r_lcdf <- pexpgrowth(
        x_inner, tc$min, tc$max, tc$r,
        log.p = TRUE
      )
      expect_equal(stan_lcdf, r_lcdf, tolerance = 1e-6)
    }
  )

  test_that(
    paste("Stan expgrowth_rng matches theoretical CDF:", tc$label),
    {
      n <- 10000
      stan_samples <- replicate(
        n, expgrowth_rng(tc$min, tc$max, tc$r)
      )

      expect_true(all(
        stan_samples >= tc$min & stan_samples <= tc$max
      ))

      empirical_cdf <- ecdf(stan_samples)
      x_values <- seq(tc$min, tc$max, length.out = 50)
      theoretical_cdf <- pexpgrowth(
        x_values, tc$min, tc$max, tc$r
      )
      expect_equal(
        empirical_cdf(x_values), theoretical_cdf,
        tolerance = 0.05
      )
    }
  )
}
