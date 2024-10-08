test_that(
  "Stan primary_censored_dist_analytical_lcdf matches R implementation for
   Gamma",
  {
    shapes <- c(0.5, 1, 2, 5)
    rates <- c(0.1, 0.5, 1, 2)
    pwindows <- c(1, 2, 5, 10)

    for (shape in shapes) {
      for (rate in rates) {
        for (pwindow in pwindows) {
          obj <- new_primary_censored_dist(
            pgamma,
            dunif, list(),
            "pgamma", "dunif",
            shape = shape, rate = rate
          )

          q_values <- seq(0, 30, by = 1)
          r_result <- primary_censored_cdf(
            obj,
            q = q_values, pwindow = pwindow, use_numeric = FALSE
          )

          stan_result <- vapply(q_values, function(q) {
            primary_censored_dist_analytical_cdf(
              q, 2, c(shape, rate), pwindow, Inf, 1, numeric(0)
            )
          }, numeric(1))

          stan_result <- ifelse(is.nan(stan_result), 1, stan_result)
          expect_equal(r_result, stan_result, tolerance = 1e-6)
        }
      }
    }
  }
)

test_that(
  "Stan primary_censored_dist_analytical_lcdf matches R implementation for
   Lognormal",
  {
    meanlogs <- c(-1, 0, 1, 2)
    sdlogs <- c(0.5, 1, 1.5)
    pwindows <- c(1, 2, 5, 8)

    for (meanlog in meanlogs) {
      for (sdlog in sdlogs) {
        for (pwindow in pwindows) {
          obj <- new_primary_censored_dist(
            plnorm,
            dunif, list(),
            "plnorm", "dunif",
            meanlog = meanlog, sdlog = sdlog
          )

          q_values <- seq(0, 30, by = 1)
          r_result <- primary_censored_cdf(
            obj,
            q = q_values, pwindow = pwindow, use_numeric = FALSE
          )

          stan_result <- vapply(q_values, function(q) {
            primary_censored_dist_analytical_cdf(
              q, 1, c(meanlog, sdlog), pwindow, Inf, 1, numeric(0)
            )
          }, numeric(1))

          stan_result <- ifelse(is.nan(stan_result), 1, stan_result)
          expect_equal(r_result, stan_result, tolerance = 1e-6)
        }
      }
    }
  }
)

test_that(
  "Stan primary_censored_dist_analytical_lcdf matches R implementation for
   Weibull",
  {
    shapes <- c(0.5, 1, 2, 3)
    scales <- c(0.5, 1, 2, 5)
    pwindows <- c(1, 2, 5, 10)

    for (shape in shapes) {
      for (scale in scales) {
        for (pwindow in pwindows) {
          obj <- new_primary_censored_dist(
            pweibull,
            dunif, list(),
            "pweibull", "dunif",
            shape = shape, scale = scale
          )

          q_values <- seq(0, 30, by = 1)
          r_result <- primary_censored_cdf(
            obj,
            q = q_values, pwindow = pwindow, use_numeric = FALSE
          )

          stan_result <- vapply(q_values, function(q) {
            primary_censored_dist_analytical_cdf(
              q, 3, c(shape, scale), pwindow, Inf, 1, numeric(0)
            )
          }, numeric(1))

          stan_result <- ifelse(is.nan(stan_result), 1, stan_result)
          expect_equal(r_result, stan_result, tolerance = 1e-6)
        }
      }
    }
  }
)
