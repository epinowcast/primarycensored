skip_on_cran()

test_that(
  "Stan primarycensored_analytical_lcdf matches R implementation for
   Gamma",
  {
    shapes <- c(0.5, 1, 2, 5)
    rates <- c(0.1, 0.5, 1, 2)
    pwindows <- c(1, 2, 5, 10)

    for (shape in shapes) {
      for (rate in rates) {
        for (pwindow in pwindows) {
          obj <- new_pcens(
            pgamma,
            dunif, list(),
            shape = shape, rate = rate
          )

          q_values <- seq(0, 30, by = 1)
          r_result <- pcens_cdf(
            obj,
            q = q_values, pwindow = pwindow, use_numeric = FALSE
          )

          stan_result <- vapply(q_values, function(q) {
            primarycensored_analytical_cdf(
              q, 2, c(shape, rate), pwindow, Inf, 0, 1, numeric(0)
            )
          }, numeric(1))

          stan_result <- ifelse(is.nan(stan_result), 1, stan_result)
          expect_equal(
            r_result, stan_result,
            tolerance = 1e-6,
            info = sprintf(
              "Mismatch for shape = %s, rate = %s, pwindow = %s",
              shape, rate, pwindow
            )
          )
        }
      }
    }
  }
)

test_that(
  "Stan primarycensored_analytical_lcdf matches R implementation for
   Lognormal",
  {
    meanlogs <- c(-1, 0, 1, 2)
    sdlogs <- c(0.5, 1, 1.5)
    pwindows <- c(1, 2, 5, 8)

    for (meanlog in meanlogs) {
      for (sdlog in sdlogs) {
        for (pwindow in pwindows) {
          obj <- new_pcens(
            plnorm,
            dunif, list(),
            meanlog = meanlog, sdlog = sdlog
          )

          q_values <- seq(0, 30, by = 1)
          r_result <- pcens_cdf(
            obj,
            q = q_values, pwindow = pwindow, use_numeric = FALSE
          )

          stan_result <- vapply(q_values, function(q) {
            primarycensored_analytical_cdf(
              q, 1, c(meanlog, sdlog), pwindow, Inf, 0, 1, numeric(0)
            )
          }, numeric(1))

          stan_result <- ifelse(is.nan(stan_result), 1, stan_result)
          expect_equal(
            r_result, stan_result,
            tolerance = 1e-6,
            info = sprintf(
              "Mismatch for meanlog = %s, sdlog = %s, pwindow = %s",
              meanlog, sdlog, pwindow
            )
          )
        }
      }
    }
  }
)

test_that(
  "Stan primarycensored_analytical_lcdf matches R implementation for
   Weibull",
  {
    shapes <- c(0.5, 1, 2, 3)
    scales <- c(0.5, 1, 2, 5)
    pwindows <- c(1, 2, 5, 10)

    for (shape in shapes) {
      for (scale in scales) {
        for (pwindow in pwindows) {
          obj <- new_pcens(
            pweibull,
            dunif, list(),
            shape = shape, scale = scale
          )

          q_values <- seq(0, 30, by = 1)
          r_result <- pcens_cdf(
            obj,
            q = q_values, pwindow = pwindow, use_numeric = FALSE
          )

          stan_result <- vapply(q_values, function(q) {
            primarycensored_analytical_cdf(
              q, 3, c(shape, scale), pwindow, Inf, 0, 1, numeric(0)
            )
          }, numeric(1))

          stan_result <- ifelse(is.nan(stan_result), 1, stan_result)
          expect_equal(
            r_result, stan_result,
            tolerance = 1e-6,
            info = sprintf(
              "Mismatch for shape = %s, scale = %s, pwindow = %s",
              shape, scale, pwindow
            )
          )
        }
      }
    }
  }
)
