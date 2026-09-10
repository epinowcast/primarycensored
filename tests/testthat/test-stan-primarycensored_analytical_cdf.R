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
              q, 2, c(shape, rate), pwindow, 0, Inf, 1, numeric(0)
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
              q, 1, c(meanlog, sdlog), pwindow, 0, Inf, 1, numeric(0)
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
              q, 3, c(shape, scale), pwindow, 0, Inf, 1, numeric(0)
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

test_that(
  "Stan primarycensored_analytical_lcdf matches R implementation for
   generalised gamma",
  {
    skip_if_not_installed("flexsurv")
    shapes <- c(0.5, 1, 2)
    scales <- c(0.5, 2, 5)
    ks <- c(0.5, 1, 3)
    pwindows <- c(1, 2, 5)

    for (shape in shapes) {
      for (scale in scales) {
        for (k in ks) {
          for (pwindow in pwindows) {
            obj <- new_pcens(
              flexsurv::pgengamma.orig,
              dunif, list(),
              shape = shape, scale = scale, k = k
            )

            q_values <- seq(0, 30, by = 1)
            r_result <- pcens_cdf(
              obj,
              q = q_values, pwindow = pwindow, use_numeric = FALSE
            )

            stan_result <- vapply(q_values, function(q) {
              primarycensored_analytical_cdf(
                q, 5, c(shape, scale, k), pwindow, 0, Inf, 1, numeric(0)
              )
            }, numeric(1))

            stan_result <- ifelse(is.nan(stan_result), 1, stan_result)
            expect_equal(
              r_result, stan_result,
              tolerance = 1e-6,
              info = sprintf(
                "Mismatch for shape = %s, scale = %s, k = %s, pwindow = %s",
                shape, scale, k, pwindow
              )
            )
          }
        }
      }
    }
  }
)

test_that(
  "Stan analytical generalised gamma matches the numeric CDF and the gamma and
   Weibull analytical solutions",
  {
    skip_if_not_installed("flexsurv")
    q_values <- seq(0, 20, by = 1)
    pwindow <- 2
    expect_identical(check_for_analytical(5L, 1L), 1L)

    analytical <- vapply(q_values, function(q) {
      primarycensored_analytical_cdf(
        q, 5, c(1.5, 2, 0.8), pwindow, 0, Inf, 1, numeric(0)
      )
    }, numeric(1))
    numeric <- vapply(q_values, function(q) {
      primarycensored_cdf(
        q, 5, c(1.5, 2, 0.8), pwindow, 0, Inf, 1, numeric(0)
      )
    }, numeric(1))
    expect_equal(analytical, numeric, tolerance = 1e-6)

    # shape = 1 is the gamma distribution with shape k and rate 1 / scale
    gamma_result <- vapply(q_values, function(q) {
      primarycensored_analytical_cdf(
        q, 2, c(0.8, 1 / 2), pwindow, 0, Inf, 1, numeric(0)
      )
    }, numeric(1))
    gengamma_result <- vapply(q_values, function(q) {
      primarycensored_analytical_cdf(
        q, 5, c(1, 2, 0.8), pwindow, 0, Inf, 1, numeric(0)
      )
    }, numeric(1))
    expect_equal(gengamma_result, gamma_result, tolerance = 1e-8)

    # k = 1 is the Weibull distribution
    weibull_result <- vapply(q_values, function(q) {
      primarycensored_analytical_cdf(
        q, 3, c(1.5, 2), pwindow, 0, Inf, 1, numeric(0)
      )
    }, numeric(1))
    gengamma_result <- vapply(q_values, function(q) {
      primarycensored_analytical_cdf(
        q, 5, c(1.5, 2, 1), pwindow, 0, Inf, 1, numeric(0)
      )
    }, numeric(1))
    expect_equal(gengamma_result, weibull_result, tolerance = 1e-8)
  }
)
