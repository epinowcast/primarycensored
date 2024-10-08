test_that("new_primary_censored_dist creates object with correct structure", {
  pdist_name <- "pgamma"
  pdist <- pgamma
  dprimary_name <- "dunif"
  dprimary <- dunif
  shape <- 2
  rate <- 1

  obj <- new_primary_censored_dist(
    pdist,
    dprimary, list(),
    pdist_name, dprimary_name,
    shape = shape, rate = rate
  )

  expect_s3_class(obj, "pcens_pgamma_dunif")
  expect_identical(obj$pdist, pgamma)
  expect_identical(obj$dprimary, dunif)
  expect_identical(obj$args, list(shape = shape, rate = rate))

  new_obj <- new_primary_censored_dist(
    pgamma, dunif, list(),
    shape = shape, rate = rate
  )
  expect_identical(obj, new_obj)
})

test_that(
  "primary_censored_cdf methods dispatch correctly to existing
   analytical solutions",
  {
    pdist_name <- "pgamma"
    pdist <- pgamma
    dprimary_name <- "dunif"
    dprimary <- dunif

    obj_gamma <- new_primary_censored_dist(
      pdist, dprimary, list(),
      pdist_name, dprimary_name,
      shape = 2, rate = 1
    )

    pdist_name <- "plnorm"
    pdist <- plnorm
    dprimary_name <- "dunif"
    dprimary <- dunif

    obj_lnorm <- new_primary_censored_dist(
      pdist, dprimary, list(),
      pdist_name, dprimary_name,
      meanlog = 0, sdlog = 1
    )

    pdist_name <- "pweibull"
    pdist <- pweibull
    dprimary_name <- "dunif"
    dprimary <- dunif

    obj_weibull <- new_primary_censored_dist(
      pdist, dprimary, list(),
      pdist_name, dprimary_name,
      shape = 2, scale = 1
    )

    expect_s3_class(obj_gamma, "pcens_pgamma_dunif")
    expect_s3_class(obj_lnorm, "pcens_plnorm_dunif")
    expect_s3_class(obj_weibull, "pcens_pweibull_dunif")

    q_values <- c(5, 10)
    pwindow <- 2

    expect_no_error(
      primary_censored_cdf(obj_gamma, q = q_values, pwindow = pwindow)
    )
    expect_no_error(
      primary_censored_cdf(obj_lnorm, q = q_values, pwindow = pwindow)
    )
    expect_no_error(
      primary_censored_cdf(obj_weibull, q = q_values, pwindow = pwindow)
    )
  }
)

test_that(
  "primary_censored_cdf errors as expected when the wrong distributional
   parameters are supplied",
  {
    pdist_name <- "pgamma"
    pdist <- pgamma
    dprimary_name <- "dunif"
    dprimary <- dunif

    obj_gamma <- new_primary_censored_dist(
      pdist, dprimary, list(),
      pdist_name, dprimary_name,
      rate = 1
    )

    expect_error(
      primary_censored_cdf(obj_gamma, q = 1, pwindow = 1),
      "shape parameter is required for Gamma distribution"
    )

    obj_gamma_no_rate <- new_primary_censored_dist(
      pdist, dprimary, list(),
      pdist_name, dprimary_name,
      shape = 2
    )

    expect_error(
      primary_censored_cdf(obj_gamma_no_rate, q = 1, pwindow = 1),
      "scale or rate parameter is required for Gamma distribution"
    )

    pdist_name <- "plnorm"
    pdist <- plnorm

    obj_lnorm_no_meanlog <- new_primary_censored_dist(
      pdist, dprimary, list(),
      pdist_name, dprimary_name,
      sdlog = 1
    )

    expect_error(
      primary_censored_cdf(obj_lnorm_no_meanlog, q = 1, pwindow = 1),
      "meanlog parameter is required for Log-Normal distribution"
    )

    obj_lnorm_no_sdlog <- new_primary_censored_dist(
      pdist, dprimary, list(),
      pdist_name, dprimary_name,
      meanlog = 0
    )

    expect_error(
      primary_censored_cdf(obj_lnorm_no_sdlog, q = 1, pwindow = 1),
      "sdlog parameter is required for Log-Normal distribution"
    )
  }
)

test_that(
  "primary_censored_cdf.default computes the same values as
   primary_censored_cdf.pcens_pgamma_dunif",
  {
    pdist_name <- "pgamma"
    pdist <- pgamma
    dprimary_name <- "dunif"
    dprimary <- dunif

    shapes <- c(0.5, 1, 2, 5)
    rates <- c(0.1, 0.5, 1, 2)
    pwindows <- c(1, 2, 5, 10)

    for (shape in shapes) {
      for (rate in rates) {
        for (pwindow in pwindows) {
          obj <- new_primary_censored_dist(
            pdist,
            dprimary, list(),
            pdist_name, dprimary_name,
            shape = shape, rate = rate
          )

          q_values <- seq(0, 30, by = 0.1)
          result_numeric <- primary_censored_cdf(
            obj,
            q = q_values, pwindow = pwindow, use_numeric = TRUE
          )
          result_analytical <- primary_censored_cdf(
            obj,
            q = q_values, pwindow = pwindow, use_numeric = FALSE
          )

          # Check properties of numeric result
          expect_type(result_numeric, "double")
          expect_length(result_numeric, length(q_values))
          expect_true(
            all(diff(result_numeric) >= 0)
          ) # Ensure CDF is non-decreasing

          # Check that analytical and numeric results are the same
          expect_equal(
            result_numeric, result_analytical,
            tolerance = 1e-5,
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
  "primary_censored_cdf.default computes the same values as
   primary_censored_cdf.pcens_plnorm_dunif",
  {
    pdist_name <- "plnorm"
    pdist <- plnorm
    dprimary_name <- "dunif"
    dprimary <- dunif

    meanlogs <- c(-1, 0, 1, 2)
    sdlogs <- c(0.5, 1, 1.5)
    pwindows <- c(1, 2, 5, 8)

    for (meanlog in meanlogs) {
      for (sdlog in sdlogs) {
        for (pwindow in pwindows) {
          obj <- new_primary_censored_dist(
            pdist,
            dprimary, list(),
            pdist_name, dprimary_name,
            meanlog = meanlog, sdlog = sdlog
          )

          q_values <- seq(0, 30, by = 0.1)
          result_numeric <- primary_censored_cdf(
            obj,
            q = q_values, pwindow = pwindow, use_numeric = TRUE
          )
          result_analytical <- primary_censored_cdf(
            obj,
            q = q_values, pwindow = pwindow, use_numeric = FALSE
          )
          # Check properties of numeric result
          expect_type(result_numeric, "double")
          expect_length(result_numeric, length(q_values))
          expect_true(
            all(diff(result_numeric) >= 0)
          ) # Ensure CDF is non-decreasing

          # Check that analytical and numeric results are the same
          expect_equal(
            result_numeric, result_analytical,
            tolerance = 1e-5,
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
  "primary_censored_cdf.default computes the same values as
   primary_censored_cdf.pcens_pweibull_dunif",
  {
    pdist_name <- "pweibull"
    pdist <- pweibull
    dprimary_name <- "dunif"
    dprimary <- dunif

    shapes <- c(0.5, 1, 2, 3)
    scales <- c(0.5, 1, 2, 5)
    pwindows <- c(1, 2, 5, 10)

    for (shape in shapes) {
      for (scale in scales) {
        for (pwindow in pwindows) {
          obj <- new_primary_censored_dist(
            pdist,
            dprimary, list(),
            pdist_name, dprimary_name,
            shape = shape, scale = scale
          )

          q_values <- seq(0, 30, by = 0.1)
          result_numeric <- primary_censored_cdf(
            obj,
            q = q_values, pwindow = pwindow, use_numeric = TRUE
          )
          result_analytical <- primary_censored_cdf(
            obj,
            q = q_values, pwindow = pwindow, use_numeric = FALSE
          )

          # Check properties of numeric result
          expect_type(result_numeric, "double")
          expect_length(result_numeric, length(q_values))
          expect_true(
            all(diff(result_numeric) >= 0)
          ) # Ensure CDF is non-decreasing

          # Check that analytical and numeric results are the same
          expect_equal(
            result_numeric, result_analytical,
            tolerance = 1e-5,
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
