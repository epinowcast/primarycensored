test_that("pcens_cdf methods dispatch correctly to existing
   analytical solutions", {
  pdist <- pgamma
  dprimary <- dunif

  obj_gamma <- new_pcens(
    pdist,
    dprimary,
    list(),
    shape = 2,
    rate = 1
  )

  pdist <- plnorm
  dprimary <- dunif

  obj_lnorm <- new_pcens(
    pdist,
    dprimary,
    list(),
    meanlog = 0,
    sdlog = 1
  )

  pdist <- pweibull
  dprimary <- dunif

  obj_weibull <- new_pcens(
    pdist,
    dprimary,
    list(),
    shape = 2,
    scale = 1
  )

  expect_s3_class(obj_gamma, "pcens_pgamma_dunif")
  expect_s3_class(obj_lnorm, "pcens_plnorm_dunif")
  expect_s3_class(obj_weibull, "pcens_pweibull_dunif")

  q_values <- c(5, 10)
  pwindow <- 2

  expect_no_error(
    pcens_cdf(obj_gamma, q = q_values, pwindow = pwindow)
  )
  expect_no_error(
    pcens_cdf(obj_lnorm, q = q_values, pwindow = pwindow)
  )
  expect_no_error(
    pcens_cdf(obj_weibull, q = q_values, pwindow = pwindow)
  )
})

test_that("pcens_cdf errors as expected when the wrong distributional
   parameters are supplied", {
  pdist <- pgamma
  dprimary <- dunif

  obj_gamma <- new_pcens(
    pdist,
    dprimary,
    list(),
    rate = 1
  )

  expect_error(
    pcens_cdf(obj_gamma, q = 1, pwindow = 1),
    "shape parameter is required for Gamma distribution"
  )

  obj_gamma_no_rate <- new_pcens(
    pdist,
    dprimary,
    list(),
    shape = 2
  )

  expect_error(
    pcens_cdf(obj_gamma_no_rate, q = 1, pwindow = 1),
    "scale or rate parameter is required for Gamma distribution"
  )

  pdist <- plnorm

  obj_lnorm_no_meanlog <- new_pcens(
    pdist,
    dprimary,
    list(),
    sdlog = 1
  )

  expect_error(
    pcens_cdf(obj_lnorm_no_meanlog, q = 1, pwindow = 1),
    "meanlog parameter is required for Log-Normal distribution"
  )

  obj_lnorm_no_sdlog <- new_pcens(
    pdist,
    dprimary,
    list(),
    meanlog = 0
  )

  expect_error(
    pcens_cdf(obj_lnorm_no_sdlog, q = 1, pwindow = 1),
    "sdlog parameter is required for Log-Normal distribution"
  )

  pdist <- pweibull

  obj_weibull_no_shape <- new_pcens(
    pdist,
    dprimary,
    list(),
    scale = 1
  )

  expect_error(
    pcens_cdf(obj_weibull_no_shape, q = 1, pwindow = 1),
    "shape parameter is required for Weibull distribution"
  )

  obj_weibull_no_scale <- new_pcens(
    pdist,
    dprimary,
    list(),
    shape = 2
  )

  expect_error(
    pcens_cdf(obj_weibull_no_scale, q = 1, pwindow = 1),
    "scale parameter is required for Weibull distribution"
  )
})

test_that("pcens_cdf.default computes the same values as
   pcens_cdf.pcens_pgamma_dunif", {
  pdist <- pgamma
  dprimary <- dunif

  shapes <- c(0.5, 1, 2, 5)
  rates <- c(0.1, 0.5, 1, 2)
  pwindows <- c(1, 2, 5, 10)

  for (shape in shapes) {
    for (rate in rates) {
      for (pwindow in pwindows) {
        obj <- new_pcens(
          pdist,
          dprimary,
          list(),
          shape = shape,
          rate = rate
        )

        q_values <- seq(0, 30, by = 0.1)
        result_numeric <- pcens_cdf(
          obj,
          q = q_values,
          pwindow = pwindow,
          use_numeric = TRUE
        )
        result_analytical <- pcens_cdf(
          obj,
          q = q_values,
          pwindow = pwindow,
          use_numeric = FALSE
        )

        # Check properties of numeric result
        expect_type(result_numeric, "double")
        expect_length(result_numeric, length(q_values))
        expect_true(
          all(diff(result_numeric) >= 0)
        ) # Ensure CDF is non-decreasing

        # Check that analytical and numeric results are the same
        expect_equal(
          result_numeric,
          result_analytical,
          tolerance = 1e-5,
          info = sprintf(
            "Mismatch for shape = %s, rate = %s, pwindow = %s",
            shape,
            rate,
            pwindow
          )
        )
      }
    }
  }
})

test_that("pcens_cdf.default computes the same values as
   pcens_cdf.pcens_plnorm_dunif", {
  pdist <- plnorm
  dprimary <- dunif

  meanlogs <- c(-1, 0, 1, 2)
  sdlogs <- c(0.5, 1, 1.5)
  pwindows <- c(1, 2, 5, 8)

  for (meanlog in meanlogs) {
    for (sdlog in sdlogs) {
      for (pwindow in pwindows) {
        obj <- new_pcens(
          pdist,
          dprimary,
          list(),
          meanlog = meanlog,
          sdlog = sdlog
        )

        q_values <- seq(0, 30, by = 0.1)
        result_numeric <- pcens_cdf(
          obj,
          q = q_values,
          pwindow = pwindow,
          use_numeric = TRUE
        )
        result_analytical <- pcens_cdf(
          obj,
          q = q_values,
          pwindow = pwindow,
          use_numeric = FALSE
        )
        # Check properties of numeric result
        expect_type(result_numeric, "double")
        expect_length(result_numeric, length(q_values))
        expect_true(
          all(diff(result_numeric) >= 0)
        ) # Ensure CDF is non-decreasing

        # Check that analytical and numeric results are the same
        expect_equal(
          result_numeric,
          result_analytical,
          tolerance = 1e-5,
          info = sprintf(
            "Mismatch for meanlog = %s, sdlog = %s, pwindow = %s",
            meanlog,
            sdlog,
            pwindow
          )
        )
      }
    }
  }
})

test_that("pcens_cdf.default computes the same values as
   pcens_cdf.pcens_pweibull_dunif", {
  pdist <- pweibull
  dprimary <- dunif

  shapes <- c(0.5, 1, 2)
  scales <- c(0.5, 1, 2)
  pwindows <- c(1, 2, 3, 4, 5)

  for (shape in shapes) {
    for (scale in scales) {
      for (pwindow in pwindows) {
        obj <- new_pcens(
          pdist,
          dprimary,
          list(),
          shape = shape,
          scale = scale
        )

        q_values <- seq(0, 30, by = 0.1)
        result_numeric <- pcens_cdf(
          obj,
          q = q_values,
          pwindow = pwindow,
          use_numeric = TRUE
        )
        result_analytical <- pcens_cdf(
          obj,
          q = q_values,
          pwindow = pwindow,
          use_numeric = FALSE
        )

        # Check properties of numeric result
        expect_type(result_numeric, "double")
        expect_length(result_numeric, length(q_values))
        expect_true(
          all(diff(result_numeric) >= 0)
        ) # Ensure CDF is non-decreasing

        # Check that analytical and numeric results are the same
        expect_equal(
          result_numeric,
          result_analytical,
          tolerance = 1e-5,
          info = sprintf(
            "Mismatch for shape = %s, scale = %s, pwindow = %s",
            shape,
            scale,
            pwindow
          )
        )
      }
    }
  }
})

test_that("new_pcens *_name deprecation is soft.", {
  pdist <- function(...) pgamma(...)
  dprimary <- function(...) dunif(...)
  shape <- 2
  rate <- 1

  neg_obj <- new_pcens(
    pdist,
    dprimary,
    list(),
    shape = shape,
    rate = rate
  )

  expect_s3_class(neg_obj, "pcens_unknown_unknown")

  ref_obj <- new_pcens(
    add_name_attribute(pdist, "pgamma"),
    add_name_attribute(dprimary, "dunif"),
    list(),
    shape = shape,
    rate = rate
  )

  lifecycle::expect_deprecated(
    obj <- new_pcens(
      # nolint
      pdist,
      add_name_attribute(dprimary, "dunif"),
      list(),
      pdist_name = "pgamma",
      shape = shape,
      rate = rate
    )
  )

  lifecycle::expect_deprecated(
    new_obj <- new_pcens(
      # nolint
      add_name_attribute(pdist, "pgamma"),
      dprimary,
      list(),
      dprimary_name = "dunif",
      shape = shape,
      rate = rate
    )
  )

  expect_identical(body(obj$pdist), body(ref_obj$pdist))
  expect_identical(body(new_obj$pdist), body(ref_obj$pdist))
  expect_identical(formals(obj$pdist), formals(ref_obj$pdist))
  expect_identical(formals(new_obj$pdist), formals(ref_obj$pdist))
  expect_identical(body(obj$dprimary), body(ref_obj$dprimary))
  expect_identical(body(new_obj$dprimary), body(ref_obj$dprimary))
  expect_identical(formals(obj$dprimary), formals(new_obj$dprimary))
  expect_identical(formals(new_obj$dprimary), formals(ref_obj$dprimary))
})

test_that("pcens_cdf returns values in [0, 1] for all analytical methods", {
  # Test case from issue #238: parameters that stress numerical precision
  # Lognormal with tight distribution where CDF approaches 1 quickly
  obj_lnorm <- new_pcens(
    plnorm,
    dunif,
    list(),
    meanlog = 0.55,
    sdlog = 0.27
  )

  q_values <- seq(0, 30, by = 1)
  cdf_lnorm <- pcens_cdf(obj_lnorm, q = q_values, pwindow = 1)
  expect_true(
    all(cdf_lnorm >= 0 & cdf_lnorm <= 1),
    info = "Lognormal CDF should be in [0, 1]"
  )

  # Gamma distribution
  obj_gamma <- new_pcens(
    pgamma,
    dunif,
    list(),
    shape = 2,
    rate = 0.5
  )

  cdf_gamma <- pcens_cdf(obj_gamma, q = q_values, pwindow = 1)
  expect_true(
    all(cdf_gamma >= 0 & cdf_gamma <= 1),
    info = "Gamma CDF should be in [0, 1]"
  )

  # Weibull distribution
  obj_weibull <- new_pcens(
    pweibull,
    dunif,
    list(),
    shape = 2,
    scale = 1
  )

  cdf_weibull <- pcens_cdf(obj_weibull, q = q_values, pwindow = 1)
  expect_true(
    all(cdf_weibull >= 0 & cdf_weibull <= 1),
    info = "Weibull CDF should be in [0, 1]"
  )
})

test_that("new_pcens works with custom function with name attribute", {
  # Create custom functions with name attributes
  custom_pdist <- function(x, shape, rate) pgamma(x, shape, rate)
  custom_dprimary <- function(x, min = 0, max = 1) dunif(x, min, max)

  named_pdist <- add_name_attribute(custom_pdist, "pgamma")
  named_dprimary <- add_name_attribute(custom_dprimary, "dunif")

  # Create pcens object with custom named functions
  obj <- new_pcens(
    named_pdist,
    named_dprimary,
    list(),
    shape = 2,
    rate = 1
  )

  # Check class is set correctly using function names
  expect_s3_class(obj, "pcens_pgamma_dunif")

  # Check functions are preserved correctly
  expect_identical(body(obj$pdist), body(custom_pdist))
  expect_identical(formals(obj$pdist), formals(custom_pdist))
  expect_identical(body(obj$dprimary), body(custom_dprimary))
  expect_identical(formals(obj$dprimary), formals(custom_dprimary))

  # Check arguments are preserved
  expect_identical(obj$args, list(shape = 2, rate = 1))
})
