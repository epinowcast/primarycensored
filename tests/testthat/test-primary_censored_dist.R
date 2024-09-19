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

test_that("primary_censored_cdf.default computes correct values", {
  pdist_name <- "pgamma"
  pdist <- pgamma
  dprimary_name <- "dunif"
  dprimary <- dunif
  shape <- 2
  rate <- 0.5

  obj <- new_primary_censored_dist(
    pdist,
    dprimary, list(),
    pdist_name, dprimary_name,
    shape = shape, rate = rate
  )

  q_values <- 0:20
  pwindow <- 2
  result_numeric <- primary_censored_cdf(
    obj,
    q = q_values, pwindow = pwindow, use_numeric = TRUE
  )
  result_analytical <- primary_censored_cdf(
    obj,
    q = q_values, pwindow = pwindow, use_numeric = FALSE
  )

  # Check properties of numeric result
  expect_true(all(result_numeric >= 0))
  expect_true(all(result_numeric <= 1))
  expect_type(result_numeric, "double")
  expect_length(result_numeric, length(q_values))
  expect_true(all(diff(result_numeric) >= 0)) # Ensure CDF is non-decreasing

  # Check that analytical and numeric results are the same
  expect_equal(result_numeric, result_analytical, tolerance = 1e-6)
})

test_that("primary_censored_cdf methods dispatch correctly", {
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

  expect_s3_class(obj_gamma, "pcens_pgamma_dunif")
  expect_s3_class(obj_lnorm, "pcens_plnorm_dunif")

  q_values <- c(5, 10)
  pwindow <- 10

  expect_no_error(
    primary_censored_cdf(obj_gamma, q = q_values, pwindow = pwindow)
  )
  expect_no_error(
    primary_censored_cdf(obj_lnorm, q = q_values, pwindow = pwindow)
  )
})

test_that("primary_censored_cdf.pcens_pgamma_dunif uses numeric method", {
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

  q_values <- c(5, 10)
  pwindow <- 10

  expect_identical(
    primary_censored_cdf(obj, q = q_values, pwindow = pwindow),
    primary_censored_cdf(
      obj,
      q = q_values, pwindow = pwindow, use_numeric = TRUE
    )
  )
})

test_that("primary_censored_cdf.pcens_plnorm_dunif uses numeric method", {
  pdist_name <- "plnorm"
  pdist <- plnorm
  dprimary_name <- "dunif"
  dprimary <- dunif
  meanlog <- 0
  sdlog <- 1

  obj <- new_primary_censored_dist(
    pdist,
    dprimary, list(),
    pdist_name, dprimary_name,
    meanlog = meanlog, sdlog = sdlog
  )

  q_values <- c(5, 10)
  pwindow <- 10

  expect_identical(
    primary_censored_cdf(obj, q = q_values, pwindow = pwindow),
    primary_censored_cdf(
      obj,
      q = q_values, pwindow = pwindow, use_numeric = TRUE
    )
  )
})
