test_that("new_primary_censored_dist creates object with correct structure", {
  pdist <- pgamma
  dprimary <- dunif
  dprimary_args <- list(min = 0, max = 10)
  shape <- 2
  rate <- 1

  obj <- new_primary_censored_dist(
    pdist, dprimary, dprimary_args,
    shape = shape, rate = rate
  )

  expect_s3_class(obj, "pcens_pgamma_dunif")
  expect_s3_class(obj, "pcens_numeric")
  expect_identical(obj$pdist, pdist)
  expect_identical(obj$dprimary, dprimary)
  expect_identical(obj$dprimary_args, dprimary_args)
  expect_identical(obj$args, list(shape = shape, rate = rate))
})

test_that("primary_censored_cdf.pcens_numeric computes correct values", {
  pdist <- pgamma
  dprimary <- dunif
  dprimary_args <- list(min = 0, max = 10)
  shape <- 2
  rate <- 1

  obj <- new_primary_censored_dist(
    pdist, dprimary, dprimary_args,
    shape = shape, rate = rate
  )

  q_values <- c(0, 5, 10, 15)
  pwindow <- 10

  result <- primary_censored_cdf(obj, q = q_values, pwindow = pwindow)

  expect_true(all(result >= 0))
  expect_true(all(result <= 1))
  expect_type(result, "double")
  expect_length(result, length(q_values))
  expect_true(all(diff(result) >= 0)) # Ensure CDF is non-decreasing
})

test_that("primary_censored_cdf methods dispatch correctly", {
  pdist_gamma <- pgamma
  pdist_lnorm <- plnorm
  dprimary <- dunif
  dprimary_args <- list(min = 0, max = 10)

  obj_gamma <- new_primary_censored_dist(
    pdist_gamma, dprimary, dprimary_args,
    shape = 2, rate = 1
  )
  obj_lnorm <- new_primary_censored_dist(
    pdist_lnorm, dprimary, dprimary_args,
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

test_that("primary_censored_cdf handles edge cases correctly", {
  pdist <- pgamma
  dprimary <- dunif
  dprimary_args <- list(min = 0, max = 10)
  shape <- 2
  rate <- 1

  obj <- new_primary_censored_dist(
    pdist, dprimary, dprimary_args,
    shape = shape, rate = rate
  )

  expect_type(primary_censored_cdf(obj, q = 0, pwindow = 10), "double")
  expect_identical(primary_censored_cdf(obj, q = 0, pwindow = 10), 0)
  expect_identical(primary_censored_cdf(obj, q = -1, pwindow = 10), 0)
  expect_lte(primary_censored_cdf(obj, q = Inf, pwindow = 10), 1)
})
