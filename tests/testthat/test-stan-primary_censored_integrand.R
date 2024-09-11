skip_on_cran()
if (!on_ci()) {
  skip_on_os("windows")
  skip_on_os("mac")
}


test_that("Stan primary_censored_integrand produces expected output", {
  x <- 4.5
  xc <- 0.5
  theta <- c(1.0, 0.5)  # Example parameters for lognormal distribution
  x_r <- c(5.0, 1.0)  # d and pwindow
  # dist_id, primary_dist_id, dist_params_len, primary_params_len
  x_i <- c(1, 1, 2, 0)

  result <- primary_censored_integrand(x, xc, theta, x_r, x_i)

  expect_type(result, "double")
  expect_false(is.na(result))
  expect_true(is.finite(result))
  expect_gt(result, 0)
})

test_that("Stan primary_censored_integrand handles edge cases", {
  theta <- c(1.0, 0.5)  # Example parameters for lognormal distribution
  x_r <- c(1.0, 1.0)  # d and pwindow
  # dist_id, primary_dist_id, dist_params_len, primary_params_len
  x_i <- c(1, 1, 2, 0)

  result_near_zero <- primary_censored_integrand(1e-10, 1e-10, theta, x_r, x_i)
  expect_true(is.finite(result_near_zero))

  result_at_d <- primary_censored_integrand(1, 0.0, theta, x_r, x_i)
  expect_true(is.finite(result_at_d))
  expect_gt(result_at_d, result_near_zero)

  result_near_d <- primary_censored_integrand(
    0.99999, 0.00001, theta, x_r, x_i
  )
  expect_true(is.finite(result_near_d))
  expect_gt(result_at_d, result_near_d)
})

test_that("Stan primary_censored_integrand is continuous", {
  theta <- c(1.5, 0.75)  # Example parameters for lognormal distribution
  x_r <- c(5.0, 5.0)  # d and pwindow
  # dist_id, primary_dist_id, dist_params_len, primary_params_len
  x_i <- c(1, 1, 2, 0)

  x_values <- seq(0.1, 4.9, by = 0.1)
  results <- vapply(
    x_values,
    primary_censored_integrand,
    xc = 0.1,
    theta = theta,
    x_r = x_r,
    x_i = x_i,
    numeric(1)
  )

  # Check if the function is continuous (no sudden jumps)
  diffs <- diff(results)
  expect_lt(max(abs(diffs)), 1e-2)
})

test_that("Stan primary_censored_integrand handles different distributions", {
  x <- 4.5
  xc <- 0.5
  x_r <- c(5.0, 1.0)  # d and pwindow

  # Test for lognormal distribution
  theta_lognormal <- c(1.0, 0.5)
  x_i_lognormal <- c(1, 1, 2, 0)
  result_lognormal <- primary_censored_integrand(
    x, xc, theta_lognormal, x_r, x_i_lognormal
  )
  expect_true(is.finite(result_lognormal))

  # Test for gamma distribution
  theta_gamma <- c(2.0, 1.0)
  x_i_gamma <- c(2, 1, 2, 0)
  result_gamma <- primary_censored_integrand(
    x, xc, theta_gamma, x_r, x_i_gamma
  )
  expect_true(is.finite(result_gamma))

  # Test for weibull distribution
  theta_weibull <- c(1.5, 2.0)
  x_i_weibull <- c(5, 1, 2, 0)
  result_weibull <- primary_censored_integrand(
    x, xc, theta_weibull, x_r, x_i_weibull
  )
  expect_true(is.finite(result_weibull))
})

test_that("Stan primary_censored_integrand handles extreme parameter values", {
  x <-4.5
  xc <- 0.5
  x_r <- c(5.0, 1.0)  # d and pwindow
  # dist_id, primary_dist_id, dist_params_len, primary_params_le
  x_i <- c(1, 1, 2, 0)
  # Test with very small scale parameter
  theta_small_scale <- c(1.0, 1e-10)
  result_small_scale <- primary_censored_integrand(
    x, xc, theta_small_scale, x_r, x_i
  )
  expect_true(is.finite(result_small_scale))

  # Test with very large scale parameter
  theta_large_scale <- c(1.0, 1e10)
  result_large_scale <- primary_censored_integrand(
    x, xc, theta_large_scale, x_r, x_i
  )
  expect_true(is.finite(result_large_scale))

  # Test with very small shape parameter (for distributions that use it)
  x_i_weibull <- c(5, 1, 2, 0)
  theta_small_shape <- c(1e-10, 1.0)
  result_small_shape <- primary_censored_integrand(
    x, xc, theta_small_shape, x_r, x_i_weibull
  )
  expect_true(is.finite(result_small_shape))

  # Test with very large shape parameter
  theta_large_shape <- c(1e10, 1.0)
  result_large_shape <- primary_censored_integrand(
    x, xc, theta_large_shape, x_r, x_i_weibull
  )
  expect_true(is.finite(result_large_shape))
})
