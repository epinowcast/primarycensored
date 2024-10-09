skip_on_cran()
if (on_ci()) {
  skip_on_os("windows")
  skip_on_os("mac")
}


test_that("Stan primarycensored_ode produces expected output", {
  t <- 4.5
  y <- 0.5
  theta <- c(1.0, 0.5) # Example parameters for lognormal distribution
  x_r <- c(5.0, 1.0) # d and pwindow
  # dist_id, primray_id, dist_params_len, primary_params_len
  x_i <- c(1, 1, 2, 0)

  result <- primarycensored_ode(t, y, theta, x_r, x_i)

  expect_type(result, "double")
  expect_length(result, 1)
  expect_false(is.na(result[1]))
  expect_true(is.finite(result[1]))
  expect_gt(result[1], 0)
})

test_that("Stan primarycensored_ode handles edge cases", {
  theta <- c(1.0, 0.5) # Example parameters for lognormal distribution
  x_r <- c(1.0, 1.0) # d and pwindow
  # dist_id, primray_id, dist_params_len, primary_params_len
  x_i <- c(1, 1, 2, 0)

  result_near_zero <- primarycensored_ode(1e-10, 1e-10, theta, x_r, x_i)
  expect_true(is.finite(result_near_zero[1]))

  result_at_d <- primarycensored_ode(1, 0.0, theta, x_r, x_i)
  expect_true(is.finite(result_at_d[1]))
  expect_gt(result_at_d[1], result_near_zero[1])

  result_near_d <- primarycensored_ode(
    0.99999, 0.00001, theta, x_r, x_i
  )
  expect_true(is.finite(result_near_d[1]))
  expect_gt(result_at_d[1], result_near_d[1])
})

test_that("Stan primarycensored_ode is continuous", {
  theta <- c(1.5, 0.75) # Example parameters for lognormal distribution
  x_r <- c(5.0, 5.0) # d and pwindow
  # dist_id, primray_id, dist_params_len, primary_params_len
  x_i <- c(1, 1, 2, 0)

  t_values <- seq(0.1, 4.9, by = 0.1)
  results <- vapply(
    t_values,
    function(t) primarycensored_ode(t, 0.1, theta, x_r, x_i)[1],
    numeric(1)
  )

  # Check if the function is continuous (no sudden jumps)
  diffs <- diff(results)
  expect_lt(max(abs(diffs)), 1e-2)
})

test_that("Stan primarycensored_ode handles different distributions", {
  t <- 4.5
  y <- 0.5
  x_r <- c(5.0, 1.0) # d and pwindow

  # Test for lognormal distribution
  theta_lognormal <- c(1.0, 0.5)
  x_i_lognormal <- c(1, 1, 2, 0)
  result_lognormal <- primarycensored_ode(
    t, y, theta_lognormal, x_r, x_i_lognormal
  )
  expect_true(is.finite(result_lognormal[1]))

  # Test for gamma distribution
  theta_gamma <- c(2.0, 1.0)
  x_i_gamma <- c(2, 1, 2, 0)
  result_gamma <- primarycensored_ode(
    t, y, theta_gamma, x_r, x_i_gamma
  )
  expect_true(is.finite(result_gamma[1]))

  # Test for weibull distribution
  theta_weibull <- c(1.5, 2.0)
  x_i_weibull <- c(5, 1, 2, 0)
  result_weibull <- primarycensored_ode(
    t, y, theta_weibull, x_r, x_i_weibull
  )
  expect_true(is.finite(result_weibull[1]))
})

test_that("Stan primarycensored_ode handles extreme parameter values", {
  t <- 4.5
  y <- 0.5
  x_r <- c(5.0, 1.0) # d and pwindow
  # dist_id, primray_id, dist_params_len, primary_params_len
  x_i <- c(1, 1, 2, 0)
  # Test with very small scale parameter
  theta_small_scale <- c(1.0, 1e-10)
  result_small_scale <- primarycensored_ode(
    t, y, theta_small_scale, x_r, x_i
  )
  expect_true(is.finite(result_small_scale[1]))

  # Test with very large scale parameter
  theta_large_scale <- c(1.0, 1e10)
  result_large_scale <- primarycensored_ode(
    t, y, theta_large_scale, x_r, x_i
  )
  expect_true(is.finite(result_large_scale[1]))

  # Test with very small shape parameter (for distributions that use it)
  x_i_weibull <- c(5, 1, 2, 0)
  theta_small_shape <- c(1e-10, 1.0)
  result_small_shape <- primarycensored_ode(
    t, y, theta_small_shape, x_r, x_i_weibull
  )
  expect_true(is.finite(result_small_shape[1]))

  # Test with very large shape parameter
  theta_large_shape <- c(1e10, 1.0)
  result_large_shape <- primarycensored_ode(
    t, y, theta_large_shape, x_r, x_i_weibull
  )
  expect_true(is.finite(result_large_shape[1]))
})
