skip_on_cran()
skip_if_not_installed("cmdstanr")
skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)))

# pcens_model.stan takes the log PMF of grouped rows from one
# primarycensored_sone_lpmf_vectorized() call per group, see
# primarycensored_sone_groups(). These tests check the grouping and compare
# the model's log density and gradient with the same model with no rows
# grouped, where every row calls primarycensored_lpmf().

test_that("primarycensored_sone_groups keeps dense swindow 1 groups", {
  rows <- data.frame(
    d = c(0, 1, 2, 3, 0, 10, 2, 1, 1, 2, -1, 4),
    swindow = c(1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1),
    pwindow = c(1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1),
    L = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1.5, -Inf, 0),
    D = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf, 20, 20, 20, Inf, Inf)
  )
  groups <- primarycensored_sone_groups( # nolint: object_usage_linter.
    as.integer(rows$d), as.integer(rows$d + rows$swindow),
    as.integer(rows$pwindow), rows$L, rows$D
  )
  # Dense: rows 1-4 and 12 (pwindow 1, L 0, D Inf) and rows 8-9 (D 20)
  expect_length(unique(groups[c(1:4, 12)]), 1)
  expect_length(unique(groups[8:9]), 1)
  expect_setequal(groups[c(1, 8)], 1:2)
  # Sparse (rows 5-6 cover 2 of 11 delays) or a group of one delay above 0
  # (row 10), another swindow (row 7) and a negative delay (row 11)
  expect_identical(groups[c(5:7, 10:11)], rep(0L, 5))
})

# Mixed swindows, pwindows, integer and non-integer L and D, plus a sparse
# group of one long delay that stays on the per-row path.
simulate_rows <- function(n_rows) {
  set.seed(101)
  delays <- data.frame(
    swindow = sample(c(1, 1, 1, 2), n_rows, replace = TRUE),
    pwindow = sample(c(1, 2), n_rows, replace = TRUE),
    relative_obs_time = sample(c(Inf, 20, 25.5), n_rows, replace = TRUE),
    start_relative_obs_time = sample(
      c(-Inf, 0, 1.5), n_rows,
      replace = TRUE
    ),
    n = sample.int(3, n_rows, replace = TRUE)
  )
  delays$delay <- pmax(
    pmin(floor(rlnorm(n_rows, 1.2, 0.5)), 15),
    ceiling(pmax(delays$start_relative_obs_time, 0))
  )
  delays <- rbind(delays, data.frame(
    swindow = 1, pwindow = 1, relative_obs_time = 60,
    start_relative_obs_time = 0, n = 1, delay = 40
  ))
  delays$delay_upper <- delays$delay + delays$swindow
  delays
}

expect_grouped_matches_per_row <- function(delays, case) {
  info <- paste(
    "dist", case$dist_id, "primary", case$primary_id,
    "reduce_sum", case$use_reduce_sum
  )
  n_params <- length(case$params)
  n_primary <- length(case$primary_params)
  stan_data <- pcd_as_stan_data(
    delays,
    dist_id = case$dist_id,
    primary_id = case$primary_id,
    param_bounds = list(
      lower = rep(0.01, n_params), upper = rep(50, n_params)
    ),
    primary_param_bounds = list(
      lower = rep(0, n_primary), upper = rep(10, n_primary)
    ),
    priors = list(location = rep(1, n_params), scale = rep(1, n_params)),
    primary_priors = list(
      location = rep(0.1, n_primary), scale = rep(1, n_primary)
    ),
    use_reduce_sum = case$use_reduce_sum,
    truncation_check_multiplier = NULL,
    dist_options = case$dist_options
  )
  init <- list(
    params = as.array(case$params),
    primary_params = as.array(case$primary_params),
    np_pmf = as.array(if (is.null(case$np_pmf)) 1 else case$np_pmf),
    np_alpha = numeric(0), np_sigma = numeric(0), np_eps = numeric(0)
  )
  # nolint start: object_usage_linter.
  grouped <- stan_log_prob_and_grad(compile_pcens_model(), stan_data, init)
  per_row <- stan_log_prob_and_grad(
    compile_pcens_model(grouped = FALSE), stan_data, init
  )
  # nolint end
  testthat::expect_true(is.finite(grouped$log_prob), info = info)
  testthat::expect_equal(
    grouped$log_prob, per_row$log_prob,
    tolerance = 1e-10, info = info
  )
  testthat::expect_equal(
    grouped$grad, per_row$grad,
    tolerance = 1e-8, info = info
  )
}

test_that("pcens_model matches per-row primarycensored_lpmf in log density
  and gradient", {
  delays <- simulate_rows(200)
  cases <- list(
    list(
      dist_id = 1L, primary_id = 1L, params = c(1.3, 0.6),
      primary_params = numeric(0), use_reduce_sum = FALSE
    ),
    list(
      dist_id = 2L, primary_id = 1L, params = c(2, 0.5),
      primary_params = numeric(0), use_reduce_sum = TRUE
    ),
    list(
      dist_id = 3L, primary_id = 2L, params = c(1.5, 5),
      primary_params = 0.2, use_reduce_sum = FALSE
    )
  )
  for (case in cases) {
    expect_grouped_matches_per_row(delays, case)
  }
})

test_that("pcens_model matches per-row primarycensored_lpmf when every row
  is grouped", {
  delays <- data.frame(
    delay = 0:12, delay_upper = 1:13, n = 13:1, pwindow = 1,
    relative_obs_time = 20
  )
  for (use_reduce_sum in c(FALSE, TRUE)) {
    expect_grouped_matches_per_row(delays, list(
      dist_id = 2L, primary_id = 1L, params = c(2, 0.5),
      primary_params = numeric(0), use_reduce_sum = use_reduce_sum
    ))
  }
})

test_that("pcens_model matches per-row primarycensored_lpmf for a
  non-parametric delay", {
  delays <- data.frame(
    delay = c(2, 3, 2, 3), delay_upper = c(3, 4, 3, 4), n = c(5, 3, 2, 1),
    pwindow = 1, relative_obs_time = c(4, 4, 6, 6)
  )
  expect_grouped_matches_per_row(delays, list(
    dist_id = 26L, primary_id = 1L, params = numeric(0),
    primary_params = numeric(0), use_reduce_sum = FALSE,
    dist_options = list(K = 6, boundaries = 0:6),
    np_pmf = c(0.1, 0.2, 0.3, 0.2, 0.1, 0.1)
  ))
})
