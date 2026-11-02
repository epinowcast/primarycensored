skip_on_cran()
skip_if_not_installed("cmdstanr")
skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)))

# primarycensored_grouped_lpmf() computes the log CDF once per unique point
# of a group. These tests check it, and the grouped likelihood of
# pcens_model.stan, against per-row primarycensored_lpmf() calls.

# Delays with repeats, secondary windows of 1, 2 and 3, and counts
group_rows <- data.frame(
  d = c(0, 0, 1, 2, 2, 2, 3, 5, 8, 1, 4),
  swindow = c(1, 1, 1, 1, 1, 2, 3, 1, 2, 2, 1),
  n = c(3, 2, 4, 1, 5, 2, 1, 2, 1, 3, 2)
)

per_row_lpmf <- function(rows, dist_id, params, pwindow, L, D,
                         primary_id = 1L, primary_params = numeric(0)) {
  lpmf <- mapply(
    function(d, d_upper) {
      primarycensored_lpmf( # nolint: object_usage_linter.
        d, dist_id, params, pwindow, d_upper, L, D, primary_id,
        primary_params
      )
    },
    rows$d, rows$d + rows$swindow
  )
  sum(rows$n * lpmf)
}

grouped_lpmf <- function(rows, dist_id, params, pwindow, L, D,
                         primary_id = 1L, primary_params = numeric(0)) {
  n_rows <- nrow(rows)
  groups <- .pcd_stan_groups( # nolint: object_usage_linter.
    rows$d, rows$d + rows$swindow, rows$n, rep(pwindow, n_rows),
    rep(L, n_rows), rep(D, n_rows)
  )
  primarycensored_grouped_lpmf( # nolint: object_usage_linter.
    groups$interval_n, groups$points, groups$interval_lower,
    groups$interval_upper, dist_id, params, pwindow, L, D, primary_id,
    primary_params
  )
}

test_that("primarycensored_grouped_lpmf matches summed primarycensored_lpmf", {
  cases <- list(
    list(dist_id = 1L, params = c(1.2, 0.6)),
    list(dist_id = 2L, params = c(2, 0.5)),
    list(dist_id = 3L, params = c(1.5, 4)),
    list(dist_id = 1L, params = c(1.2, 0.6), primary_id = 2L,
         primary_params = 0.2)
  )
  bounds <- list(c(-Inf, Inf), c(0, Inf), c(0, 12), c(-Inf, 10.5))
  for (case in cases) {
    primary_id <- if (is.null(case$primary_id)) 1L else case$primary_id
    for (pwindow in c(1, 2)) {
      for (b in bounds) {
        args <- list(
          group_rows, case$dist_id, case$params, pwindow, b[1], b[2],
          primary_id, as.numeric(case$primary_params)
        )
        expect_equal(
          do.call(grouped_lpmf, args), do.call(per_row_lpmf, args),
          tolerance = 1e-10,
          info = paste("dist", case$dist_id, "primary", primary_id,
                       "pwindow", pwindow, "L", b[1], "D", b[2])
        )
      }
    }
  }
})

test_that("primarycensored_grouped_lpmf matches for a delay on the reals", {
  rows <- data.frame(d = c(-2, -1, -1, 0, 3), swindow = c(1, 1, 1, 2, 1),
                     n = c(1, 2, 3, 4, 5))
  for (b in list(c(-Inf, Inf), c(-3, 6))) {
    args <- list(rows, 18L, c(0.5, 2), 1, b[1], b[2])
    expect_equal(
      do.call(grouped_lpmf, args), do.call(per_row_lpmf, args),
      tolerance = 1e-10
    )
  }
})

test_that("primarycensored_grouped_lpmf rejects points outside [L, D]", {
  expect_error(primarycensored_grouped_lpmf( # nolint: object_usage_linter.
    1L, c(0, 1), 1L, 2L, 1L, c(1, 0.5), 1, 0, 0.5, 1L, numeric(0)
  ))
})

# Mixed pwindows, secondary windows, L and D with heavy repetition
simulate_delays <- function(n_rows) {
  set.seed(101)
  delays <- data.frame(
    swindow = sample(c(1, 1, 1, 2), n_rows, replace = TRUE),
    pwindow = sample(c(1, 2), n_rows, replace = TRUE),
    relative_obs_time = sample(c(Inf, 20, 25.5), n_rows, replace = TRUE),
    start_relative_obs_time = sample(c(-Inf, 0, 1.5), n_rows, replace = TRUE),
    n = sample.int(3, n_rows, replace = TRUE)
  )
  delays$delay <- pmax(
    pmin(floor(rlnorm(n_rows, 1.2, 0.5)), 15),
    ceiling(pmax(delays$start_relative_obs_time, 0))
  )
  delays$delay_upper <- delays$delay + delays$swindow
  delays
}

# The package model with model methods, compiled from a copy so it does not
# replace the executable used by pcd_cmdstan_model()
model_path <- file.path(tempdir(), "pcensgrouped.stan")
file.copy(
  system.file("stan", "pcens_model.stan", package = "primarycensored"),
  model_path,
  overwrite = TRUE
)
grouped_model <- suppressMessages(suppressWarnings(cmdstanr::cmdstan_model(
  model_path,
  include_paths = pcd_stan_path(), compile_model_methods = TRUE
)))

expect_grouped_model_matches <- function(delays, case) {
  info <- paste(
    "dist", case$dist_id, "primary", case$primary_id,
    "reduce_sum", case$use_reduce_sum
  )
  n_params <- length(case$params)
  n_primary <- length(case$primary_params)
  stan_data <- pcd_as_stan_data( # nolint: object_usage_linter.
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
  grouped <- stan_log_prob_and_grad(grouped_model, stan_data, init)
  stan_data$use_grouped <- 0L
  per_row <- stan_log_prob_and_grad(grouped_model, stan_data, init)
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

test_that("pcens_model grouped likelihood matches the per-row likelihood", {
  delays <- simulate_delays(200)
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
    expect_grouped_model_matches(delays, case)
  }
})

test_that("pcens_model grouped likelihood matches for a non-parametric
  delay", {
  delays <- data.frame(
    delay = c(2, 3, 2, 3), delay_upper = c(3, 4, 3, 5), n = c(5, 3, 2, 1),
    pwindow = 1, relative_obs_time = c(5, 5, 6, 6)
  )
  expect_grouped_model_matches(delays, list(
    dist_id = 26L, primary_id = 1L, params = numeric(0),
    primary_params = numeric(0), use_reduce_sum = FALSE,
    dist_options = list(K = 6, boundaries = 0:6),
    np_pmf = c(0.1, 0.2, 0.3, 0.2, 0.1, 0.1)
  ))
})
