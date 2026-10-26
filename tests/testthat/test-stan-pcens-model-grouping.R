skip_on_cran()
skip_if_not_installed("cmdstanr")
skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)))

# pcens_model.stan takes the log PMF of rows with swindow = 1 from one
# primarycensored_sone_lpmf_vectorized() call per (pwindow, L, D) group. These
# tests compare its log density and gradient with a reference model that calls
# primarycensored_lpmf() for every row.

reference_model <- function() {
  functions <- pcd_load_stan_functions(
    wrap_in_block = TRUE, write_to_file = FALSE
  )
  code <- paste0(
    functions, "\n",
    "data {\n",
    "  int N;\n",
    "  array[N] int d;\n",
    "  array[N] int d_upper;\n",
    "  array[N] int n;\n",
    "  array[N] int pwindow;\n",
    "  array[N] real L;\n",
    "  array[N] real D;\n",
    "  int dist_id;\n",
    "  int primary_id;\n",
    "  int n_params;\n",
    "  int n_primary_params;\n",
    "  vector[n_params] param_lower_bounds;\n",
    "  vector[n_params] param_upper_bounds;\n",
    "  vector[n_primary_params] primary_param_lower_bounds;\n",
    "  vector[n_primary_params] primary_param_upper_bounds;\n",
    "  vector[n_params] prior_location;\n",
    "  vector[n_params] prior_scale;\n",
    "  vector[n_primary_params] primary_prior_location;\n",
    "  vector[n_primary_params] primary_prior_scale;\n",
    "}\n",
    "parameters {\n",
    "  vector<lower=param_lower_bounds, upper=param_upper_bounds>",
    "[n_params] params;\n",
    "  vector<lower=primary_param_lower_bounds,",
    " upper=primary_param_upper_bounds>[n_primary_params] primary_params;\n",
    "}\n",
    "model {\n",
    "  params ~ normal(prior_location, prior_scale);\n",
    "  primary_params ~ normal(primary_prior_location, primary_prior_scale);\n",
    "  for (i in 1:N) {\n",
    "    target += n[i] * primarycensored_lpmf(\n",
    "      d[i] | dist_id, to_array_1d(params), pwindow[i], d_upper[i],\n",
    "      L[i], D[i], primary_id, to_array_1d(primary_params)\n",
    "    );\n",
    "  }\n",
    "}\n"
  )
  path <- file.path(tempdir(), "pcd_reference_model.stan")
  writeLines(code, path)
  suppressMessages(suppressWarnings(
    cmdstanr::cmdstan_model(path, compile_model_methods = TRUE)
  ))
}

# Mixed swindows, pwindows, lower and upper truncation, plus a sparse group
# of one long delay that stays on the per-row path.
simulate_rows <- function(n_rows) {
  set.seed(101)
  delays <- data.frame(
    swindow = sample(c(1, 1, 1, 2), n_rows, replace = TRUE),
    pwindow = sample(c(1, 2), n_rows, replace = TRUE),
    relative_obs_time = sample(c(Inf, 20), n_rows, replace = TRUE),
    start_relative_obs_time = sample(c(-Inf, 0, 2), n_rows, replace = TRUE),
    n = sample.int(3, n_rows, replace = TRUE)
  )
  delays$delay <- pmax(
    pmin(floor(rlnorm(n_rows, 1.2, 0.5)), 15),
    pmax(delays$start_relative_obs_time, 0)
  )
  delays <- rbind(delays, data.frame(
    swindow = 1, pwindow = 1, relative_obs_time = 60,
    start_relative_obs_time = 0, n = 1, delay = 40
  ))
  delays$delay_upper <- delays$delay + delays$swindow
  delays
}

test_that("pcens_model matches per-row primarycensored_lpmf in log density
  and gradient", {
  model <- suppressMessages(suppressWarnings(
    pcd_cmdstan_model(compile_model_methods = TRUE, force_recompile = TRUE)
  ))
  reference <- reference_model()
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
    info <- paste("dist", case$dist_id, "primary", case$primary_id)
    n_primary <- length(case$primary_params)
    stan_data <- pcd_as_stan_data(
      delays,
      dist_id = case$dist_id,
      primary_id = case$primary_id,
      param_bounds = list(lower = c(0.01, 0.01), upper = c(50, 50)),
      primary_param_bounds = list(
        lower = rep(0, n_primary), upper = rep(10, n_primary)
      ),
      priors = list(location = c(1, 1), scale = c(1, 1)),
      primary_priors = list(
        location = rep(0.1, n_primary), scale = rep(1, n_primary)
      ),
      use_reduce_sum = case$use_reduce_sum,
      truncation_check_multiplier = NULL
    )
    init <- list(
      params = case$params, primary_params = as.array(case$primary_params)
    )
    grouped <- stan_log_prob_and_grad( # nolint: object_usage_linter.
      model, stan_data,
      c(init, list(
        np_pmf = as.array(1), np_alpha = numeric(0), np_sigma = numeric(0),
        np_eps = numeric(0)
      ))
    )
    per_row <- stan_log_prob_and_grad( # nolint: object_usage_linter.
      reference, stan_data, init
    )
    expect_equal(
      grouped$log_prob, per_row$log_prob,
      tolerance = 1e-10, info = info
    )
    expect_equal(grouped$grad, per_row$grad, tolerance = 1e-8, info = info)
  }
})
