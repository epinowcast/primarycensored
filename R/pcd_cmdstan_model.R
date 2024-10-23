#' Create a CmdStanModel with primarycensored Stan functions
#'
#' This function creates a CmdStanModel object using the Stan model and
#' functions from primarycensored and optionally includes additional
#' user-specified Stan files.
#'
#' @param include_paths Character vector of paths to include for Stan
#'  compilation. Defaults to the result of `pcd_stan_path()`.
#'
#' @param ... Additional arguments passed to cmdstanr::cmdstan_model().
#'
#' @return A CmdStanModel object.
#'
#' @details
#' The underlying Stan model (`pcens_model.stan`) supports various features:
#' - Multiple probability distributions for modeling delays
#' - Primary and secondary censoring
#' - Truncation
#' - Optional use of reduce_sum for improved performance (via
#'   within chain parallelism).
#' - Flexible prior specifications
#' - Optional computation of log-likelihood for model comparison
#'
#' @export
#' @family modelhelpers
#'
#' @examplesIf requireNamespace("cmdstanr", quietly = TRUE)
#' if (!is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
#'   model <- pcd_cmdstan_model(compile = FALSE)
#'   model
#' }
pcd_cmdstan_model <- function(
    include_paths = primarycensored::pcd_stan_path(),
    ...) {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("Package 'cmdstanr' is required but not installed for this function.")
  }

  pcd_stan_model <- system.file(
    "stan", "pcens_model.stan",
    package = "primarycensored"
  )

  cmdstanr::cmdstan_model(
    pcd_stan_model,
    include_paths = include_paths,
    ...
  )
}

#' Prepare data for primarycensored Stan model
#'
#' This function takes in delay data and prepares it for use with the
#' primarycensored Stan model.
#'
#' @param data A data frame containing the delay data.
#'
#' @param delay Column name for observed delays (default: "delay")
#'
#' @param delay_upper Column name for upper bound of delays
#' (default: "delay_upper")
#'
#' @param n Column name for count of observations (default: "n")
#'
#' @param pwindow Column name for primary window (default: "pwindow")
#'
#' @param relative_obs_time Column name for relative observation time
#' (default: "relative_obs_time")
#'
#' @param dist_id Integer identifying the delay distribution:
#'   1 = Lognormal, 2 = Gamma, 3 = Weibull, 4 = Exponential,
#'   5 = Generalized Gamma, 6 = Negative Binomial, 7 = Poisson,
#'   8 = Bernoulli, 9 = Beta, 10 = Binomial, 11 = Categorical, 12 = Cauchy,
#'   13 = Chi-square, 14 = Dirichlet, 15 = Gumbel, 16 = Inverse Gamma,
#'   17 = Logistic
#'
#' @param primary_id Integer identifying the primary distribution:
#'   1 = Uniform, 2 = Exponential growth
#'
#' @param param_bounds A list with elements `lower` and `upper`, each a numeric
#'   vector specifying bounds for the delay distribution parameters.
#'
#' @param primary_param_bounds A list with elements `lower` and `upper`, each a
#'   numeric vector specifying bounds for the primary distribution parameters.
#'
#' @param priors A list with elements `location` and `scale`, each a numeric
#'   vector specifying priors for the delay distribution parameters.
#'
#' @param primary_priors A list with elements `location` and `scale`, each a
#'   numeric vector specifying priors for the primary distribution parameters.
#'
#' @param compute_log_lik Logical; compute log likelihood? (default: FALSE)
#'
#' @param use_reduce_sum Logical; use reduce_sum for performance?
#'  (default: FALSE)
#'
#' @param truncation_check_multiplier Numeric multiplier to use for checking
#'   if the truncation time D is appropriate relative to the maximum delay
#'   for each unique D value. Set to NULL to skip the check. Default is 2.
#'
#' @return A list containing the data formatted for use with
#' [pcd_cmdstan_model()]
#'
#' @export
#' @family modelhelpers
#'
#' @examples
#' data <- data.frame(
#'   delay = c(1, 2, 3),
#'   delay_upper = c(2, 3, 4),
#'   n = c(10, 20, 15),
#'   pwindow = c(1, 1, 2),
#'   relative_obs_time = c(10, 10, 10)
#' )
#' stan_data <- pcd_as_stan_data(
#'   data,
#'   dist_id = 1,
#'   primary_id = 1,
#'   param_bounds = list(lower = c(0, 0), upper = c(10, 10)),
#'   primary_param_bounds = list(lower = numeric(0), upper = numeric(0)),
#'   priors = list(location = c(1, 1), scale = c(1, 1)),
#'   primary_priors = list(location = numeric(0), scale = numeric(0))
#' )
pcd_as_stan_data <- function(
    data, delay = "delay", delay_upper = "delay_upper",
    n = "n", pwindow = "pwindow",
    relative_obs_time = "relative_obs_time",
    dist_id, primary_id,
    param_bounds, primary_param_bounds,
    priors, primary_priors,
    compute_log_lik = FALSE,
    use_reduce_sum = FALSE,
    truncation_check_multiplier = 2) {
  required_cols <- c(delay, delay_upper, n, pwindow, relative_obs_time)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns: ", toString(missing_cols), "\n",
      "Please ensure your data frame contains these columns or set the",
      " corresponding arguments:\n",
      "delay = '", delay, "'\n",
      "delay_upper = '", delay_upper, "'\n",
      "n = '", n, "'\n",
      "pwindow = '", pwindow, "'\n",
      "relative_obs_time = '", relative_obs_time, "'"
    )
  }

  if (!is.null(truncation_check_multiplier)) {
    unique_D <- unique(data[[relative_obs_time]])
    for (D in unique_D) {
      delays_subset <- data[[delay]][data[[relative_obs_time]] == D]
      check_truncation(
        delays = delays_subset,
        D = D,
        multiplier = truncation_check_multiplier
      )
    }
  }

  stan_data <- list(
    N = nrow(data),
    d = data[[delay]],
    d_upper = data[[delay_upper]],
    n = data[[n]],
    pwindow = data[[pwindow]],
    D = data[[relative_obs_time]],
    dist_id = dist_id,
    primary_id = primary_id,
    n_params = length(param_bounds$lower),
    n_primary_params = length(primary_param_bounds$lower),
    compute_log_lik = as.integer(compute_log_lik),
    use_reduce_sum = as.integer(use_reduce_sum),
    param_lower_bounds = param_bounds$lower,
    param_upper_bounds = param_bounds$upper,
    primary_param_lower_bounds = primary_param_bounds$lower,
    primary_param_upper_bounds = primary_param_bounds$upper,
    prior_location = priors$location,
    prior_scale = priors$scale,
    primary_prior_location = primary_priors$location,
    primary_prior_scale = primary_priors$scale
  )

  return(stan_data)
}
