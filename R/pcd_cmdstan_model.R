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
    stop(
      "Package 'cmdstanr' is required but not installed for this function.",
      call. = FALSE
    )
  }

  pcd_stan_model <- system.file(
    "stan",
    "pcens_model.stan",
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
#' @param start_relative_obs_time Column name for start of relative observation
#'  time, used as the lower truncation point L. If the column is not present in
#'  data, L = 0 is assumed for all observations.
#'  (default: "start_relative_obs_time")
#'
#' @param relative_obs_time Column name for relative observation time, used as
#'  the upper truncation point D (default: "relative_obs_time")
#'
#' @param dist_id Integer identifying the delay distribution:
#'   You can use [pcd_stan_dist_id()] to get the dist ID for a
#'   distribution or look at the [pcd_distributions] data set.
#'
#' @param primary_id Integer identifying the primary distribution:
#'   You can use [pcd_stan_dist_id()] to get the primary dist ID for a
#'   distribution (make sure to select the "primary" type) or look at the
#'   [pcd_primary_distributions] data set.
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
#'   (default: FALSE)
#'
#' @param truncation_check_multiplier Numeric multiplier to use for checking
#'   if the truncation time D is appropriate relative to the maximum delay
#'   for each unique D value. Set to NULL to skip the check. Default is 2.
#'
#' @param nonparametric Optional list configuring a non-parametric delay
#'   (`dist_id` 26 or 27). When `NULL` (the default) the model uses the
#'   parametric path and all `np_*` Stan data fields are sized 0. When
#'   supplied, expects a list with elements:
#'   * `K`: integer, number of bins.
#'   * `boundaries`: numeric vector of length `K + 1`.
#'   * `paramtype`: `"simplex"` (Dirichlet on the PMF, `dist_id` 26),
#'     `"hazard"` / `"hazardrw"` (random walk on logit hazards,
#'     `dist_id` 27), or `"hazardre"` (IID logit random effects on the
#'     hazards around a mean, `dist_id` 27 with `np_paramtype = 3`).
#'   * `dirichlet_alpha`: optional numeric vector of length `K`,
#'     defaults to `rep(1, K)`.
#'   * `hazard_priors`: optional list with elements `alpha_mean`,
#'     `alpha_sd`, `log_sigma_mean`, `log_sigma_sd`. Defaults to
#'     `list(alpha_mean = 0, alpha_sd = 5, log_sigma_mean = 0,
#'     log_sigma_sd = 1)`.
#'   When supplied, `param_bounds`, `priors`, and the parametric
#'   `dist_id` argument are ignored (`dist_id` is set from `paramtype`).
#'
#' @return A list containing the data formatted for use with
#'   [pcd_cmdstan_model()]
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
    data,
    delay = "delay",
    delay_upper = "delay_upper",
    n = "n",
    pwindow = "pwindow",
    start_relative_obs_time = "start_relative_obs_time",
    relative_obs_time = "relative_obs_time",
    dist_id,
    primary_id,
    param_bounds,
    primary_param_bounds,
    priors,
    primary_priors,
    compute_log_lik = FALSE,
    use_reduce_sum = FALSE,
    truncation_check_multiplier = 2,
    nonparametric = NULL) {
  np <- .build_np_stan_fields(nonparametric)
  if (!is.null(nonparametric)) {
    dist_id <- np$dist_id
    param_bounds <- list(lower = numeric(0), upper = numeric(0))
    priors <- list(location = numeric(0), scale = numeric(0))
  }
  required_cols <- c(delay, delay_upper, n, pwindow, relative_obs_time)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns: ",
      toString(missing_cols),
      "\n",
      "Please ensure your data frame contains these columns or set the",
      " corresponding arguments:\n",
      "delay = '",
      delay,
      "'\n",
      "delay_upper = '",
      delay_upper,
      "'\n",
      "n = '",
      n,
      "'\n",
      "pwindow = '",
      pwindow,
      "'\n",
      "relative_obs_time = '",
      relative_obs_time,
      "'",
      call. = FALSE
    )
  }

  # Handle start_relative_obs_time column: if not present, default to 0
  if (!start_relative_obs_time %in% names(data)) {
    data[[start_relative_obs_time]] <- 0
  }

  # Validate truncation bounds: L must be less than D
  .check_truncation_bounds_df(data, start_relative_obs_time, relative_obs_time)

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
    L = data[[start_relative_obs_time]],
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
    primary_prior_scale = primary_priors$scale,
    nonparametric = np$nonparametric,
    K_np = np$K_np,
    np_boundaries = np$np_boundaries,
    np_paramtype = np$np_paramtype,
    np_dirichlet_alpha = np$np_dirichlet_alpha,
    np_alpha_mean = np$np_alpha_mean,
    np_alpha_sd = np$np_alpha_sd,
    np_log_sigma_mean = np$np_log_sigma_mean,
    np_log_sigma_sd = np$np_log_sigma_sd
  )

  return(stan_data)
}

# Build the np_* Stan data fields. When `nonparametric` is NULL all fields
# are sized 0 and `np_paramtype` defaults to 1 (a no-op for the parametric
# path). When supplied, validates the list and returns Stan-shaped fields
# plus the inferred `dist_id` (26 for "simplex", 27 for "hazard").
.build_np_stan_fields <- function(nonparametric) {
  if (is.null(nonparametric)) {
    # When the non-parametric path is off all np_* fields collapse to the
    # smallest valid size for their declared types. `np_boundaries` is
    # `vector[K_np + 1]` so it cannot be empty even when `K_np = 0`.
    return(list(
      nonparametric = 0L,
      K_np = 0L,
      np_boundaries = 0,
      np_paramtype = 1L,
      np_dirichlet_alpha = numeric(0),
      np_alpha_mean = 0,
      np_alpha_sd = 1,
      np_log_sigma_mean = 0,
      np_log_sigma_sd = 1,
      dist_id = NA_integer_
    ))
  }
  required <- c("K", "boundaries", "paramtype")
  missing_fields <- setdiff(required, names(nonparametric))
  if (length(missing_fields) > 0) {
    stop(
      "`nonparametric` is missing required elements: ",
      toString(missing_fields),
      call. = FALSE
    )
  }
  K <- as.integer(nonparametric$K)
  if (!is.finite(K) || K < 1) {
    stop("`nonparametric$K` must be a positive integer.", call. = FALSE)
  }
  boundaries <- as.numeric(nonparametric$boundaries)
  if (length(boundaries) != K + 1L) {
    stop(
      "`nonparametric$boundaries` must have length K + 1 (got ",
      length(boundaries), ", expected ", K + 1L, ").",
      call. = FALSE
    )
  }
  paramtype <- match.arg(
    nonparametric$paramtype,
    choices = c("simplex", "hazard", "hazardrw", "hazardre")
  )
  np_paramtype <- switch(paramtype,
    simplex = 1L,
    hazard = 2L,
    hazardrw = 2L,
    hazardre = 3L
  )
  dist_id <- if (paramtype == "simplex") 26L else 27L
  dirichlet_alpha <- nonparametric$dirichlet_alpha %||% rep(1, K)
  if (length(dirichlet_alpha) != K) {
    stop(
      "`nonparametric$dirichlet_alpha` must have length K.",
      call. = FALSE
    )
  }
  hazard_priors <- nonparametric$hazard_priors %||% list()
  hp <- utils::modifyList(
    list(
      alpha_mean = 0,
      alpha_sd = 5,
      log_sigma_mean = 0,
      log_sigma_sd = 1
    ),
    hazard_priors
  )
  list(
    nonparametric = 1L,
    K_np = K,
    np_boundaries = boundaries,
    np_paramtype = np_paramtype,
    np_dirichlet_alpha = as.numeric(dirichlet_alpha),
    np_alpha_mean = as.numeric(hp$alpha_mean),
    np_alpha_sd = as.numeric(hp$alpha_sd),
    np_log_sigma_mean = as.numeric(hp$log_sigma_mean),
    np_log_sigma_sd = as.numeric(hp$log_sigma_sd),
    dist_id = dist_id
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a
