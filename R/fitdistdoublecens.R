#' Fit a distribution to doubly censored data
#'
#' This function wraps the custom approach for fitting distributions to doubly
#' censored data using fitdistrplus and primarycensored. It handles primary
#' censoring (when the primary event time is not known exactly), secondary
#' censoring (when the secondary event time is interval-censored), and
#' truncation (when events are only observed within a delay range \[L, D\]).
#'
#' @details
#' ## How distribution functions are resolved
#'
#' The `distr` argument names a distribution. The function looks up the
#' density and CDF functions by prepending `d` and `p` to the name (e.g.
#' `distr = "gamma"` resolves to `dgamma()` and `pgamma()`). Custom
#' distributions can be used as long as the corresponding `d<distr>()` and
#' `p<distr>()` functions are defined.
#'
#' ## Non-parametric distributions
#'
#' Two non-parametric distributions are supported. They share a common
#' fitting machinery: the dist function carries a `vector_param` attribute
#' (`"pmf"` for [pdiscretestep()]/[ddiscretestep()], `"hazards"` for
#' [pdiscretehazard()]/[ddiscretehazard()]) that drives this function to
#' build a closure mapping flat scalar parameters into the underlying
#' vector argument.
#'
#' - `distr = "discretestep"`: free parameters `p1, ..., p_{K-1}` (in
#'   `[0, 1]`); the last bin probability is `1 - sum(p1, ..., p_{K-1})`.
#'   See [pdiscretestep()] for parameterisation details and the soft
#'   simplex penalty applied when probabilities are infeasible.
#' - `distr = "discretehazard"`: free parameters `alpha`, `log_sigma`,
#'   `eps_1, ..., eps_{K-1}`. The hazard form parameterises the same
#'   family of step distributions as `"discretestep"`, but its free
#'   parameters drive a Gaussian random walk on the logit hazard. The
#'   smoothing regularises the recovered PMF against over-fitting in
#'   sparse data and replaces the simplex constraint with an
#'   unconstrained optimisation, which tends to be more stable. See
#'   [pdiscretehazard()] for full parameterisation details and the
#'   MAP-equivalent prior penalty applied during fitting; pass `prior`
#'   to override the default prior settings.
#'
#' For non-parametric distributions `K` is implied by `length(start)`:
#' `K = length(start) + 1` for `"discretestep"` and
#' `K = length(start) - 1` for `"discretehazard"`. `start` is therefore
#' required.
#'
#' @param censdata A data frame with columns 'left' and 'right' representing
#'  the lower and upper bounds of the censored observations. Unlike
#'  [fitdistrplus::fitdistcens()] `NA` is not supported for either the
#'  upper or lower bounds.
#'
#' @param distr A character string naming the distribution to be fitted.
#'  Special values `"discretestep"` and `"discretehazard"` select the
#'  non-parametric step-distribution fitting; see Details.
#'
#' @param left Column name for lower bound of observed values (default:
#'  "left").
#'
#' @param right Column name for upper bound of observed values (default:
#'  "right").
#'
#' @param pwindow Column name for primary window (default: "pwindow").
#'
#' @param L Column name for minimum delay (lower truncation point). For any
#'  finite L the distribution is left-truncated at L; use `L = -Inf` for no
#'  left truncation. This is useful for modelling generation intervals where
#'  day 0 is excluded, particularly when used in renewal models. (default:
#'  "L"). If the column is not present in censdata, L = -Inf is assumed.
#'
#' @param D Column name for maximum delay (upper truncation point). If finite,
#'  the distribution is truncated at D. If set to Inf, no upper truncation is
#'  applied. (default: "D"). Observations whose secondary censoring interval
#'  straddles `D` (`left < D <= right`) are accepted: the upper endpoint is
#'  internally clipped to `D` and the likelihood becomes
#'  `P(X in [left, min(right, D)] | L <= X <= D)`. This is a no-op for the
#'  standard parametric case where `right <= D`. Observations with
#'  `left >= D` are rejected because under truncation at `D` no event with
#'  latent value `>= D` is observable.
#'
#' @inheritParams pprimarycensored
#'
#' @param prior Optional list of prior settings used by the dist function's
#'   `fit_penalty` attribute (currently only `"discretehazard"`). Each
#'   element is itself a list with `mean` and `sd` entries. Defaults are
#'   used for any component not supplied. See [pdiscretehazard()] for
#'   the default values.
#'
#' @param ... Additional arguments to be passed to [fitdistrplus::fitdist()].
#'   For non-parametric distributions, `start` is required and determines
#'   the number of bins; pass `boundaries` here to override the
#'   default `0:K` unit-width bins.
#'
#' @param truncation_check_multiplier Numeric multiplier to use for checking
#'   if the truncation time D is appropriate relative to the maximum delay.
#'   Set to NULL to skip the check. Default is 2.
#'
#' @return An object of class "fitdist" as returned by fitdistrplus::fitdist.
#'
#' @export
#' @family modelhelpers
#' @seealso [pdiscretestep()] [pdiscretehazard()]
#' @examplesIf requireNamespace("fitdistrplus", quietly = TRUE)
#' # Example with normal distribution
#' set.seed(123)
#' n <- 1000
#' true_mean <- 5
#' true_sd <- 2
#' pwindow <- 2
#' swindow <- 2
#' D <- 10
#' samples <- rprimarycensored(
#'   n, rnorm,
#'   mean = true_mean, sd = true_sd,
#'   pwindow = pwindow, swindow = swindow, D = D
#' )
#'
#' delay_data <- data.frame(
#'   left = samples,
#'   right = samples + swindow,
#'   pwindow = rep(pwindow, n),
#'   D = rep(D, n)
#' )
#'
#' fit_norm <- fitdistdoublecens(
#'   delay_data,
#'   distr = "norm",
#'   start = list(mean = 0, sd = 1)
#' )
#'
#' summary(fit_norm)
#'
#' \donttest{
#' # Example with discretestep (non-parametric PMF) distribution
#' set.seed(42)
#' true_pmf <- c(0.1, 0.3, 0.4, 0.15, 0.05)
#' step_samples <- rprimarycensored(
#'   500, rdiscretestep,
#'   boundaries = 0:5, pmf = true_pmf,
#'   pwindow = 1, swindow = 1, D = 6
#' )
#' step_data <- data.frame(
#'   left = step_samples,
#'   right = step_samples + 1,
#'   pwindow = rep(1, 500),
#'   D = rep(6, 500)
#' )
#' fit_step <- fitdistdoublecens(
#'   step_data,
#'   distr = "discretestep",
#'   boundaries = 0:5,
#'   start = as.list(setNames(rep(0.2, 4), paste0("p", 1:4)))
#' )
#'
#' # Example with discretehazard (logit-hazard random walk) distribution
#' fit_haz <- fitdistdoublecens(
#'   step_data,
#'   distr = "discretehazard",
#'   boundaries = 0:5,
#'   start = c(
#'     list(alpha = -2, log_sigma = log(0.5)),
#'     as.list(setNames(rep(0, 4), paste0("eps_", 1:4)))
#'   )
#' )
#' }
fitdistdoublecens <- function(
    censdata,
    distr,
    left = "left",
    right = "right",
    pwindow = "pwindow",
    L = "L",
    D = "D",
    dprimary = stats::dunif,
    primary_args = NULL,
    pprimary = NULL,
    dprimary_args = NULL,
    truncation_check_multiplier = 2,
    prior = NULL,
    ...) {
  if (!requireNamespace("fitdistrplus", quietly = TRUE)) {
    stop(
      "Package 'fitdistrplus' is required but not installed for this function.",
      call. = FALSE
    )
  }
  if (!requireNamespace("withr", quietly = TRUE)) {
    stop(
      "Package 'withr' is required but not installed for this function.",
      call. = FALSE
    )
  }

  primary_args <- .resolve_primary_args(
    primary_args, dprimary_args, "fitdistdoublecens"
  )

  # Handle L column: if not present, default to -Inf (no left truncation).
  if (!L %in% names(censdata)) {
    censdata[[L]] <- -Inf
  }

  .check_truncation_bounds_df(censdata, L, D)
  invalid_obs <- which(censdata[[left]] < censdata[[L]])
  if (length(invalid_obs) > 0) {
    stop(
      "Observations must be >= L. Found ", length(invalid_obs),
      " observation(s) where ", left, " < L. First invalid row: ",
      invalid_obs[1], " (", left, " = ", censdata[[left]][invalid_obs[1]],
      ", L = ", censdata[[L]][invalid_obs[1]], ")",
      call. = FALSE
    )
  }
  invalid_upper <- which(
    is.finite(censdata[[D]]) & censdata[[left]] >= censdata[[D]]
  )
  if (length(invalid_upper) > 0) {
    stop(
      "Observations must be < D. Found ", length(invalid_upper),
      " observation(s) where ", left, " >= D. First invalid row: ",
      invalid_upper[1], " (", left, " = ",
      censdata[[left]][invalid_upper[1]],
      ", D = ", censdata[[D]][invalid_upper[1]], ")",
      ". Under truncation at D no event with latent value >= D is ",
      "observable.",
      call. = FALSE
    )
  }
  required_cols <- c(left, right, pwindow, D)
  missing_cols <- setdiff(required_cols, names(censdata))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns: ",
      toString(missing_cols),
      call. = FALSE
    )
  }
  if (!is.null(truncation_check_multiplier)) {
    unique_D <- unique(censdata[[D]])
    for (d in unique_D) {
      delays_subset <- censdata[[left]][censdata[[D]] == d]
      check_truncation(
        delays = delays_subset,
        D = d,
        multiplier = truncation_check_multiplier
      )
    }
  }

  pdist_name <- paste0("p", distr)
  ddist_name <- paste0("d", distr)
  pdist <- add_name_attribute(get(pdist_name), pdist_name)
  ddist <- get(ddist_name)

  params <- data.frame(
    swindow = censdata[[right]] - censdata[[left]],
    pwindow = censdata[[pwindow]],
    L = censdata[[L]],
    D = censdata[[D]]
  )
  delays <- censdata[[left]]
  N <- length(delays)

  vector_param <- attr(ddist, "vector_param")
  fit_penalty <- attr(ddist, "fit_penalty")
  param_transform <- attr(ddist, "param_transform")

  dots <- list(...)

  # Separate any extra distribution-level args (e.g. `boundaries` for the
  # step distribution) from arguments destined for fitdistrplus::fitdist.
  fitdist_arg_names <- c(
    "start", "fix.arg", "lower", "upper", "method", "optim.method",
    "custom.optim", "discrete", "weights", "silent", "calcvcov",
    "checkstartfix", "keepdata", "keepdata.nb", "control"
  )
  pdist_extras <- dots[setdiff(names(dots), fitdist_arg_names)]
  dots <- dots[intersect(names(dots), fitdist_arg_names)]

  closures <- .build_pcens_closures(
    pdist = pdist,
    ddist = ddist,
    params = params,
    dprimary = dprimary,
    primary_args = primary_args,
    pprimary = pprimary,
    vector_param = vector_param,
    param_transform = param_transform,
    fit_penalty = fit_penalty,
    prior = prior,
    N = N,
    start = dots$start,
    pdist_extras = pdist_extras
  )

  # Apply default bounds for non-parametric fits when not supplied.
  dots <- .nonparametric_defaults(
    dots = dots, vector_param = vector_param, par_names = closures$par_names
  )

  fit_env <- new.env(parent = emptyenv())
  fit_env$delays <- delays
  fit_env$dpcens_dist <- closures$dpcens_dist
  fit_env$ppcens_dist <- closures$ppcens_dist

  fit_args <- c(
    list(delays, distr = "pcens_dist"),
    dots
  )

  withr::with_environment(
    fit_env,
    do.call(fitdistrplus::fitdist, fit_args)
  )
}

# ---- closure builder -------------------------------------------------------

#' Build the (dpcens_dist, ppcens_dist) closure pair for fitdistdoublecens
#'
#' Single path: scalar named parameters are gathered from the closure's
#' environment, optionally folded into the vector argument named by
#' \code{vector_param}, and dispatched through \code{.dpcens}/\code{.ppcens}.
#' Formals on the closures are derived from the supplied \code{start} list
#' (parametric) or from a vector_param-aware naming convention
#' (non-parametric).
#'
#' @keywords internal
.build_pcens_closures <- function(
    pdist,
    ddist,
    params,
    dprimary,
    primary_args,
    pprimary = NULL,
    vector_param,
    param_transform = NULL,
    fit_penalty,
    prior,
    N,
    start,
    pdist_extras = list()) {
  if (is.null(vector_param)) {
    # Parametric path: parameter names come from start, or fall back to
    # the formals of the d<distr> function.
    if (!is.null(start)) {
      par_names <- names(start)
    } else {
      par_names <- setdiff(names(formals(ddist)), c("x", "log"))
    }
  } else {
    if (is.null(start)) {
      stop(
        "Argument `start` must be supplied for distr with ",
        "vector_param = '", vector_param, "'.",
        call. = FALSE
      )
    }
    par_names <- names(start)
  }

  # The closure body looks up scalar named args from its environment,
  # builds either a numeric vector (non-parametric) or a flat list
  # (parametric), and dispatches through .dpcens / .ppcens.
  build_call_args <- function(env_args) {
    if (is.null(vector_param)) {
      # Parametric: pass through scalar args by name.
      return(env_args[par_names])
    }
    par_named <- env_args[par_names]
    if (!is.null(param_transform)) {
      # Distribution-supplied mapping from named scalar args to a numeric
      # vector. Used by the logit-hazard random walk, where the free
      # parameters (alpha, log_sigma, eps_*) are not the hazards themselves.
      full_vec <- param_transform(par_named)
    } else if (vector_param == "pmf") {
      free <- unlist(par_named, use.names = FALSE)
      full_vec <- c(free, 1 - sum(free))
    } else {
      full_vec <- unlist(par_named, use.names = FALSE)
    }
    stats::setNames(list(full_vec), vector_param)
  }

  dpcens_dist <- function() {
    env_args <- as.list(environment())
    extra <- build_call_args(env_args)
    base_dens <- do.call(
      .dpcens,
      c(
        list(
          x = env_args$x,
          params = params,
          pdist = pdist,
          dprimary = dprimary,
          primary_args = primary_args,
          pprimary = pprimary
        ),
        pdist_extras,
        extra
      )
    )
    if (!is.null(fit_penalty)) {
      par_named <- env_args[par_names]
      penalty_per_obs <- fit_penalty(par_named, N = N, prior_settings = prior)
      base_dens <- pmax(
        base_dens * exp(-penalty_per_obs),
        .Machine$double.eps
      )
    }
    base_dens
  }

  ppcens_dist <- function() {
    env_args <- as.list(environment())
    extra <- build_call_args(env_args)
    do.call(
      .ppcens,
      c(
        list(
          q = env_args$q,
          params = params,
          pdist = pdist,
          dprimary = dprimary,
          primary_args = primary_args,
          pprimary = pprimary
        ),
        pdist_extras,
        extra
      )
    )
  }

  free_formals <- vector("list", length(par_names))
  names(free_formals) <- par_names
  formals(dpcens_dist) <- c(alist(x = ), free_formals)
  formals(ppcens_dist) <- c(alist(q = ), free_formals)

  list(
    dpcens_dist = dpcens_dist,
    ppcens_dist = ppcens_dist,
    par_names = par_names
  )
}

#' Apply default bounds for non-parametric fits
#'
#' Sets sensible defaults for \code{lower} and \code{upper} when the user
#' has not supplied them, based on the \code{vector_param} kind.
#'
#' @keywords internal
.nonparametric_defaults <- function(dots, vector_param, par_names) {
  if (is.null(vector_param)) {
    return(dots)
  }
  if (vector_param == "pmf") {
    if (is.null(dots$lower)) dots$lower <- rep(0, length(par_names))
    if (is.null(dots$upper)) dots$upper <- rep(1, length(par_names))
    return(dots)
  }
  if (vector_param == "hazards") {
    # Expect par_names = c("alpha", "log_sigma", eps_1, ..., eps_{K-1}).
    n_eps <- length(par_names) - 2L
    if (is.null(dots$lower)) {
      dots$lower <- c(-Inf, -6, rep(-5, n_eps))
    }
    if (is.null(dots$upper)) {
      dots$upper <- c(Inf, 4, rep(5, n_eps))
    }
  }
  dots
}

# ---- low-level wrappers ----------------------------------------------------

#' Define a fitdistrplus compatible wrapper around dprimarycensored
#' @inheritParams dprimarycensored
#'
#' @param params A data frame with columns 'swindow', 'pwindow', 'L', and 'D'
#' corresponding to the secondary window sizes, primary window sizes, upper
#' truncation times, and lower truncation times for each element in x.
#' @keywords internal
.dpcens <- function(
    x,
    params,
    pdist,
    dprimary,
    primary_args,
    pprimary = NULL,
    ...) {
  # Wrap in `suppressMessages` so the per-call upper-clip notice from
  # dprimarycensored() is not emitted on every fitdistrplus iteration.
  suppressMessages(tryCatch(
    {
      unique_params <- unique(params)
      if (nrow(unique_params) == 1) {
        dprimarycensored(
          x,
          pdist,
          pwindow = unique_params$pwindow[1],
          swindow = unique_params$swindow[1],
          L = unique_params$L[1],
          D = unique_params$D[1],
          dprimary = dprimary,
          primary_args = primary_args,
          pprimary = pprimary,
          ...
        )
      } else {
        result <- numeric(length(x))
        for (i in seq_len(nrow(unique_params))) {
          sw <- unique_params$swindow[i]
          pw <- unique_params$pwindow[i]
          Ds <- unique_params$D[i] # nolint
          Ls <- unique_params$L[i] # nolint
          mask <- params$swindow == sw &
            params$pwindow == pw &
            params$D == Ds &
            params$L == Ls
          result[mask] <- dprimarycensored(
            x[mask],
            pdist,
            pwindow = pw,
            swindow = sw,
            L = Ls,
            D = Ds,
            dprimary = dprimary,
            primary_args = primary_args,
            pprimary = pprimary,
            ...
          )
        }
        result
      }
    },
    error = function(e) {
      rep(NaN, length(x))
    }
  ))
}

#' Define a fitdistrplus compatible wrapper around pprimarycensored
#' @inheritParams pprimarycensored
#' @keywords internal
.ppcens <- function(q, params, pdist, dprimary, primary_args, pprimary = NULL,
                    ...) {
  tryCatch(
    {
      mapply(
        function(q_i, pw, L_i, D_i) {
          pprimarycensored(
            q_i,
            pdist,
            pwindow = pw,
            L = L_i,
            D = D_i,
            dprimary = dprimary,
            primary_args = primary_args,
            pprimary = pprimary,
            ...
          )
        },
        q,
        params$pwindow,
        params$L,
        params$D,
        SIMPLIFY = TRUE
      )
    },
    error = function(e) {
      rep(NaN, length(q))
    }
  )
}
