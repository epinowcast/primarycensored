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
#' Parametric distributions are fitted in the parameterisation named by
#' `start`, and the returned estimates and covariance matrix use the same
#' names. For example, gamma can be fitted with either
#' `start = list(shape = , rate = )` or `start = list(shape = , scale = )`.
#' Parameters can be held fixed by passing `fix.arg` to
#' [fitdistrplus::fitdist()] through `...`, either as a named list or as a
#' function of the delays returning one.
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
#'   parameters drive either a Gaussian random walk on the logit hazard
#'   (`hazard_model = "rw"`, the default,
#'   `logit(h_i) = alpha + sigma * cumsum(eps)`) or an IID logit
#'   random-effect transform (`hazard_model = "re"`,
#'   `logit(h_i) = alpha + sigma * eps_i` with `eps_i ~ N(0, 1)`). The
#'   smoothing of the random walk regularises the recovered PMF against
#'   over-fitting in sparse data and replaces the simplex constraint
#'   with an unconstrained optimisation; the random-effect variant
#'   models hazards as independent draws around `alpha` rather than a
#'   smoothed trajectory. See [pdiscretehazard()] for full
#'   parameterisation details and the MAP-equivalent prior penalty
#'   applied during fitting; pass `prior` to override the default prior
#'   settings.
#'
#' For non-parametric distributions `K` is implied by `length(start)`:
#' `K = length(start) + 1` for `"discretestep"` and
#' `K = length(start) - 1` for `"discretehazard"`. `start` is therefore
#' required.
#'
#' ## Exact observations
#'
#' Rows with `pwindow = 0` have an exactly known primary event time. Rows
#' with `left == right` have an exactly known secondary event time and
#' contribute a density rather than a probability to the likelihood (see
#' [dprimarycensored()]). Rows of different types can be mixed in one fit,
#' as for the exact, single interval censored and doubly interval censored
#' observations of `coarseDataTools::dic.fit()`. Data with the primary
#' event in \[`EL`, `ER`\] and the secondary event in \[`SL`, `SR`\] map
#' to `left = SL - EL`, `right = SR - EL` and `pwindow = ER - EL`.
#'
#' @param censdata A data frame with columns 'left' and 'right' representing
#'  the lower and upper bounds of the censored observations. Unlike
#'  [fitdistrplus::fitdistcens()] `NA` is not supported for either the
#'  upper or lower bounds. Use `left == right` for an exactly observed
#'  secondary event.
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
#' @param pwindow Column name for primary window (default: "pwindow"). Use
#'  a primary window of 0 for an exactly observed primary event.
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
#' @param hazard_model One of `"rw"` (default) or `"re"`. Only consulted
#'   when `distr = "discretehazard"`. `"rw"` selects the random-walk
#'   transform `logit(h_i) = alpha + sigma * cumsum(eps)`; `"re"`
#'   selects the IID logit random-effect transform
#'   `logit(h_i) = alpha + sigma * eps_i`. See Details.
#'
#' @param ... Additional arguments to be passed to [fitdistrplus::fitdist()].
#'
#' @param check Logical; if `TRUE` (the default) `pdist` is validated with
#'   [check_pdist()] and `dprimary` with [check_dprimary()]. Neither changes
#'   across a fit, so validation runs on the first likelihood evaluation
#'   only rather than on every one. Set to `FALSE` to skip it entirely.
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
    dprimary = dunif,
    primary_args = NULL,
    pprimary = NULL,
    dprimary_args = NULL,
    truncation_check_multiplier = 2,
    prior = NULL,
    hazard_model = c("rw", "re"),
    check = TRUE,
    ...) {
  hazard_model <- match.arg(hazard_model)
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
    bad_left <- censdata[[left]][invalid_upper[1]]
    bad_D <- censdata[[D]][invalid_upper[1]]
    stop(
      "Upper truncation point is greater than D. Maximum ", left, " is ",
      max(censdata[[left]][invalid_upper]),
      " and D is ", bad_D,
      ". Under truncation at D no event with latent value >= D is ",
      "observable; resolve this by filtering ", left,
      " to values strictly less than D. Found ", length(invalid_upper),
      " observation(s) where ", left, " >= D; first invalid row: ",
      invalid_upper[1], " (", left, " = ", bad_left,
      ", D = ", bad_D, ").",
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
  fit_bounds <- attr(ddist, "fit_bounds")

  # For the hazard family the transform is selected at fit time from
  # `hazard_model`; for other families it is whatever the dist function
  # carries as a `param_transform` attribute (currently only the simplex
  # transform on the discretestep family).
  if (identical(vector_param, "hazards")) {
    param_transform <- .make_hazard_transform(hazard_model)
  }

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

  # fitdistrplus accepts `fix.arg` as a list or as a function of the data
  # returning one. Either way the fixed names must be arguments of the
  # synthetic density.
  fix_arg <- dots$fix.arg
  if (is.function(fix_arg)) {
    fix_arg <- fix_arg(delays)
  }
  fix_names <- names(fix_arg)

  # `pdist` and `dprimary` are fixed across the fit, so validate on the first
  # likelihood evaluation and skip it thereafter. Revalidating on every
  # evaluation costs four extra `pdist` calls each time and advances the RNG
  # stream, which makes seeded runs depend on the number of evaluations.
  validation <- new.env(parent = emptyenv())
  validation$pending <- isTRUE(check)
  check_once <- function() {
    if (!validation$pending) {
      return(FALSE)
    }
    validation$pending <- FALSE
    TRUE
  }

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
    fix_names = fix_names,
    pdist_extras = pdist_extras,
    check_once = check_once
  )

  # If the dist function carries a `fit_bounds` attribute, use it to fill
  # in any `lower`/`upper` the caller did not supply.
  if (!is.null(fit_bounds)) {
    fb <- fit_bounds(closures$par_names)
    if (is.null(dots$lower)) dots$lower <- fb$lower
    if (is.null(dots$upper)) dots$upper <- fb$upper
  }

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

# Closure builder (`.build_pcens_closures`) lives in
# R/nonparametric_helpers.R alongside the rest of the non-parametric
# machinery. Default `lower`/`upper` bounds are now declared by the
# dist function itself via the `fit_bounds` attribute (see
# `pdiscretestep()` and `pdiscretehazard()`).

# ---- low-level wrappers ----------------------------------------------------

#' Define a fitdistrplus compatible wrapper around dprimarycensored
#' @inheritParams dprimarycensored
#'
#' @param params A data frame with columns 'swindow', 'pwindow', 'L', and 'D'
#' corresponding to the secondary window sizes, primary window sizes, upper
#' truncation times, and lower truncation times for each element in x.
#'
#' @param pcens_cache Optional environment shared across calls with the same
#'   `params`, `pdist` and `dprimary`, as made by [.build_pcens_closures()].
#'   The `pcens` object and the grouping of `params` are built on the first
#'   call and kept in it, and later calls only [update()][update.pcens()] the
#'   parameters. `NULL` (the default) builds them on every call.
#' @keywords internal
.dpcens <- function(
    x,
    params,
    pdist,
    dprimary,
    primary_args,
    pprimary = NULL,
    check = TRUE,
    pcens_cache = NULL,
    ...) {
  # Wrap in `suppressMessages` so the per-call upper-clip notice from
  # pcens_pmf() is not emitted on every fitdistrplus iteration.
  suppressMessages(tryCatch(
    {
      # Validate once for the whole vector. `pdist` and `dprimary` are the
      # same for every observation.
      if (isTRUE(check)) {
        check_pdist(pdist, D = max(params$D), ...)
        for (pw in unique(params$pwindow)) {
          check_dprimary(dprimary, pw, primary_args)
        }
      }

      state <- .fit_pcens_state(
        pcens_cache, pdist, dprimary, primary_args, pprimary, list(...)
      )
      if (is.null(state$dgroups)) {
        state$dgroups <- .param_groups(
          params, c("swindow", "pwindow", "L", "D")
        )
      }
      groups <- state$dgroups
      if (length(groups) == 1L) {
        g <- groups[[1L]]
        pcens_pmf(
          state$obj, x, g$pwindow,
          swindow = g$swindow, L = g$L, D = g$D
        )
      } else {
        result <- numeric(length(x))
        for (g in groups) {
          result[g$mask] <- pcens_pmf(
            state$obj, x[g$mask], g$pwindow,
            swindow = g$swindow, L = g$L, D = g$D
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
#' @inheritParams .dpcens
#' @keywords internal
.ppcens <- function(q, params, pdist, dprimary, primary_args, pprimary = NULL,
                    check = TRUE, pcens_cache = NULL, ...) {
  tryCatch(
    {
      # Validate once for the whole vector. `pdist` and `dprimary` are the
      # same for every observation.
      if (isTRUE(check)) {
        check_pdist(pdist, D = max(params$D), ...)
        for (pw in unique(params$pwindow)) {
          check_dprimary(dprimary, pw, primary_args)
        }
      }

      state <- .fit_pcens_state(
        pcens_cache, pdist, dprimary, primary_args, pprimary, list(...)
      )
      obj <- state$obj
      cdf <- function(q_i, pw, L_i, D_i) {
        .check_truncation_bounds(L_i, D_i)
        # Evaluate the CDF first, as .normalise_cdf() may not use it
        result <- pcens_cdf(obj, q_i, pw)
        .normalise_cdf(result, q_i, L_i, D_i, obj, pw)
      }

      if (length(q) != nrow(params)) {
        # Recycle as mapply() does
        return(mapply(
          cdf, q, params$pwindow, params$L, params$D,
          SIMPLIFY = TRUE
        ))
      }
      if (is.null(state$pgroups)) {
        state$pgroups <- .param_groups(params, c("pwindow", "L", "D"))
      }
      result <- numeric(length(q))
      for (g in state$pgroups) {
        result[g$mask] <- cdf(q[g$mask], g$pwindow, g$L, g$D)
      }
      names(result) <- names(q)
      result
    },
    error = function(e) {
      rep(NaN, length(q))
    }
  )
}

#' Get the pcens object for a likelihood evaluation
#'
#' Builds a `pcens` object with [.build_pcens()], or, when `cache` already
#' holds one, updates its delay parameters with [update()][update.pcens()].
#'
#' @inheritParams .dpcens
#'
#' @param cache Environment to keep the object in, or `NULL`.
#'
#' @param args Named list of delay distribution parameters.
#'
#' @return An environment with the `pcens` object in `obj`. This is `cache`
#'   when it is not `NULL`.
#'
#' @keywords internal
.fit_pcens_state <- function(cache, pdist, dprimary, primary_args, pprimary,
                             args) {
  if (!is.null(cache) && !is.null(cache$obj)) {
    cache$obj <- do.call(update, c(list(cache$obj), args))
    return(cache)
  }
  state <- cache
  if (is.null(state)) {
    state <- new.env(parent = emptyenv())
  }
  if (is.null(primary_args)) {
    primary_args <- list()
  }
  state$obj <- .build_pcens(
    pdist, dprimary, primary_args, pprimary, args,
    pwindow = NULL, D = NULL, check = FALSE
  )
  state
}

#' Group observations that share censoring and truncation settings
#'
#' @param params A data frame of per-observation settings.
#'
#' @param cols Names of the columns to group by.
#'
#' @return A list with one element per unique combination of `cols`. Each
#'   element is a list of the values of `cols` and a logical `mask` selecting
#'   the rows of `params` with those values.
#'
#' @keywords internal
.param_groups <- function(params, cols) {
  keys <- unique(params[cols])
  lapply(seq_len(nrow(keys)), function(i) {
    group <- lapply(keys[cols], `[[`, i)
    mask <- rep(TRUE, nrow(params))
    for (col in cols) {
      mask <- mask & params[[col]] == group[[col]]
    }
    group$mask <- mask
    group
  })
}
