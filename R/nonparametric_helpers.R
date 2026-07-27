# Internal helpers for the non-parametric delay families
# (`discretestep` and `discretehazard`) and their integration with
# `fitdistdoublecens()`. None of these are exported.

#' Test whether a numeric vector lies on the simplex
#'
#' Returns \code{TRUE} when the vector contains no missing values, no
#' negative entries, and sums to 1 within \eqn{10^{-8}}. Used by the
#' \code{discretestep} family to apply a soft penalty inside fitting
#' closures rather than erroring.
#'
#' @keywords internal
.is_valid_simplex <- function(p, tol = 1e-8) {
  if (!is.numeric(p) || anyNA(p)) {
    return(FALSE)
  }
  if (any(p < 0)) {
    return(FALSE)
  }
  abs(sum(p) - 1) <= tol
}

#' Validate boundaries and pmf for step distribution functions
#' @keywords internal
.validate_discretestep_args <- function(boundaries, pmf) {
  if (!is.numeric(boundaries)) {
    stop("boundaries must be numeric.", call. = FALSE)
  }
  if (!is.numeric(pmf)) {
    stop("pmf must be numeric.", call. = FALSE)
  }
  if (anyNA(boundaries) || anyNA(pmf)) {
    stop(
      "boundaries and pmf must not contain NA values.",
      call. = FALSE
    )
  }
  if (length(boundaries) != length(pmf) + 1L) {
    stop(
      "length(boundaries) must equal length(pmf) + 1.",
      call. = FALSE
    )
  }
  if (any(diff(boundaries) <= 0)) {
    stop("boundaries must be strictly increasing.", call. = FALSE)
  }
  invisible(NULL)
}

#' Sort eps_* parameter names by trailing numeric suffix
#'
#' Lexicographic sorting of \code{c("eps_1", ..., "eps_10")} returns
#' \code{eps_10} before \code{eps_2}, which would scramble the random-walk
#' innovations when \eqn{K > 10}. This helper extracts the trailing integer
#' and orders by it.
#'
#' @keywords internal
.sort_eps_names <- function(eps_names) {
  if (length(eps_names) <= 1L) {
    return(eps_names)
  }
  idx <- suppressWarnings(as.integer(sub("^eps_", "", eps_names)))
  if (anyNA(idx)) {
    return(eps_names)
  }
  eps_names[order(idx)]
}

#' Build a hazard-vector transform for a given model
#'
#' Returns a function that maps a named list of scalar parameters
#' (\code{alpha}, \code{log_sigma}, \code{eps_1}, ..., \code{eps_{K-1}})
#' to a length-\eqn{K} hazard vector with the final entry pinned to 1.
#' `model = "rw"` uses a Gaussian random walk on the logit hazard;
#' `model = "re"` treats the innovations as IID logit random effects
#' around the intercept.
#'
#' @keywords internal
.make_hazard_transform <- function(model = c("rw", "re")) {
  model <- match.arg(model)
  function(par_named) {
    alpha_val <- par_named[["alpha"]]
    log_sigma_val <- par_named[["log_sigma"]]
    eps_names <- .sort_eps_names(
      grep("^eps_", names(par_named), value = TRUE)
    )
    eps_vals <- if (length(eps_names) > 0) {
      unlist(par_named[eps_names])
    } else {
      numeric(0)
    }
    sigma <- exp(log_sigma_val)
    K <- length(eps_vals) + 1L
    deltas <- if (model == "rw") {
      cumsum(c(0, eps_vals))
    } else {
      c(0, eps_vals)
    }
    logit_h <- alpha_val + sigma * deltas
    h <- 1 / (1 + exp(-logit_h))
    h[K] <- 1
    h
  }
}

#' Shared MAP penalty for the logit-hazard families
#'
#' Both RW and RE variants share priors on \code{alpha}, \code{log_sigma},
#' and \code{eps_*}; only how the hazards are built from those scalars
#' differs (see \code{.make_hazard_transform()}).
#'
#' @keywords internal
.hazard_fit_penalty <- function(par_named, N, prior_settings = NULL) {
  alpha_val <- par_named[["alpha"]]
  log_sigma_val <- par_named[["log_sigma"]]
  eps_names <- .sort_eps_names(
    grep("^eps_", names(par_named), value = TRUE)
  )
  eps_vals <- if (length(eps_names) > 0) {
    unlist(par_named[eps_names])
  } else {
    numeric(0)
  }
  prior <- .resolve_hazard_prior(prior_settings)
  # Negative log prior (so adding to -loglik amounts to negative log
  # posterior). Amortised across N observations.
  nlp <- 0.5 * sum(eps_vals^2) +
    (alpha_val - prior$alpha[["mean"]])^2 /
      (2 * prior$alpha[["sd"]]^2) +
    (log_sigma_val - prior$log_sigma[["mean"]])^2 /
      (2 * prior$log_sigma[["sd"]]^2)
  nlp / N
}

#' Resolve the prior settings for the hazard fit_penalty
#'
#' Defaults: \code{alpha ~ N(0, 5)}, \code{log_sigma ~ N(log(0.5), 1)}.
#' User overrides are merged in: each component of \code{prior_settings}
#' may itself be a list with \code{mean} and \code{sd} entries.
#'
#' @keywords internal
.resolve_hazard_prior <- function(prior_settings = NULL) {
  defaults <- list(
    alpha = c(mean = 0, sd = 5),
    log_sigma = c(mean = log(0.5), sd = 1)
  )
  if (is.null(prior_settings)) {
    return(defaults)
  }
  for (nm in names(prior_settings)) {
    if (!nm %in% names(defaults)) next
    user <- prior_settings[[nm]]
    if (!is.null(user[["mean"]])) {
      defaults[[nm]][["mean"]] <- user[["mean"]]
    }
    if (!is.null(user[["sd"]])) {
      defaults[[nm]][["sd"]] <- user[["sd"]]
    }
  }
  defaults
}

#' Build the (dpcens_dist, ppcens_dist) closure pair for fitdistdoublecens
#'
#' Single path: scalar named parameters are gathered from the closure's
#' environment, optionally folded into the vector argument named by
#' \code{vector_param} via the distribution's \code{param_transform}, and
#' dispatched through \code{.dpcens}/\code{.ppcens}. Formals on the closures
#' are derived from the supplied \code{start} list (parametric) or from a
#' \code{vector_param}-aware naming convention (non-parametric).
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

  # The closure body looks up scalar named args from its environment and
  # dispatches through .dpcens/.ppcens. For non-parametric distributions
  # the named scalars are first folded into a single numeric vector via
  # the distribution's `param_transform`; the result is passed as the
  # `vector_param` argument.
  build_call_args <- function(env_args) {
    if (is.null(vector_param)) {
      # Parametric: pass through scalar args by name.
      return(env_args[par_names])
    }
    par_named <- env_args[par_names]
    full_vec <- param_transform(par_named)
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
