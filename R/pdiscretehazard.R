#' Hazard-parameterised piecewise-constant CDF
#'
#' Returns the CDF of a discrete distribution specified by its bin-wise
#' discrete-time hazards. Converts \code{hazards} to a PMF via
#' [hazards_to_pmf()] then delegates to [pdiscretestep()]: the two
#' parameterisations describe the same family of step distributions and
#' yield identical CDFs.
#'
#' ## When to use the hazard parameterisation
#'
#' Outside fitting, [pdiscretehazard()] is a deterministic wrapper around
#' [pdiscretestep()]. For evaluating a CDF, sampling, or computing a PMF
#' there is no advantage over the direct-PMF parameterisation.
#'
#' The hazard parameterisation exists to support fitting via
#' [fitdistdoublecens()]. There the free parameters are not the hazards
#' directly but a Gaussian random walk on the logit hazard. This smooths
#' the recovered PMF, regularises against over-fitting in sparse data,
#' and turns a constrained simplex problem into an unconstrained
#' optimisation, all of which make the fit more stable. The PMF is
#' guaranteed to be valid by construction.
#'
#' ## Use with \code{fitdistdoublecens()}
#'
#' This function carries the attribute \code{vector_param = "hazards"} so
#' that [fitdistdoublecens()] can drive it from the logit-hazard random
#' walk described above. The free parameters are \code{alpha} (level),
#' \code{log_sigma} (log random-walk scale), and
#' \code{eps_1, ..., eps_{K-1}} (standardised innovations). The
#' logit-hazard for bin \eqn{i} is
#' \code{alpha + sigma * cumsum(c(0, eps_1, ..., eps_{K-1}))[i]}, where
#' \code{sigma = exp(log_sigma)}. The final hazard is forced to 1 so the
#' PMF sums to 1.
#'
#' A MAP-style prior penalty is amortised across observations so that
#' \code{fitdistrplus::fitdist} (which minimises \code{-sum(log d_i)})
#' implicitly optimises the posterior mode. The penalty is supplied via the
#' \code{fit_penalty} attribute on this function. The default corresponds
#' to independent priors:
#' \itemize{
#'   \item \code{eps_j ~ N(0, 1)} (random-walk innovations);
#'   \item \code{alpha ~ N(0, 5)} (weakly informative level prior);
#'   \item \code{log_sigma ~ N(log(0.5), 1)} (weakly informative scale).
#' }
#' Pass \code{prior} to [fitdistdoublecens()] to override the default
#' settings.
#'
#' @param q Numeric vector of quantiles.
#' @param boundaries Numeric vector of length \eqn{K+1} defining the bin
#'   edges. Must be strictly increasing. Defaults to \code{0:K}.
#' @param hazards Numeric vector of length \eqn{K} (or \eqn{K-1}, in which
#'   case a trailing 1 is appended) giving the discrete-time conditional
#'   hazard for each bin. Values must be in \eqn{[0, 1]}.
#'
#' @return Numeric vector of CDF values, the same length as \code{q}.
#'
#' @family pdiscretehazard
#' @seealso [pdiscretestep()], [hazards_to_pmf()], [fitdistdoublecens()]
#' @export
#' @examples
#' # Equivalent to pdiscretestep with the corresponding PMF
#' hazards <- c(0.3, 0.5, 1)
#' boundaries <- 0:3
#' pdiscretehazard(c(0.5, 1, 2, 3),
#'   boundaries = boundaries,
#'   hazards = hazards
#' )
pdiscretehazard <- function(q, boundaries = NULL, hazards) {
  if (is.null(boundaries)) {
    boundaries <- seq.int(0L, length(hazards))
  }
  pmf <- hazards_to_pmf(hazards) # nolint: object_usage_linter
  pdiscretestep(q, boundaries, pmf) # nolint: object_usage_linter
}
attr(pdiscretehazard, "name") <- "pdiscretehazard"
attr(pdiscretehazard, "vector_param") <- "hazards"
# Transform a named list of scalar parameters into the vector argument
# named by `vector_param`. For the logit-hazard random walk the inputs
# are alpha, log_sigma, eps_1, ..., eps_{K-1}; the output is a length-K
# hazards vector with the final entry forced to 1.
attr(pdiscretehazard, "param_transform") <- function(par_named) {
  alpha_val <- par_named[["alpha"]]
  log_sigma_val <- par_named[["log_sigma"]]
  eps_names <- grep("^eps_", names(par_named), value = TRUE)
  eps_vals <- if (length(eps_names) > 0) {
    unlist(par_named[eps_names])
  } else {
    numeric(0)
  }
  sigma <- exp(log_sigma_val)
  K <- length(eps_vals) + 1L
  logit_h <- alpha_val + sigma * cumsum(c(0, eps_vals))
  h <- 1 / (1 + exp(-logit_h))
  h[K] <- 1
  h
}
attr(pdiscretehazard, "fit_penalty") <- function(par_named, N,
                                                 prior_settings = NULL) {
  alpha_val <- par_named[["alpha"]]
  log_sigma_val <- par_named[["log_sigma"]]
  eps_names <- grep("^eps_", names(par_named), value = TRUE)
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

#' Hazard-parameterised piecewise-constant PMF
#'
#' Returns the probability mass for each value in \code{x}. Converts
#' \code{hazards} to a PMF via [hazards_to_pmf()] then delegates to
#' [ddiscretestep()].
#'
#' @inheritParams pdiscretehazard
#' @param x Numeric vector of values at which to evaluate the PMF.
#'
#' @return Numeric vector of PMF values, the same length as \code{x}.
#'
#' @family pdiscretehazard
#' @seealso [ddiscretestep()], [hazards_to_pmf()]
#' @export
#' @examples
#' hazards <- c(0.3, 0.5, 1)
#' ddiscretehazard(c(1, 2, 3), boundaries = 0:3, hazards = hazards)
ddiscretehazard <- function(x, boundaries = NULL, hazards) {
  if (is.null(boundaries)) {
    boundaries <- seq.int(0L, length(hazards))
  }
  pmf <- hazards_to_pmf(hazards) # nolint: object_usage_linter
  ddiscretestep(x, boundaries, pmf) # nolint: object_usage_linter
}
attr(ddiscretehazard, "name") <- "ddiscretehazard"
attr(ddiscretehazard, "vector_param") <- "hazards"
attr(ddiscretehazard, "fit_penalty") <- attr(pdiscretehazard, "fit_penalty")
attr(ddiscretehazard, "param_transform") <- attr(
  pdiscretehazard, "param_transform"
)

#' Sample from a hazard-parameterised step distribution
#'
#' Draws \code{n} independent samples from the discrete distribution defined
#' by \code{hazards}. Converts \code{hazards} to a PMF via
#' [hazards_to_pmf()] then delegates to [rdiscretestep()].
#'
#' @param n Integer. Number of samples to draw.
#' @inheritParams pdiscretehazard
#'
#' @return Numeric vector of length \code{n}.
#'
#' @family pdiscretehazard
#' @seealso [rdiscretestep()], [hazards_to_pmf()]
#' @export
#' @examples
#' set.seed(42)
#' rdiscretehazard(10, boundaries = 0:3, hazards = c(0.3, 0.5, 1))
rdiscretehazard <- function(n, boundaries = NULL, hazards) {
  if (is.null(boundaries)) {
    boundaries <- seq.int(0L, length(hazards))
  }
  pmf <- hazards_to_pmf(hazards) # nolint: object_usage_linter
  rdiscretestep(n, boundaries, pmf) # nolint: object_usage_linter
}
attr(rdiscretehazard, "name") <- "rdiscretehazard"
attr(rdiscretehazard, "vector_param") <- "hazards"

# ---- internal helpers -------------------------------------------------------

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
