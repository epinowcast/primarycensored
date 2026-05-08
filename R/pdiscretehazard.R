#' Hazard-parameterised piecewise-constant CDF
#'
#' Returns the CDF of a discrete distribution specified by its bin-wise
#' discrete-time hazards. Converts \code{hazards} to a PMF via
#' [hazards_to_pmf()] then delegates to [pdiscretestep()].
#'
#' Outside fitting it is a deterministic wrapper around [pdiscretestep()].
#' It earns its keep as a fitting parameterisation in [fitdistdoublecens()]
#' because the random walk on the logit hazard smooths the recovered PMF.
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
#' @concept pdiscretehazard
#' @seealso [pdiscretestep()], [hazards_to_pmf()], [pmf_to_hazards()],
#'   [fitdistdoublecens()]
#' @export
#' @examples
#' hazards <- c(0.3, 0.5, 1)
#' pdiscretehazard(c(0.5, 1, 2, 3), boundaries = 0:3, hazards = hazards)
pdiscretehazard <- function(q, boundaries = NULL, hazards) {
  if (is.null(boundaries)) {
    boundaries <- seq.int(0L, length(hazards))
  }
  pmf <- hazards_to_pmf(hazards)
  pdiscretestep(q, boundaries, pmf)
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
  logit_h <- alpha_val + sigma * cumsum(c(0, eps_vals))
  h <- 1 / (1 + exp(-logit_h))
  h[K] <- 1
  h
}
attr(pdiscretehazard, "fit_penalty") <- function(par_named, N,
                                                 prior_settings = NULL) {
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

#' Hazard-parameterised piecewise-constant PMF
#'
#' Returns the probability mass for each value in \code{x}. Converts
#' \code{hazards} to a PMF via [hazards_to_pmf()] then delegates to
#' [ddiscretestep()].
#'
#' Outside fitting it is a deterministic wrapper around [ddiscretestep()].
#' It earns its keep as a fitting parameterisation in [fitdistdoublecens()]
#' because the random walk on the logit hazard smooths the recovered PMF.
#'
#' @inheritParams pdiscretehazard
#' @param x Numeric vector of values at which to evaluate the PMF.
#'
#' @return Numeric vector of PMF values, the same length as \code{x}.
#'
#' @family pdiscretehazard
#' @concept pdiscretehazard
#' @seealso [ddiscretestep()], [hazards_to_pmf()], [pmf_to_hazards()],
#'   [fitdistdoublecens()]
#' @export
#' @examples
#' hazards <- c(0.3, 0.5, 1)
#' ddiscretehazard(c(1, 2, 3), boundaries = 0:3, hazards = hazards)
ddiscretehazard <- function(x, boundaries = NULL, hazards) {
  if (is.null(boundaries)) {
    boundaries <- seq.int(0L, length(hazards))
  }
  pmf <- hazards_to_pmf(hazards)
  ddiscretestep(x, boundaries, pmf)
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
#' Outside fitting it is a deterministic wrapper around [rdiscretestep()].
#' It earns its keep as a fitting parameterisation in [fitdistdoublecens()]
#' because the random walk on the logit hazard smooths the recovered PMF.
#'
#' @param n Integer. Number of samples to draw.
#' @inheritParams pdiscretehazard
#'
#' @return Numeric vector of length \code{n}.
#'
#' @family pdiscretehazard
#' @concept pdiscretehazard
#' @seealso [rdiscretestep()], [hazards_to_pmf()], [pmf_to_hazards()],
#'   [fitdistdoublecens()]
#' @export
#' @examples
#' set.seed(42)
#' rdiscretehazard(10, boundaries = 0:3, hazards = c(0.3, 0.5, 1))
rdiscretehazard <- function(n, boundaries = NULL, hazards) {
  if (is.null(boundaries)) {
    boundaries <- seq.int(0L, length(hazards))
  }
  pmf <- hazards_to_pmf(hazards)
  rdiscretestep(n, boundaries, pmf)
}
attr(rdiscretehazard, "name") <- "rdiscretehazard"
attr(rdiscretehazard, "vector_param") <- "hazards"

# Internal helpers (`.sort_eps_names`, `.resolve_hazard_prior`) now live
# in R/nonparametric_helpers.R alongside the rest of the non-parametric
# machinery.
