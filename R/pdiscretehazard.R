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
# `param_transform` for the hazard families is resolved at fit time from
# the `hazard_model` argument of [fitdistdoublecens()]; both the random
# walk and the IID random-effect variants map named scalar free
# parameters (alpha, log_sigma, eps_1, ..., eps_{K-1}) into a length-K
# hazard vector with the final entry pinned to 1.
attr(pdiscretehazard, "fit_penalty") <- .hazard_fit_penalty
# Default `lower`/`upper` bounds used by `fitdistdoublecens()` for the
# logit-hazard random walk: `alpha` unbounded, `log_sigma` in a narrow
# strip around 0, `eps_*` clipped to keep the optimiser away from
# absorbing barriers.
attr(pdiscretehazard, "fit_bounds") <- function(par_names) {
  n_eps <- length(par_names) - 2L
  list(
    lower = c(-Inf, -6, rep(-5, n_eps)),
    upper = c(Inf, 4, rep(5, n_eps))
  )
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
attr(ddiscretehazard, "fit_bounds") <- attr(pdiscretehazard, "fit_bounds")

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

#' Start values for the logit-hazard parameterisation
#'
#' Builds a named list of starting values for [fitdistdoublecens()] with
#' \code{distr = "discretehazard"}. The free parameters are the logit
#' intercept \code{alpha}, the log random-walk (or random-effect) scale
#' \code{log_sigma}, and the \code{K - 1} innovations \code{eps_1, ...,
#' eps_{K-1}}. The final bin hazard is pinned to 1 inside the
#' parameterisation so the implied PMF sums to 1.
#'
#' @param K Integer, number of bins in the hazard parameterisation.
#' @param alpha Numeric, start value for the logit intercept.
#' @param log_sigma Numeric, start value for the log scale.
#' @param eps Numeric, scalar or length-\code{K - 1} vector of start
#'   values for the innovations \code{eps_1, ..., eps_{K-1}}.
#'
#' @return Named list suitable for the \code{start} argument of
#'   [fitdistdoublecens()].
#'
#' @family pdiscretehazard
#' @concept pdiscretehazard
#' @export
#' @examples
#' discretehazard_start(K = 5)
discretehazard_start <- function(K, alpha = -2, log_sigma = log(1),
                                 eps = 0) {
  K <- as.integer(K)
  if (!is.finite(K) || K < 1L) {
    stop("`K` must be a positive integer.", call. = FALSE)
  }
  n_eps <- K - 1L
  if (length(eps) == 1L) {
    eps_vec <- rep(eps, n_eps)
  } else if (length(eps) == n_eps) {
    eps_vec <- eps
  } else {
    stop(
      "`eps` must be a scalar or have length K - 1 (got ",
      length(eps), ", expected ", n_eps, ").",
      call. = FALSE
    )
  }
  eps_list <- if (n_eps > 0L) {
    stats::setNames(
      as.list(eps_vec),
      paste0("eps_", seq_len(n_eps))
    )
  } else {
    list()
  }
  c(list(alpha = alpha, log_sigma = log_sigma), eps_list)
}

# Internal helpers (`.sort_eps_names`, `.resolve_hazard_prior`) now live
# in R/nonparametric_helpers.R alongside the rest of the non-parametric
# machinery.
