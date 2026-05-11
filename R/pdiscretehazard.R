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
# Map named scalar free parameters (alpha, log_sigma, eps_1, ..., eps_{K-1})
# into a length-K hazard vector with the final entry pinned to 1. The
# random-walk model places `logit(h) = alpha + sigma * cumsum(eps)`.
attr(pdiscretehazard, "param_transform") <- .make_hazard_transform("rw")
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
attr(ddiscretehazard, "param_transform") <- attr(
  pdiscretehazard, "param_transform"
)
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

#' Hazard-parameterised step CDF aliases (RW and RE variants)
#'
#' Thin aliases of [pdiscretehazard()]/[ddiscretehazard()]/[rdiscretehazard()].
#' Outside fitting they are byte-identical: the dist family is fully
#' specified by the hazard vector. They differ at the fitting layer in
#' [fitdistdoublecens()]:
#'
#' - `pdiscretehazardrw` uses the random-walk transform
#'   `logit(h_i) = alpha + sigma * cumsum(eps)` (default for
#'   `distr = "discretehazard"`).
#' - `pdiscretehazardre` uses an IID logit random-effect transform
#'   `logit(h_i) = alpha + sigma * eps_i`.
#'
#' @inheritParams pdiscretehazard
#'
#' @return Numeric vector of CDF (or PMF, or sample) values; see the
#'   underlying function.
#'
#' @family pdiscretehazard
#' @export
#' @examples
#' hazards <- c(0.3, 0.5, 1)
#' pdiscretehazardre(c(0.5, 1, 2, 3), boundaries = 0:3, hazards = hazards)
pdiscretehazardrw <- pdiscretehazard
attr(pdiscretehazardrw, "name") <- "pdiscretehazardrw"

#' @rdname pdiscretehazardrw
#' @export
ddiscretehazardrw <- ddiscretehazard
attr(ddiscretehazardrw, "name") <- "ddiscretehazardrw"

#' @rdname pdiscretehazardrw
#' @export
rdiscretehazardrw <- rdiscretehazard
attr(rdiscretehazardrw, "name") <- "rdiscretehazardrw"

#' @rdname pdiscretehazardrw
#' @export
pdiscretehazardre <- pdiscretehazard
attr(pdiscretehazardre, "name") <- "pdiscretehazardre"
attr(pdiscretehazardre, "param_transform") <- .make_hazard_transform("re")

#' @rdname pdiscretehazardrw
#' @export
ddiscretehazardre <- ddiscretehazard
attr(ddiscretehazardre, "name") <- "ddiscretehazardre"
attr(ddiscretehazardre, "param_transform") <- .make_hazard_transform("re")

#' @rdname pdiscretehazardrw
#' @export
rdiscretehazardre <- rdiscretehazard
attr(rdiscretehazardre, "name") <- "rdiscretehazardre"

# Internal helpers (`.sort_eps_names`, `.resolve_hazard_prior`) now live
# in R/nonparametric_helpers.R alongside the rest of the non-parametric
# machinery.
