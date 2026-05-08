#' Step (piecewise-constant) CDF
#'
#' Returns the CDF of a discrete distribution whose mass is concentrated at
#' the right boundary of each bin. The CDF is right-continuous and piecewise
#' constant with jumps at \code{boundaries[2], ..., boundaries[K+1]}.
#'
#' Below \code{boundaries[1]} the function returns 0. At
#' \code{boundaries[i+1]} (the right edge of bin \emph{i}), F jumps by
#' \code{pmf[i]}, so \eqn{F(boundaries[i+1]) = \sum_{j=1}^{i} pmf_j}.
#' For \code{q} in \eqn{[boundaries[i], boundaries[i+1])}, F equals
#' \eqn{\sum_{j=1}^{i-1} pmf_j}. At or above \code{boundaries[K+1]}
#' the function returns 1.
#'
#' ## Use with \code{fitdistdoublecens()}
#'
#' This function carries the attribute \code{vector_param = "pmf"} so that
#' \code{\link{fitdistdoublecens}} can drive it from a flat list of scalar
#' parameters \code{p1, ..., p_{K-1}}. The free parameters are the first
#' \eqn{K-1} bin probabilities; the last is set to
#' \code{1 - sum(p1, ..., p_{K-1})}. When the implied probabilities violate
#' the simplex (any negative entry, or sum departing from 1 by more than
#' \eqn{10^{-8}}), the function returns 0 (or near-zero density in
#' \code{\link{ddiscretestep}}) rather than erroring; this drives the
#' optimiser back to the feasible region.
#'
#' @param q Numeric vector of quantiles.
#' @param boundaries Numeric vector of length \eqn{K+1} defining the bin
#'   edges. Must be strictly increasing. Defaults to \code{0:K} (unit-width
#'   daily bins) where \code{K} is inferred from \code{length(pmf)}.
#' @param pmf Numeric vector of length \eqn{K} giving the probability mass
#'   for each bin. Must be non-negative and sum to approximately 1; if
#'   either condition is violated the function returns a vector of zeros
#'   (a soft simplex penalty for use inside optimisation).
#'
#' @return Numeric vector of CDF values, the same length as \code{q}.
#'
#' @family pdiscretestep
#' @concept pdiscretestep
#' @export
#' @examples
#' # Two-bin PMF: mass 0.3 at x=1, mass 0.7 at x=2
#' pdiscretestep(c(0.5, 1, 1.5, 2), boundaries = 0:2, pmf = c(0.3, 0.7))
pdiscretestep <- function(q, boundaries = NULL, pmf) {
  if (is.null(boundaries)) {
    boundaries <- seq.int(0L, length(pmf))
  }
  .validate_discretestep_args(boundaries, pmf)
  if (!.is_valid_simplex(pmf)) {
    return(rep(0, length(q)))
  }
  K <- length(pmf)
  cum_pmf <- cumsum(pmf)
  # Mass pmf[i] is located at boundaries[i+1] (right edge of bin i).
  # CDF is right-continuous: F(q) = sum(pmf[j] for boundaries[j+1] <= q).
  # Use findInterval on the right edges boundaries[2..K+1].
  right_edges <- boundaries[-1L]
  # findInterval with left.open=FALSE: idx counts right_edges <= q
  idx <- findInterval(q, right_edges, left.open = FALSE)
  # Clamp to [1, K] for safe vector indexing (R drops idx == 0)
  safe_idx <- pmin(pmax(idx, 1L), K)
  result <- ifelse(idx == 0L, 0, cum_pmf[safe_idx])
  result
}
attr(pdiscretestep, "name") <- "pdiscretestep"
attr(pdiscretestep, "vector_param") <- "pmf"

#' Step (piecewise-constant) PMF
#'
#' Returns the probability mass for each value in \code{x}. Mass
#' \code{pmf[i]} is located at \code{boundaries[i+1]} (the right edge of
#' bin \emph{i}); all other values return 0.
#'
#' Like \code{\link{pdiscretestep}}, this function applies a soft simplex
#' penalty: if \code{pmf} contains negative entries or fails to sum to 1
#' (within \eqn{10^{-8}}), the function returns near-zero density
#' (\code{.Machine$double.eps}) rather than erroring. This makes it safe to
#' call from inside \code{fitdistrplus::fitdist()} closures driven by
#' \code{\link{fitdistdoublecens}}.
#'
#' @inheritParams pdiscretestep
#' @param x Numeric vector of values at which to evaluate the PMF.
#' @param pmf Numeric vector of length \eqn{K} giving the probability mass
#'   for each bin. Must be non-negative and sum to approximately 1; if
#'   either condition is violated the function returns a vector of
#'   \code{.Machine$double.eps} rather than 0 (a soft simplex penalty
#'   that keeps log-density finite inside fitting closures).
#'
#' @return Numeric vector of PMF values, the same length as \code{x}.
#'
#' @family pdiscretestep
#' @concept pdiscretestep
#' @export
#' @examples
#' ddiscretestep(c(0, 1, 2, 3), boundaries = 0:3, pmf = c(0.2, 0.5, 0.3))
ddiscretestep <- function(x, boundaries = NULL, pmf) {
  if (is.null(boundaries)) {
    boundaries <- seq.int(0L, length(pmf))
  }
  .validate_discretestep_args(boundaries, pmf)
  if (!.is_valid_simplex(pmf)) {
    return(rep(.Machine$double.eps, length(x)))
  }
  right_edges <- boundaries[-1L]
  # Vectorised lookup: match each x against the right-edge vector and pull
  # the corresponding pmf entry; non-matches return 0.
  idx <- match(x, right_edges)
  result <- ifelse(is.na(idx), 0, pmf[idx])
  result
}
attr(ddiscretestep, "name") <- "ddiscretestep"
attr(ddiscretestep, "vector_param") <- "pmf"

#' Sample from a step distribution
#'
#' Draws \code{n} independent samples from the discrete distribution with
#' probability mass \code{pmf[i]} at \code{boundaries[i+1]}.
#'
#' @param n Integer. Number of samples to draw.
#' @inheritParams pdiscretestep
#'
#' @return Numeric vector of length \code{n}.
#'
#' @family pdiscretestep
#' @concept pdiscretestep
#' @export
#' @examples
#' set.seed(42)
#' rdiscretestep(10, boundaries = 0:3, pmf = c(0.2, 0.5, 0.3))
rdiscretestep <- function(n, boundaries = NULL, pmf) {
  if (is.null(boundaries)) {
    boundaries <- seq.int(0L, length(pmf))
  }
  .validate_discretestep_args(boundaries, pmf)
  right_edges <- boundaries[-1]
  sample(right_edges, size = n, replace = TRUE, prob = pmf)
}
attr(rdiscretestep, "name") <- "rdiscretestep"
attr(rdiscretestep, "vector_param") <- "pmf"

#' Convert discrete-time hazards to a PMF
#'
#' Given a vector of discrete-time hazards \eqn{h_1, \ldots, h_K}, computes
#' the corresponding PMF:
#' \deqn{pmf_i = h_i \prod_{j < i} (1 - h_j)}
#'
#' The last hazard must equal 1 (to ensure the PMF sums to 1). If a vector
#' of length \eqn{K-1} is supplied (all hazards except the final exit
#' hazard), 1 is appended automatically.
#'
#' @param hazards Numeric vector of hazards in \eqn{[0, 1]}. The last entry
#'   must equal 1 (within \eqn{10^{-8}}), or a vector of length \eqn{K-1}
#'   may be supplied and the trailing 1 will be appended.
#'
#' @return Numeric vector of PMF values.
#'
#' @family pdiscretestep
#' @concept pdiscretestep
#' @export
#' @examples
#' hazards_to_pmf(c(0.2, 0.3, 1))
hazards_to_pmf <- function(hazards) {
  if (!is.numeric(hazards) || anyNA(hazards) ||
    any(hazards < 0) || any(hazards > 1)) {
    stop("hazards must be numeric values in [0, 1].", call. = FALSE)
  }
  K <- length(hazards)
  if (K == 0L) stop("hazards must have length >= 1.", call. = FALSE)
  tol <- 1e-8
  if (abs(hazards[K] - 1) > tol) {
    # Assume user omitted the final absorbing hazard
    hazards <- c(hazards, 1)
    K <- K + 1L
  }
  survival <- c(1, cumprod(1 - hazards[-K]))
  pmf <- hazards * survival
  pmf
}

#' Convert a PMF to discrete-time hazards
#'
#' Inverts [hazards_to_pmf()]: given a PMF, computes the discrete-time
#' conditional hazard at each time point.
#' \deqn{h_i = pmf_i / (1 - \sum_{j < i} pmf_j)}
#'
#' The returned vector has the same length as \code{pmf}, with the last
#' entry equal to 1.
#'
#' @param pmf Numeric vector of probabilities. Must be non-negative and
#'   sum to approximately 1.
#'
#' @return Numeric vector of hazards in \eqn{[0, 1]}.
#'
#' @family pdiscretestep
#' @concept pdiscretestep
#' @export
#' @examples
#' pmf_to_hazards(c(0.2, 0.3, 0.5))
pmf_to_hazards <- function(pmf) {
  if (!is.numeric(pmf) || anyNA(pmf) || any(pmf < 0)) {
    stop("pmf must be a non-negative numeric vector.", call. = FALSE)
  }
  tol <- 1e-6
  if (abs(sum(pmf) - 1) > tol) {
    stop(
      "pmf must sum to 1 (got ", round(sum(pmf), 6), ").",
      call. = FALSE
    )
  }
  K <- length(pmf)
  survival <- 1 - c(0, cumsum(pmf[-K]))
  # Trailing zeros in `pmf` give `survival = 0` and `pmf / survival = NaN`.
  # Define the conditional hazard at exhausted-mass bins as 0 instead.
  hazards <- ifelse(survival > 0, pmf / survival, 0)
  # Clamp last entry to exactly 1 (avoid floating-point deviation)
  hazards[K] <- 1
  hazards
}

# Internal helpers (`.is_valid_simplex`, `.validate_discretestep_args`)
# now live in R/nonparametric_helpers.R alongside the rest of the
# non-parametric machinery.
