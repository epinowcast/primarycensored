#' Compute the primary event censored PMF for delays
#'
#'
#' This function computes the primary event censored probability mass function
#' (PMF) for a given set of quantiles. It adjusts the PMF of the primary event
#' distribution by accounting for the delay distribution and potential
#' truncation at a maximum delay (D) and minimum delay (L). The function allows
#' for custom primary event distributions and delay distributions.
#'
#' @inheritParams pprimarycensored
#'
#' @param x Vector of quantiles
#'
#' @param swindow Secondary event window (default: 1)
#'
#' @param log Logical; if TRUE, probabilities p are given as log(p)
#'
#' @param ... Additional arguments to be passed to the distribution function
#'
#' @return Vector of primary event censored PMFs, normalized over \[L, D\] if
#' truncation is applied
#'
#' @aliases dpcens
#'
#' @importFrom stats dunif
#'
#' @export
#'
#' @details
#' The primary event censored PMF is computed by taking the difference of the
#' primary event censored cumulative distribution function (CDF) at two points,
#' \eqn{d + \text{swindow}} and \eqn{d}. The primary event censored PMF,
#' \eqn{f_{\text{cens}}(d)}, is given by:
#' \deqn{
#' f_{\text{cens}}(d) = F_{\text{cens}}(d + \text{swindow}) - F_{\text{cens}}(d)
#' }
#' where \eqn{F_{\text{cens}}} is the primary event censored CDF.
#'
#' The function first computes the CDFs for all unique points (including both
#' \eqn{d} and \eqn{d + \text{swindow}}) using [pprimarycensored()]. It then
#' creates a lookup table for these CDFs to efficiently calculate the PMF for
#' each input value. For delays less than L, the function returns 0.
#'
#' When the secondary censoring interval extends past the upper truncation
#' point (\eqn{d + \text{swindow} > D}) but the lower endpoint satisfies
#' \eqn{d < D}, the upper endpoint is internally clipped to \eqn{D} before
#' evaluating the CDF. The likelihood for such an observation is
#' \eqn{P(X \in [d, \min(d + \text{swindow}, D)] \mid L \le X \le D)}, which
#' equals the usual interval probability when \eqn{d + \text{swindow} \le D}.
#' This avoids erroring when an observation's secondary window straddles the
#' truncation point (relevant for non-parametric delays such as
#' [pdiscretestep()]).
#'
#' Observations with \eqn{d \ge D} are rejected with an error: under the
#' truncation \eqn{X \le D}, no event with latent value \eqn{d \ge D} is
#' observable, and accepting such inputs would otherwise yield a 0/0
#' likelihood.
#'
#' The PMF is normalised to
#' ensure it sums to 1 over the range \[L, D\). This normalization uses:
#' \deqn{
#' f_{\text{cens,norm}}(d) = \frac{f_{\text{cens}}(d)}{
#'   F_{\text{cens}}(D) - F_{\text{cens}}(L)}
#' }
#' where \eqn{f_{\text{cens,norm}}(d)} is the normalized PMF. For the
#' explanation and mathematical details of the CDF, refer to the documentation
#' of [pprimarycensored()].
#'
#' @family primarycensored
#'
#' @importFrom stats setNames
#'
#' @examples
#' # Example: Weibull distribution with uniform primary events
#' dprimarycensored(c(0.1, 0.5, 1), pweibull, shape = 1.5, scale = 2.0)
#'
#' # Example: Weibull distribution with exponential growth primary events
#' dprimarycensored(
#'   c(0.1, 0.5, 1), pweibull,
#'   dprimary = dexpgrowth,
#'   dprimary_args = list(r = 0.2), shape = 1.5, scale = 2.0
#' )
#'
#' # Example: Left-truncated distribution (e.g., for generation intervals)
#' dprimarycensored(1:9, pweibull, L = 1, D = 10, shape = 1.5, scale = 2.0)
dprimarycensored <- function(
    x,
    pdist,
    pwindow = 1,
    swindow = 1,
    L = -Inf,
    D = Inf,
    dprimary = stats::dunif,
    primary_args = NULL,
    pprimary = NULL,
    dprimary_args = NULL,
    log = FALSE,
    ...) {
  .check_truncation_bounds(L, D)

  primary_args <- .resolve_primary_args( # nolint: object_usage_linter
    primary_args, dprimary_args, "dprimarycensored"
  )
  pdist <- .resolve_pdist(pdist, type = "p") # nolint: object_usage_linter
  pprimary <- .resolve_pprimary( # nolint: object_usage_linter
    dprimary, pprimary
  )

  check_pdist(pdist, D = D, ...)
  check_dprimary(dprimary, pwindow, primary_args)

  if (min(x) < L) {
    stop(
      "Some values of x are below L. Minimum x is ",
      min(x),
      " and L is ",
      L,
      ". Resolve this by filtering x to only include values >= L.",
      call. = FALSE
    )
  }

  if (is.finite(D) && max(x) >= D) {
    stop(
      "Upper truncation point is greater than D. Maximum x is ",
      max(x),
      " and D is ",
      D,
      ". Under truncation at D no event with latent value >= D is ",
      "observable; resolve this by filtering x to values strictly less than D.",
      call. = FALSE
    )
  }

  # Clip the upper end of each secondary interval at D so observations with
  # `x + swindow > D` (legitimate when the secondary censoring interval
  # straddles D) are still valid. The likelihood becomes
  # `P(X in [x, min(x + swindow, D)] | L <= X <= D)`, which equals the usual
  # interval probability when `x + swindow <= D` (the parametric default) and
  # captures the residual mass between `x` and `D` otherwise.
  upper_raw <- x + swindow
  upper <- pmin(upper_raw, D)
  if (is.finite(D) && any(upper_raw > D)) {
    message(
      "Upper truncation point is greater than D. It is ",
      max(upper_raw),
      " and D is ",
      D,
      "; clipping the upper end of secondary intervals at D."
    )
  }

  # Compute CDFs for all unique points
  unique_points <- sort(unique(c(x, upper)))
  if (length(unique_points) == 0) {
    return(rep(0, length(x)))
  }

  # Compute raw (unnormalised) CDFs via `L = -Inf, D = Inf` so PMF differences
  # below can be normalised with the truncation-aware F_cens(L) and F_cens(D).
  cdfs <- pprimarycensored(
    unique_points,
    pdist,
    pwindow = pwindow,
    L = -Inf,
    D = Inf,
    dprimary = dprimary,
    primary_args = primary_args,
    pprimary = pprimary,
    ...
  )

  # Create a lookup table for CDFs
  cdf_lookup <- stats::setNames(cdfs, as.character(unique_points))

  result <- vapply(
    seq_along(x),
    function(i) {
      cdf_upper <- cdf_lookup[as.character(upper[i])]
      cdf_lower <- cdf_lookup[as.character(x[i])]
      return(cdf_upper - cdf_lower)
    },
    numeric(1)
  )

  # Fast path: with no truncation on either side the raw PMF needs no
  # renormalisation, so skip the two extra `pprimarycensored` lookups below.
  if (!(is.infinite(L) && is.infinite(D))) {
    # F_cens(D). Reuse the existing `unique_points` lookup when D lands on
    # one of them; otherwise compute on demand.
    if (is.infinite(D)) {
      cdf_D <- 1
    } else if (D %in% unique_points) {
      cdf_D <- cdf_lookup[[as.character(D)]]
    } else {
      cdf_D <- pprimarycensored(
        D,
        pdist,
        pwindow = pwindow,
        L = -Inf,
        D = Inf,
        dprimary = dprimary,
        primary_args = primary_args,
        pprimary = pprimary,
        ...
      )
    }

    # F_cens(L). `L = -Inf` is the "no left truncation" sentinel so we skip
    # the integral; otherwise reuse the lookup when L is already in it.
    if (is.infinite(L)) {
      cdf_L <- 0
    } else if (L %in% unique_points) {
      cdf_L <- cdf_lookup[[as.character(L)]]
    } else {
      cdf_L <- pprimarycensored(
        L,
        pdist,
        pwindow = pwindow,
        L = -Inf,
        D = Inf,
        dprimary = dprimary,
        primary_args = primary_args,
        pprimary = pprimary,
        ...
      )
    }

    # Divide by (F(D) - F(L)). Skip the division when the normaliser is 1
    # (e.g. a finite `L` that sits below the support of the delay, so
    # `F_cens(L) = 0`, paired with `D = Inf` where `F_cens(D) = 1`).
    normaliser <- cdf_D - cdf_L
    if (normaliser != 1) {
      result <- result / normaliser
    }
  }

  # Ensure non-negative values (can become slightly negative due to
  # floating-point precision when computing CDF differences)
  result <- pmax(0, result)

  if (log) {
    return(log(result))
  } else {
    return(result)
  }
}

#' @rdname dprimarycensored
#' @export
dpcens <- dprimarycensored
