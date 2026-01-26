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
#' If truncation is applied (finite D or L > 0), the PMF is normalized to
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
    L = 0,
    D = Inf,
    dprimary = stats::dunif,
    dprimary_args = list(),
    log = FALSE,
    ...) {
  .check_truncation_bounds(L, D)

  check_pdist(pdist, D, ...)
  check_dprimary(dprimary, pwindow, dprimary_args)

  if (max(x + swindow) > D) {
    stop(
      "Upper truncation point is greater than D. It is ",
      max(x + swindow),
      " and D is ",
      D,
      ". Resolve this by increasing D to be the maximum",
      " of x + swindow.",
      call. = FALSE
    )
  }

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

  # Compute CDFs for all unique points
  unique_points <- sort(unique(c(x, x + swindow)))
  unique_points <- unique_points[unique_points > 0]
  if (length(unique_points) == 0) {
    return(rep(0, length(x)))
  }

  cdfs <- pprimarycensored(
    unique_points,
    pdist,
    pwindow,
    Inf,
    L = 0,
    dprimary,
    dprimary_args,
    ...
  )

  # Create a lookup table for CDFs
  cdf_lookup <- stats::setNames(cdfs, as.character(unique_points))

  result <- vapply(
    x,
    function(d) {
      if (d == 0) {
        # Special case for d = 0
        cdf_upper <- cdf_lookup[as.character(swindow)]
        return(cdf_upper)
      } else {
        cdf_upper <- cdf_lookup[as.character(d + swindow)]
        cdf_lower <- cdf_lookup[as.character(d)]
        return(cdf_upper - cdf_lower)
      }
    },
    numeric(1)
  )

  # Apply truncation normalization
  if (is.finite(D) || L > 0) {
    # Get CDF at upper truncation point D
    if (max(unique_points) == D) {
      cdf_D <- max(cdfs)
    } else if (is.infinite(D)) {
      cdf_D <- 1
    } else {
      cdf_D <- pprimarycensored(
        D,
        pdist,
        pwindow,
        Inf,
        L = 0,
        dprimary,
        dprimary_args,
        ...
      )
    }

    # Get CDF at lower truncation point L
    if (L == 0) {
      cdf_L <- 0
    } else if (L %in% unique_points) {
      cdf_L <- cdf_lookup[as.character(L)]
    } else {
      cdf_L <- pprimarycensored(
        L,
        pdist,
        pwindow,
        Inf,
        L = 0,
        dprimary,
        dprimary_args,
        ...
      )
    }

    # Normalize by (F(D) - F(L))
    normaliser <- cdf_D - cdf_L
    result <- result / normaliser
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
