#' Compute the primary event censored PMF for delays
#'
#'
#' This function computes the primary event censored probability mass function
#' (PMF) for a given set of quantiles. It adjusts the PMF of the primary event
#' distribution by accounting for the delay distribution and potential
#' truncation at a maximum delay (D). The function allows for custom primary
#' event distributions and delay distributions.
#'
#' @inheritParams pprimarycensoreddist
#'
#' @param x Vector of quantiles
#'
#' @param swindow Secondary event window (default: 1)
#'
#' @param log Logical; if TRUE, probabilities p are given as log(p)
#'
#' @param ... Additional arguments to be passed to the distribution function
#'
#' @return Vector of primary event censored PMFs, normalized by D if finite
#' (truncation adjustment)
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
#' \eqn{d} and \eqn{d + \text{swindow}}) using [pprimarycensoreddist()]. It then
#' creates a lookup table for these CDFs to efficiently calculate the PMF for
#' each input value. For non-positive delays, the function returns 0.
#'
#' If a finite maximum delay \eqn{D} is specified, the PMF is normalized to
#' ensure it sums to 1 over the range \[0, D\]. This normalization can be
#' expressed as:
#' \deqn{
#' f_{\text{cens,norm}}(d) = \frac{f_{\text{cens}}(d)}{\sum_{i=0}^{D-1}
#'  f_{\text{cens}}(i)}
#' }
#' where \eqn{f_{\text{cens,norm}}(d)} is the normalized PMF and
#' \eqn{f_{\text{cens}}(d)} is the unnormalized PMF. For the explanation and
#' mathematical details of the CDF, refer to the documentation of
#' [pprimarycensoreddist()].
#'
#' @family primarycensoreddist
#'
#' @importFrom stats setNames
#'
#' @examples
#' # Example: Weibull distribution with uniform primary events
#' dprimarycensoreddist(c(0.1, 0.5, 1), pweibull, shape = 1.5, scale = 2.0)
#'
#' # Example: Weibull distribution with exponential growth primary events
#' dprimarycensoreddist(
#'   c(0.1, 0.5, 1), pweibull,
#'   dprimary = dexpgrowth,
#'   dprimary_args = list(r = 0.2), shape = 1.5, scale = 2.0
#' )
dprimarycensoreddist <- function(
    x, pdist, pwindow = 1, swindow = 1,
    D = Inf, dprimary = stats::dunif,
    dprimary_args = list(), log = FALSE,
    pdist_name = NULL, dprimary_name = NULL, ...) {
  check_pdist(pdist, D, ...)
  check_dprimary(dprimary, pwindow, dprimary_args)

  if (max(x + swindow) > D) {
    stop(
      "Upper truncation point is greater than D. It is ", max(x + swindow),
      " and D is ", D, ". Resolve this by increasing D to be the maximum",
      " of x + swindow."
    )
  }

  if (is.null(pdist_name)) {
    pdist_name <- .extract_function_name(substitute(pdist))
  }
  if (is.null(dprimary_name)) {
    dprimary_name <- .extract_function_name(substitute(dprimary))
  }

  # Compute CDFs for all unique points
  unique_points <- sort(unique(c(x, x + swindow)))
  unique_points <- unique_points[unique_points > 0]
  if (length(unique_points) == 0) {
    return(rep(0, length(x)))
  }

  cdfs <- pprimarycensoreddist(
    unique_points, pdist, pwindow, Inf, dprimary, dprimary_args,
    pdist_name = pdist_name, dprimary_name = dprimary_name, ...
  )

  # Create a lookup table for CDFs
  cdf_lookup <- stats::setNames(cdfs, as.character(unique_points))

  result <- vapply(x, function(d) {
    if (d < 0) {
      return(0) # Return 0 for negative delays
    } else if (d == 0) {
      # Special case for d = 0
      cdf_upper <- cdf_lookup[as.character(swindow)]
      return(cdf_upper)
    } else {
      cdf_upper <- cdf_lookup[as.character(d + swindow)]
      cdf_lower <- cdf_lookup[as.character(d)]
      return(cdf_upper - cdf_lower)
    }
  }, numeric(1))

  if (is.finite(D)) {
    if (max(unique_points) == D) {
      cdf_D <- max(cdfs)
    } else {
      cdf_D <- pprimarycensoreddist(
        D, pdist, pwindow, Inf, dprimary, dprimary_args, ...
      )
    }
    result <- result / cdf_D
  }

  if (log) {
    return(log(result))
  } else {
    return(result)
  }
}

#' @rdname dprimarycensoreddist
#' @export
dpcens <- dprimarycensoreddist
