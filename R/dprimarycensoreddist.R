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
#' where \eqn{F_{\text{cens}}} is the primary event censored CDF. For the
#' explanation and mathematical details of the CDF, refer to the documentation
#' of [pprimarycensoreddist()].
#'
#' @family primarycensoreddist
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
    dprimary_args = list(), log = FALSE, ...) {
  check_pdist(pdist, D, ...)
  check_dprimary(dprimary, pwindow, dprimary_args)

  if (max(x + swindow) > D) {
    stop(
      "Upper truncation point is greater than D. It is ", max(x + swindow),
      " and D is ", D, ". Resolve this by increasing D to be the maximum",
      " of x + swindow."
    )
  }

  result <- vapply(x, function(d) {
    if (d < 0) {
      return(0) # Return log(0) for non-positive delays
    } else {
      pprimarycensoreddist(
        d + swindow, pdist, pwindow, D, dprimary, dprimary_args, ...
      ) -
        pprimarycensoreddist(
          d, pdist, pwindow, D, dprimary, dprimary_args, ...
        )
    }
  }, numeric(1))

  if (log) {
    return(log(result))
  } else {
    return(result)
  }
}

#' @rdname dprimarycensoreddist
#' @export
dpcens <- dprimarycensoreddist
