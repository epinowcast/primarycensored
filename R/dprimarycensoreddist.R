#' Compute the primary event censored PMF for delays
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
    x, dist_func, pwindow = 1, swindow = 1,
    D = Inf, dprimary = dunif,
    dprimary_args = list(), log = FALSE, ...) {
  result <- vapply(x, function(d) {
    if (d <= 0) {
      return(-Inf) # Return log(0) for non-positive delays
    } else {
      pprimarycensoreddist(
        d + swindow, dist_func, pwindow, D, dprimary, dprimary_args, ...
      ) -
        pprimarycensoreddist(
          d, dist_func, pwindow, D, dprimary, dprimary_args, ...
        )
    }
  }, numeric(1))

  if (log) {
    return(result)
  } else {
    return(exp(result))
  }
}

#' @rdname dprimarycensoreddist
#' @export
dpcens <- dprimarycensoreddist
