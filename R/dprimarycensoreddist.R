#' Compute the primary event censored PMF for delays
#'
#' @inheritParams pprimarycensoreddist
#' @param x Vector of quantiles
#' @param swindow Secondary event window (default: 1)
#' @param log Logical; if TRUE, probabilities p are given as log(p)
#' @param ... Additional arguments to be passed to the distribution function
#'
#' @return Vector of primary event censored PMFs, normalized by D if finite
#' (truncation adjustment)
#'
#' @aliases dpcens
#'
#' @examples
#' # Example: Weibull distribution without truncation
#' x <- c(1.0, 2.0, 3.0)
#' pwindow <- 6.0
#' swindow <- 0.5
#' pmf <- dprimarycensoreddist(x, pweibull, shape = 1.5, scale = 2.0)
dprimarycensoreddist <- function(x, dist_func, pwindow = 1, swindow = 1, D,
                                 log = FALSE, ...) {
  result <- vapply(x, function(d) {
    if (d <= 0) {
      return(-Inf) # Return log(0) for non-positive delays
    } else {
      pprimarycensoreddist(d + swindow, dist_func, pwindow, D, ...) -
        pprimarycensoreddist(d, dist_func, pwindow, D, ...)
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
