#' Compute the primary event censored CDF for delays
#'
#' @param q Vector of quantiles
#' @param dist_func Distribution function (CDF)
#' @param pwindow Primary event window
#' @param D Maximum delay (truncation point). If finite, the distribution is
#' truncated at D. If set to Inf, no truncation is applied.
#' @param ... Additional arguments to be passed to dist_func
#'
#' @return Vector of primary event censored CDFs, normalized by D if finite
#' (truncation adjustment)
#'
#' @aliases ppcens
#'
#' @examples
#' # Example: Lognormal distribution with truncation
#' q <- c(1.0, 2.0, 3.0)
#' pwindow <- 7.0
#' D <- 10.0 # truncation point
#' cdf <- pprimarycensoreddist(q, plnorm, pwindow, D, meanlog = 0, sdlog = 1)
pprimarycensoreddist <- function(q, dist_func, pwindow = 1, D, ...) {
  result <- vapply(q, function(d) {
    if (d <= 0) {
      return(0) # Return 0 for non-positive delays
    } else {
      if (is.finite(D)) {
        integrand <- function(p) {
          d_adj <- d - p
          D_adj <- D - p
          dist_func(d_adj, ...) / dist_func(D_adj, ...)
        }
      } else {
        integrand <- function(p) {
          d_adj <- d - p
          dist_func(d_adj, ...)
        }
      }

      stats::integrate(integrand, lower = 0, upper = pwindow)$value
    }
  }, numeric(1))

  return(result)
}

#' @rdname pprimarycensoreddist
#' @export
ppcens <- pprimarycensoreddist
