#' Compute the primary event censored CDF for delays
#'
#' @param q Vector of quantiles
#'
#' @param dist_func Distribution function (CDF)
#'
#' @param pwindow Primary event window
#'
#' @param D Maximum delay (truncation point). If finite, the distribution is
#' truncated at D. If set to Inf, no truncation is applied. Defaults to Inf.
#'
#' @param dprimary Function to generate the probability density function
#' (PDF) of primary event times. This function should take a value `x` and a
#' `pwindow` parameter, and return a probability density. It should be
#' normalized to integrate to 1 over [0, pwindow]. Defaults to a uniform
#' distribution over [0, pwindow]. Users can provide custom functions or use
#' helper functions like `dexpgrowth` for an exponential growth distribution.
#' See `primary_dists.R` for examples.
#'
#' @param dprimary_args List of additional arguments to be passed to
#' dprimary. For example, when using `dexpgrowth`, you would
#' pass `list(min = 0, max = pwindow, r = 0.2)` to set the minimum, maximum,
#' and rate parameters.
#'
#' @param ... Additional arguments to be passed to dist_func
#'
#' @return Vector of primary event censored CDFs, normalized by D if finite
#' (truncation adjustment)
#'
#' @aliases ppcens
#'
#' @examples
#' # Example: Lognormal distribution with uniform primary events
#' pprimarycensoreddist(c(0.1, 0.5, 1), plnorm, meanlog = 0, sdlog = 1)
#'
#' # Example: Lognormal distribution with exponential growth primary events
#' pprimarycensoreddist(
#'   c(0.1, 0.5, 1), plnorm,
#'   dprimary = dexpgrowth,
#'   dprimary_args = list(r = 0.2), meanlog = 0, sdlog = 1
#' )
pprimarycensoreddist <- function(
    q, dist_func, pwindow = 1, D = Inf, dprimary = dunif,
    dprimary_args = list(), ...) {
  check_dist_func(dist_func, D, swindow, ...)

  result <- vapply(q, function(d) {
    if (d <= 0) {
      return(0) # Return 0 for non-positive delays
    } else {
      if (is.infinite(D)) {
        integrand <- function(p) {
          d_adj <- d - p
          dist_func(d_adj, ...) *
            do.call(
              dprimary, c(list(x = p, min = 0, max = pwindow), dprimary_args)
            )
        }
      } else {
        integrand <- function(p) {
          d_adj <- d - p
          D_adj <- D - p
          dist_func(d_adj, ...) / dist_func(D_adj, ...) *
            do.call(
              dprimary, c(list(x = p, min = 0, max = pwindow), dprimary_args)
            )
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
