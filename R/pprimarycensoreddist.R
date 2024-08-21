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
#' @param primary_dist Function to generate the probability density function
#' (PDF) of primary event times. This function should take a value `x` and a
#' `pwindow` parameter, and return a probability density. It should be
#' normalized to integrate to 1 over [0, pwindow]. Defaults to a uniform
#' distribution over [0, pwindow]. Users can provide custom functions or use
#' helper functions like `exp_primary_dist` for an exponential distribution.
#' These functions typically implement the inverse transform sampling method
#' for efficient random number generation. See `primary_dists.R` for examples.
#'
#' @param primary_args List of additional arguments to be passed to
#' primary_dist. For example, when using `exp_primary_dist`, you would
#' pass `list(rate = 0.2)` to set the rate parameter.
#'
#' @param ... Additional arguments to be passed to dist_func
#'
#' @return Vector of primary event censored CDFs, normalized by D if finite
#' (truncation adjustment)
#'
#' @aliases ppcens
#'
#' @examples
#' # Example: Lognormal distribution with truncation and uniform primary events
#' q <- c(1.0, 2.0, 3.0)
#' pwindow <- 7.0
#' D <- 10.0 # truncation point
#' cdf <- pprimarycensoreddist(q, plnorm, pwindow, D, meanlog = 0, sdlog = 1)
#'
#' # Example: Lognormal distribution with exponential primary events
#' cdf_exp <- pprimarycensoreddist(q, plnorm, pwindow, D,
#'   primary_dist = exp_primary_dist,
#'   primary_args = list(rate = 0.2),
#'   meanlog = 0, sdlog = 1
#' )
pprimarycensoreddist <- function(q, dist_func, pwindow = 1, D = Inf,
  primary_dist = primarycensoreddist::unif_primary_dist,
  primary_args = list(), ...) {

  result <- vapply(q, function(d) {
    if (d <= 0) {
      return(0) # Return 0 for non-positive delays
    } else {
      if (is.infinite(D)) {
        integrand <- function(p) {
          d_adj <- d - p
          dist_func(d_adj, ...) *
            do.call(
              primary_dist, c(list(x = p, pwindow = pwindow), primary_args)
            )
        }
      } else {
        integrand <- function(p) {
          d_adj <- d - p
          D_adj <- D - p
          dist_func(d_adj, ...) / dist_func(D_adj, ...) *
            do.call(
              primary_dist, c(list(x = p, pwindow = pwindow), primary_args)
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
