#' Generate random samples from a primary event censored distribution
#'
#' @inheritParams pprimarycensoreddist
#' @param n Number of random samples to generate
#' @param rprimary Function to generate random samples from the primary
#'   distribution (default is \code{runif})
#' @param rprimary_args List of additional arguments to be passed to rprimary
#' @param ... Additional arguments to be passed to the distribution function
#'
#' @return Vector of random samples from the primary event censored distribution
#'
#' @aliases rpcens
#'
#' @examples
#' # Example: Lognormal distribution with uniform primary events
#' rprimarycensoreddist(10, rlnorm, meanlog = 0, sdlog = 1)
#'
#' # Example: Lognormal distribution with exponential growth primary events
#' rprimarycensoreddist(
#'   10, rlnorm,
#'   rprimary = rexpgrowth, rprimary_args = list(r = 0.2),
#'   meanlog = 0, sdlog = 1
#' )
rprimarycensoreddist <- function(n, dist_func, pwindow = 1, swindow = 1,
                                 D = Inf, rprimary = stats::runif,
                                 rprimary_args = list(), ...) {
  samples <- numeric(n)
  i <- 1

  while (i <= n) {
    # Generate primary event time
    p <- do.call(
      rprimary, c(list(1), rprimary_args, list(min = 0, max = pwindow))
    )

    # Generate delay from the specified distribution
    delay <- dist_func(1, ...)

    # Calculate total delay
    total_delay <- p + delay

    # Round to the nearest swindow
    rounded_delay <- ceiling(total_delay / swindow) * swindow

    # Accept the sample if it's within the valid range (apply truncation)
    if (rounded_delay > 0 && rounded_delay <= D) {
      samples[i] <- rounded_delay
      i <- i + 1
    }
  }

  return(samples)
}

#' @rdname rprimarycensoreddist
#' @export
rpcens <- rprimarycensoreddist
