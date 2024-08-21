#' Generate random samples from a primary event censored distribution
#'
#' @inheritParams pprimarycensoreddist
#' @param n Number of random samples to generate
#' @inheritParams pprimarycensoreddist
#' @inheritParams dprimarycensoreddist
#' @param primary_dist Primary distribution function (default is \code{runif})
#' @param primary_args List of additional arguments to be passed to the primary
#'   distribution function
#' @param ... Additional arguments to be passed to the distribution function
#'
#' @return Vector of random samples from the primary event censored distribution
#'
#' @aliases rpcens
#'
#' @examples
#' # Example: Lognormal distribution with truncation
#' n <- 1000
#' D <- 10.0 # truncation point
#' samples <- rprimarycensoreddist(n, rlnorm, D = 10, meanlog = 0, sdlog = 1)
rprimarycensoreddist <- function(n, dist_func, pwindow = 1, swindow = 1,
                                 D = Inf, primary_dist = stats::runif,
                                 primary_args = list(), ...) {
  samples <- numeric(n)
  i <- 1

  while (i <= n) {
    # Generate primary event time
    p <- do.call(primary_dist, c(
      list(1), primary_args,
      list(min = 0, max = pwindow)
    ))

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
