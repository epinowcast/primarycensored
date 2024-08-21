#' Uniform primary distribution function
#'
#' @param x Vector of values from `runif(n, 0, pwindow)`
#' @param pwindow Primary event window
#'
#' @return Vector of probabilities
#' @export
#'
#' @examples
#' x <- runif(10, 0, 7)
#' probs <- unif_primary_dist(x, pwindow = 7)
unif_primary_dist <- function(x, pwindow) {
  rep(1 / pwindow, length(x))
}

#' Exponential primary distribution function
#'
#' @param x Vector of values i.e. sampled from `runif(n, 0, pwindow)`
#' @param pwindow Primary event window
#' @param rate Rate parameter for the exponential distribution
#'
#' @return Vector of probabilities
#' @export
#'
#' @examples
#' x <- runif(10, 0, 7)
#' probs <- exp_primary_dist(x, pwindow = 7, rate = 0.2)
exp_primary_dist <- function(x, pwindow, rate) {
  if (rate == 0) {
    return(rep(1 / pwindow, length(x)))
  }
  rate * exp(rate * x) / (exp(rate * pwindow) - 1)
}

#' Linear primary distribution function
#'
#' @param x Vector of values from `runif(n, 0, pwindow)`
#' @param pwindow Primary event window
#' @param slope Slope of the linear function (default is 1)
#'
#' @return Vector of probabilities
#' @export
#'
#' @examples
#' x <- runif(10, 0, 7)
#' probs <- linear_primary_dist(x, pwindow = 7, slope = 0.5)
linear_primary_dist <- function(x, pwindow, slope = 1) {
  if (slope == 0) {
    return(rep(1 / pwindow, length(x)))
  }
  (2 * slope * x + 2 - slope * pwindow) / (pwindow^2 * slope)
}

#' Normal primary distribution function
#'
#' @param x Vector of values from `runif(n, 0, pwindow)`
#' @param pwindow Primary event window
#' @param mean Mean of the normal distribution (as a fraction of pwindow)
#' @param sd Standard deviation of the normal distribution (as a fraction of
#' pwindow)
#'
#' @return Vector of probabilities
#' @export
#'
#' @examples
#' x <- runif(10, 0, 10)
#' probs <- normal_primary_dist(x, pwindow = 10, mean = 0.5, sd = 0.2)
normal_primary_dist <- function(x, pwindow, mean, sd) {
  dnorm(x, mean = mean * pwindow, sd = sd * pwindow) /
    (pnorm(pwindow, mean = mean * pwindow, sd = sd * pwindow) -
      pnorm(0, mean = mean * pwindow, sd = sd * pwindow))
}
