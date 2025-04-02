#' Exponential growth distribution functions
#'
#' Density, distribution function, and random generation for the exponential
#' growth distribution.
#'
#' @param x,q Vector of quantiles.
#'
#' @param n Number of observations. If `length(n) > 1`, the length is taken to
#'  be the number required.
#'
#' @param min Minimum value of the distribution range. Default is 0.
#'
#' @param max Maximum value of the distribution range. Default is 1.
#'
#' @param r Rate parameter for the exponential growth.
#'
#' @param log,log.p Logical; if TRUE, probabilities p are given as log(p).
#'
#' @param lower.tail Logical; if TRUE (default), probabilities are P\[X <= x\],
#'  otherwise, P\[X > x\].
#'
#' @return `dexpgrowth` gives the density, `pexpgrowth` gives the distribution
#' function, and `rexpgrowth` generates random deviates.
#'
#' The length of the result is determined by `n` for `rexpgrowth`, and is the
#' maximum of the lengths of the numerical arguments for the other functions.
#'
#' @details
#' The exponential growth distribution is defined on the interval \[min, max\]
#' with rate parameter (r). Its probability density function (PDF) is:
#'
#' \deqn{f(x) = \frac{r \cdot \exp(r \cdot (x - min))}{\exp(r \cdot max) -
#'  \exp(r \cdot min)}}
#'
#' The cumulative distribution function (CDF) is:
#'
#' \deqn{F(x) = \frac{\exp(r \cdot (x - min)) - \exp(r \cdot min)}{
#'  \exp(r \cdot max) - \exp(r \cdot min)}}
#'
#' For random number generation, we use the inverse transform sampling method:
#' 1. Generate \eqn{u \sim \text{Uniform}(0,1)}
#' 2. Set \eqn{F(x) = u} and solve for \eqn{x}:
#'    \deqn{
#'    x = min + \frac{1}{r} \cdot \log(u \cdot (\exp(r \cdot max) -
#'   \exp(r \cdot min)) + \exp(r \cdot min))
#'    }
#'
#' This method works because of the probability integral transform theorem,
#' which states that if \eqn{X} is a continuous random variable with CDF
#' \eqn{F(x)}, then \eqn{Y = F(X)} follows a \eqn{\text{Uniform}(0,1)}
#' distribution. Conversely, if \eqn{U} is a \eqn{\text{Uniform}(0,1)} random
#' variable, then \eqn{F^{-1}(U)} has the same distribution as \eqn{X}, where
#' \eqn{F^{-1}} is the inverse of the CDF.
#'
#' In our case, we generate \eqn{u} from \eqn{\text{Uniform}(0,1)}, then solve
#' \eqn{F(x) = u} for \eqn{x} to get a sample from our exponential growth
#' distribution. The formula for \eqn{x} is derived by algebraically solving
#' the equation:
#'
#' \deqn{
#' u = \frac{\exp(r \cdot (x - min)) - \exp(r \cdot min)}{\exp(r \cdot max) -
#'  \exp(r \cdot min)}
#' }
#'
#' When \eqn{r} is very close to 0 (\eqn{|r| < 1e-10}), the distribution
#' approximates a uniform distribution on \[min, max\], and we use a simpler
#' method to generate samples directly from this uniform distribution.
#'
#' @family primaryeventdistributions
#'
#' @examples
#' x <- seq(0, 1, by = 0.1)
#' probs <- dexpgrowth(x, r = 0.2)
#' cumprobs <- pexpgrowth(x, r = 0.2)
#' samples <- rexpgrowth(100, r = 0.2)
#'
#' @name expgrowth
NULL

#' @rdname expgrowth
#' @export
dexpgrowth <- function(x, min = 0, max = 1, r, log = FALSE) {
  if (abs(r) < 1e-10) {
    result <- rep(1 / (max - min), length(x))
  } else {
    result <- r * exp(r * (x - min)) / (exp(r * max) - exp(r * min))
  }
  result[x < min | x > max] <- 0
  if (log) {
    return(log(result))
  } else {
    return(result)
  }
}

#' @rdname expgrowth
#' @export
pexpgrowth <- function(
    q,
    min = 0,
    max = 1,
    r,
    lower.tail = TRUE,
    log.p = FALSE) {
  cumulative <- numeric(length(q))
  in_range <- q >= min & q <= max

  if (abs(r) < 1e-10) {
    cumulative[in_range] <- (q[in_range] - min) / (max - min)
  } else {
    cumulative[in_range] <- (exp(r * (q[in_range] - min)) - exp(r * min)) /
      (exp(r * max) - exp(r * min))
  }

  cumulative[q > max] <- 1
  cumulative[q < min] <- 0

  if (!lower.tail) cumulative <- 1 - cumulative
  if (log.p) {
    return(log(cumulative))
  } else {
    return(cumulative)
  }
}

#' @rdname expgrowth
#' @importFrom stats runif
#' @export
rexpgrowth <- function(n, min = 0, max = 1, r) {
  u <- runif(n)
  if (abs(r) < 1e-10) {
    samples <- min + u * (max - min)
  } else {
    samples <- log(u * (exp(r * max) - exp(r * min)) + exp(r * min)) / r
  }
  return(samples)
}
