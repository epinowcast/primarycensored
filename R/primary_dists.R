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
#' @param x Vector of values from [0, pwindow]
#' @param pwindow Primary event window
#' @param rate Rate parameter for the exponential distribution
#'
#' @return Vector of probabilities
#' @export
#'
#' @examples
#' x <- seq(0, 1, by = 0.1)
#' probs <- exp_primary_dist(x, pwindow = 1, rate = 0.2)
exp_primary_dist <- function(x, pwindow, rate) {
  if (rate == 0) {
    return(rep(1 / pwindow, length(x)))
  }
  rate * exp(rate * x) / (exp(rate * pwindow) - 1)
}

#' Linear primary distribution function
#'
#' @param x Vector of values from [0, pwindow]
#' @param pwindow Primary event window
#' @param slope Slope of the linear function (default is 1)
#'
#' @return Vector of probabilities
#' @export
#'
#' @examples
#' x <- seq(0, 7, by = 7)
#' probs <- linear_primary_dist(x, pwindow = 7, slope = 0.5)
linear_primary_dist <- function(x, pwindow, slope = 1) {
  if (slope == 0) {
    return(rep(1 / pwindow, length(x)))
  }
  (2 * slope * x + 2 - slope * pwindow) / (pwindow^2 * slope)
}

#' Normal primary distribution function
#'
#' @param x Vector of values from [0, pwindow]
#' @param pwindow Primary event window
#' @param mean Mean of the normal distribution (as a fraction of pwindow)
#' @param sd Standard deviation of the normal distribution (as a fraction of
#' pwindow)
#'
#' @return Vector of probabilities
#' @export
#'
#' @examples
#' x <- seq(0, 10, by = 1)
#' probs <- normal_primary_dist(x, pwindow = 10, mean = 0.5, sd = 0.2)
normal_primary_dist <- function(x, pwindow, mean, sd) {
  dnorm(x, mean = mean * pwindow, sd = sd * pwindow) /
    (pnorm(pwindow, mean = mean * pwindow, sd = sd * pwindow) -
      pnorm(0, mean = mean * pwindow, sd = sd * pwindow))
}

rboundedgrowth <- function(n, r, pwindow) {
  # Generate uniform random numbers
  u <- runif(n)

  # Use the inverse transform sampling method
  # P(x| 0, pwindow) = r * exp(rx) / (exp(r * pwindow) - 1) # nolint
  if (r != 0) {
    samples <- log(u * (exp(r * pwindow) - 1) + 1) / r
  } else {
    samples <- u * pwindow
  }
  return(samples)
}

dboundedgrowth <- function(x, r, pwindow, log = FALSE) {
  # Check if arguments are within valid ranges
  if (!is.numeric(x) || !is.numeric(r) || !is.numeric(pwindow)) {
    stop("All arguments must be numeric")
  }
  if (pwindow <= 0) {
    stop("pwindow must be positive")
  }
  
  # Use vectorized operations for efficiency
  pdf <- numeric(length(x))
  in_range <- x >= 0 & x <= pwindow
  
  if (abs(r) > 1e-10) {  # Use a small threshold to check if r is close to 0
    pdf[in_range] <- r / (1 - exp(-r * pwindow)) * exp(-r * x[in_range])
  } else {
    pdf[in_range] <- 1 / pwindow
  }
  
  # Set PDF to 0 for x outside [0, pwindow]
  pdf[!in_range] <- 0
  
  # Return log-PDF if requested
  if (log) {
    return(log(pdf))
  } else {
    return(pdf)
  }
}

pboundedgrowth <- function(x, r, pwindow) {
  if (!is.numeric(x) || !is.numeric(r) || !is.numeric(pwindow)) {
    stop("All arguments must be numeric")
  }
  if (pwindow <= 0) {
    stop("pwindow must be positive")
  }
  
  # Use vectorized operations for efficiency
  cdf <- numeric(length(x))
  in_range <- x >= 0 & x <= pwindow
  
  if (r != 0) {
    cdf[in_range] <- (exp(r * x[in_range]) - 1) / (exp(r * pwindow) - 1)
  } else {
    cdf[in_range] <- x[in_range] / pwindow
  }
  
  # Ensure CDF is 1 for x > pwindow
  cdf[x > pwindow] <- 1
  
  return(cdf)
}