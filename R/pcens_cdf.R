#' Compute primary event censored CDF
#'
#' This function dispatches to either analytical solutions (if available) or
#' numerical integration via the default method. To see which combinations have
#' analytical solutions implemented, use `methods(pcens_cdf)`. For example,
#' `pcens_cdf.gamma_unif` indicates an analytical solution exists for gamma
#' delay with uniform primary event distributions.
#'
#' @inheritParams pprimarycensored
#'
#' @param object A `primarycensored` object as created by
#' [new_pcens()].
#'
#' @param use_numeric Logical, if TRUE forces use of numeric integration
#' even for distributions with analytical solutions. This is primarily
#' useful for testing purposes or for settings where the analytical solution
#' breaks down.
#'
#' @return Vector of computed primary event censored CDFs
#'
#' @family pcens
#'
#' @export
pcens_cdf <- function(
    object,
    q,
    pwindow,
    use_numeric = FALSE) {
  UseMethod("pcens_cdf")
}

#' Default method for computing primary event censored CDF
#'
#' This method serves as a fallback for combinations of delay and primary
#' event distributions that don't have specific implementations. It uses
#' a numeric integration method.
#'
#' @inheritParams pcens_cdf
#' @inheritParams pprimarycensored
#'
#' @details
#' This method implements the numerical integration approach for computing
#' the primary event censored CDF. It uses the same mathematical formulation
#' as described in the details section of [pprimarycensored()], but
#' applies numerical integration instead of analytical solutions.
#'
#' @seealso [pprimarycensored()] for the mathematical details of the
#'  primary event censored CDF computation.
#'
#' @family pcens
#'
#' @inherit pcens_cdf return
#'
#' @export
#' @examples
#' # Create a primarycensored object with gamma delay and uniform primary
#' pcens_obj <- new_pcens(
#'   pdist = pgamma,
#'   dprimary = dunif,
#'   dprimary_args = list(min = 0, max = 1),
#'   shape = 3,
#'   scale = 2
#' )
#'
#' # Compute CDF for a single value
#' pcens_cdf(pcens_obj, q = 9, pwindow = 1)
#'
#' # Compute CDF for multiple values
#' pcens_cdf(pcens_obj, q = c(4, 6, 8), pwindow = 1)
pcens_cdf.default <- function(
    object,
    q,
    pwindow,
    use_numeric = FALSE) {
  result <- vapply(
    q,
    function(d) {
      if (d <= 0) {
        return(0) # Return 0 for non-positive delays
      }
      integrand <- function(p) {
        d_adj <- d - p
        do.call(object$pdist, c(list(q = d_adj), object$args)) *
          do.call(
            object$dprimary,
            c(list(x = p, min = 0, max = pwindow), object$dprimary_args)
          )
      }
      return(stats::integrate(integrand, lower = 0, upper = pwindow)$value)
    },
    numeric(1)
  )

  # Ensure the result is in [0, 1] (accounts for numerical errors)
  result <- pmin(1, pmax(0, result))

  return(result)
}

#' Method for Gamma delay with uniform primary
#'
#' @inheritParams pcens_cdf
#'
#' @family pcens
#'
#' @inherit pcens_cdf return
#'
#' @export
pcens_cdf.pcens_pgamma_dunif <- function(
    object,
    q,
    pwindow,
    use_numeric = FALSE) {
  if (isTRUE(use_numeric)) {
    return(
      pcens_cdf.default(object, q, pwindow, use_numeric)
    )
  }
  # Extract Gamma distribution parameters
  shape <- object$args$shape
  scale <- object$args$scale
  rate <- object$args$rate
  # if we don't have scale get fromm rate
  if (is.null(scale) && !is.null(rate)) {
    scale <- 1 / rate
  }
  if (is.null(shape)) {
    stop("shape parameter is required for Gamma distribution", call. = FALSE)
  }
  if (is.null(scale)) {
    stop(
      "scale or rate parameter is required for Gamma distribution",
      call. = FALSE
    )
  }

  partial_pgamma <- function(q) {
    stats::pgamma(q, shape = shape, scale = scale)
  }
  partial_pgamm_k_1 <- function(q) {
    stats::pgamma(q, shape = shape + 1, scale = scale)
  }
  # Adjust q so that we have [q-pwindow, q]
  q <- q - pwindow
  # Handle cases where q + pwindow <= 0
  zero_cases <- q + pwindow <= 0
  result <- ifelse(zero_cases, 0, NA)

  # Process non-zero cases only if there are any
  if (!all(zero_cases)) {
    non_zero_q <- q[!zero_cases]
    d <- non_zero_q + pwindow

    # Compute delay CDF at the interval endpoints and at the shifted (k+1)
    # distribution for the mean-shift term E[T] = shape * scale.
    F_T_q <- partial_pgamma(non_zero_q)
    F_T_d <- partial_pgamma(d)
    F_T_q_kp1 <- partial_pgamm_k_1(non_zero_q)
    F_T_d_kp1 <- partial_pgamm_k_1(d)

    E_T <- shape * scale

    # Direct CDF form:
    #   F_{S+}(d) = ( d F_T(d) - q F_T(q) - E_T (F~_T(d) - F~_T(q)) ) / w_P
    non_zero_result <-
      (d * F_T_d - non_zero_q * F_T_q -
        E_T * (F_T_d_kp1 - F_T_q_kp1)) / pwindow

    # Assign non-zero results back to the main result vector
    result[!zero_cases] <- non_zero_result
  }

  # Ensure the result is in [0, 1] (accounts for numerical errors)
  result <- pmin(1, pmax(0, result))

  return(result)
}

#' Method for Log-Normal delay with uniform primary
#'
#' @inheritParams pcens_cdf
#'
#' @family pcens
#'
#' @inherit pcens_cdf return
#'
#' @export
pcens_cdf.pcens_plnorm_dunif <- function(
    object,
    q,
    pwindow,
    use_numeric = FALSE) {
  if (isTRUE(use_numeric)) {
    return(
      pcens_cdf.default(object, q, pwindow, use_numeric)
    )
  }

  # Extract Log-Normal distribution parameters
  mu <- object$args$meanlog
  sigma <- object$args$sdlog
  if (is.null(mu)) {
    stop(
      "meanlog parameter is required for Log-Normal distribution",
      call. = FALSE
    )
  }
  if (is.null(sigma)) {
    stop(
      "sdlog parameter is required for Log-Normal distribution",
      call. = FALSE
    )
  }

  partial_plnorm <- function(q) {
    stats::plnorm(q, meanlog = mu, sdlog = sigma)
  }
  partial_plnorm_sigma2 <- function(q) {
    stats::plnorm(q, meanlog = mu + sigma^2, sdlog = sigma)
  }
  # Adjust q so that we have [q-pwindow, q]
  q <- q - pwindow

  # Handle cases where q + pwindow <= 0
  zero_cases <- q + pwindow <= 0
  result <- ifelse(zero_cases, 0, NA)

  # Process non-zero cases only if there are any
  if (!all(zero_cases)) {
    non_zero_q <- q[!zero_cases]
    d <- non_zero_q + pwindow

    # Compute delay CDF at the interval endpoints and at the shifted
    # (meanlog + sigma^2) distribution for the mean-shift term
    # E[T] = exp(mu + sigma^2 / 2).
    F_T_q <- partial_plnorm(non_zero_q)
    F_T_d <- partial_plnorm(d)
    F_T_q_shift <- partial_plnorm_sigma2(non_zero_q)
    F_T_d_shift <- partial_plnorm_sigma2(d)

    E_T <- exp(mu + 0.5 * sigma^2)

    # Direct CDF form:
    #   F_{S+}(d) = ( d F_T(d) - q F_T(q) - E_T (F~_T(d) - F~_T(q)) ) / w_P
    non_zero_result <-
      (d * F_T_d - non_zero_q * F_T_q -
        E_T * (F_T_d_shift - F_T_q_shift)) / pwindow

    # Assign non-zero results back to the main result vector
    result[!zero_cases] <- non_zero_result
  }

  # Ensure the result is in [0, 1] (accounts for numerical errors)
  result <- pmin(1, pmax(0, result))

  return(result)
}

#' Method for Weibull delay with uniform primary
#'
#' @inheritParams pcens_cdf
#'
#' @family pcens
#'
#' @inherit pcens_cdf return
#'
#' @export
pcens_cdf.pcens_pweibull_dunif <- function(
    object,
    q,
    pwindow,
    use_numeric = FALSE) {
  if (isTRUE(use_numeric)) {
    return(
      pcens_cdf.default(object, q, pwindow, use_numeric)
    )
  }

  # Extract Weibull distribution parameters
  shape <- object$args$shape
  scale <- object$args$scale
  if (is.null(shape)) {
    stop("shape parameter is required for Weibull distribution", call. = FALSE)
  }
  if (is.null(scale)) {
    stop("scale parameter is required for Weibull distribution", call. = FALSE)
  }

  partial_pweibull <- function(q) {
    stats::pweibull(q, shape = shape, scale = scale)
  }

  # Precompute constants
  inv_shape <- 1 / shape
  inv_scale <- 1 / scale
  a <- 1 + inv_shape
  lgamma_a <- lgamma(a)

  # Lower incomplete gamma gamma(a, x) via the regularised form from
  # stats::pgamma, which is numerically stable for large x where an
  # unregularised series expansion would overflow.
  g <- function(t) {
    x <- (t * inv_scale)^shape
    exp(stats::pgamma(x, shape = a, scale = 1, log.p = TRUE) + lgamma_a)
  }

  # Adjust q so that we have [q-pwindow, q]
  q <- q - pwindow

  # Handle cases where q + pwindow <= 0
  zero_cases <- q + pwindow <= 0
  result <- ifelse(zero_cases, 0, NA)

  # Process non-zero cases only if there are any
  if (!all(zero_cases)) {
    non_zero_q <- q[!zero_cases]
    d <- non_zero_q + pwindow
    # Clamp to zero for evaluating F_T and g (both undefined / zero on R_-).
    # The products q * F_T(q) and scale * g(q) are then zero when q < 0,
    # matching F_T(q) = 0 and g(q) = 0 for q <= 0.
    q_pos <- pmax(non_zero_q, 0)
    d_pos <- pmax(d, 0)

    # Compute delay CDF and helper g at the interval endpoints.
    F_T_q <- partial_pweibull(q_pos)
    F_T_d <- partial_pweibull(d_pos)
    g_q <- g(q_pos)
    g_d <- g(d_pos)

    # Direct CDF form (with E[T] = scale and the shifted-CDF role played
    # by g / scale, so that E[T] * (F~_T(d) - F~_T(q)) = scale * (g(d) - g(q))):
    #   F_{S+}(d) = ( d F_T(d) - q F_T(q) - scale (g(d) - g(q)) ) / w_P
    non_zero_result <-
      (d_pos * F_T_d - q_pos * F_T_q - scale * (g_d - g_q)) / pwindow

    # Assign non-zero results back to the main result vector
    result[!zero_cases] <- non_zero_result
  }

  # Ensure the result is in [0, 1] (accounts for numerical errors)
  result <- pmin(1, pmax(0, result))

  return(result)
}
