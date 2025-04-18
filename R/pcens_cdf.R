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

  # Ensure the result is not greater than 1 (accounts for numerical errors)
  result <- pmin(1, result)

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

    # Compute necessary survival and distribution functions
    pgamma_q <- partial_pgamma(non_zero_q)
    pgamma_q_pwindow <- partial_pgamma(non_zero_q + pwindow)
    pgamma_q_1 <- partial_pgamm_k_1(non_zero_q)
    pgamma_q_pwindow_1 <- partial_pgamm_k_1(non_zero_q + pwindow)

    Q_T <- 1 - pgamma_q_pwindow
    Delta_F_T_kp1 <- pgamma_q_pwindow_1 - pgamma_q_1
    Delta_F_T_k <- pgamma_q_pwindow - pgamma_q

    # Calculate Q_Splus using the analytical formula
    Q_Splus <- Q_T +
      (shape * scale / pwindow) * Delta_F_T_kp1 -
      (non_zero_q / pwindow) * Delta_F_T_k

    # Compute the CDF as 1 - Q_Splus
    non_zero_result <- 1 - Q_Splus

    # Assign non-zero results back to the main result vector
    result[!zero_cases] <- non_zero_result
  }

  # Ensure the result is not greater than 1 (accounts for numerical errors)
  result <- pmin(1, result)

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

    # Compute necessary survival and distribution functions
    plnorm_q <- partial_plnorm(non_zero_q)
    plnorm_q_pwindow <- partial_plnorm(non_zero_q + pwindow)
    plnorm_q_sigma2 <- partial_plnorm_sigma2(non_zero_q)
    plnorm_q_pwindow_sigma2 <- partial_plnorm_sigma2(
      non_zero_q + pwindow
    )

    Q_T <- 1 - plnorm_q_pwindow
    Delta_F_T_mu_sigma <- plnorm_q_pwindow_sigma2 - plnorm_q_sigma2
    Delta_F_T <- plnorm_q_pwindow - plnorm_q

    # Calculate Q_Splus using the analytical formula
    Q_Splus <- Q_T +
      (exp(mu + 0.5 * sigma^2) / pwindow) * Delta_F_T_mu_sigma -
      (non_zero_q / pwindow) * Delta_F_T

    # Compute the CDF as 1 - Q_Splus
    non_zero_result <- 1 - Q_Splus

    # Assign non-zero results back to the main result vector
    result[!zero_cases] <- non_zero_result
  }

  # Ensure the result is not greater than 1 (accounts for numerical errors)
  result <- pmin(1, result)

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

  if (pwindow > 3) {
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

  g <- function(t) {
    # Use the lower incomplete gamma function
    scaled_t <- (t * inv_scale)^shape
    g_out <- vapply(
      scaled_t,
      function(x) {
        a <- 1 + inv_shape
        if (abs(-x + a * log(x)) > 700 || abs(a) > 170) {
          return(0)
        }
        result <- pracma::gammainc(x, a)["lowinc"]
        return(result)
      },
      numeric(1)
    )
    return(g_out)
  }

  # Adjust q so that we have [q-pwindow, q]
  q <- q - pwindow

  # Handle cases where q + pwindow <= 0
  zero_cases <- q + pwindow <= 0
  result <- ifelse(zero_cases, 0, NA)

  # Process non-zero cases only if there are any
  if (!all(zero_cases)) {
    non_zero_q <- q[!zero_cases]
    q_pwindow <- pmax(non_zero_q + pwindow, 0)
    non_zero_q_pos <- pmax(non_zero_q, 0)

    # Compute necessary survival and distribution functions
    pweibull_q <- partial_pweibull(non_zero_q_pos)
    pweibull_q_pwindow <- partial_pweibull(q_pwindow)
    g_q <- g(non_zero_q_pos)
    g_q_pwindow <- g(q_pwindow)

    Q_T <- 1 - pweibull_q_pwindow
    Delta_g <- g_q_pwindow - g_q
    Delta_F_T <- pweibull_q_pwindow - pweibull_q

    # Calculate Q_Splus using the analytical formula
    Q_Splus <- Q_T +
      (scale / pwindow) * Delta_g -
      (non_zero_q / pwindow) * Delta_F_T

    # Compute the CDF as 1 - Q_Splus
    non_zero_result <- 1 - Q_Splus

    # Assign non-zero results back to the main result vector
    result[!zero_cases] <- non_zero_result
  }

  # Ensure the result is not greater than 1 (accounts for numerical errors)
  result <- pmin(1, result)

  return(result)
}
