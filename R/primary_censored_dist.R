#' S3 class for primary event censored distribution computation
#'
#' @inheritParams pprimarycensoreddist
#'
#' @return An object of class primary_censored_cdf
#'
#' @family primary_censored_dist
#'
#' @export
new_primary_censored_dist <- function(
    pdist, dprimary, dprimary_args,
    pdist_name = NULL,
    dprimary_name = NULL, ...) {
  if (is.null(pdist_name)) {
    pdist_name <- .extract_function_name(substitute(pdist))
  }
  if (is.null(dprimary_name)) {
    dprimary_name <- .extract_function_name(substitute(dprimary))
  }

  structure(
    list(
      pdist = pdist,
      dprimary = dprimary,
      dprimary_args = dprimary_args,
      args = list(...)
    ),
    class = c(
      paste0(
        "pcens_",
        pdist_name, "_",
        dprimary_name
      )
    )
  )
}

#' Compute primary event censored CDF
#'
#' @inheritParams pprimarycensoreddist
#'
#' @param object A `primary_censored_dist` object as created by
#' [new_primary_censored_dist()].
#'
#' @param use_numeric Logical, if TRUE forces use of numeric integration
#' even for distributions with analytical solutions. This is primarily
#' useful for testing purposes or for settings where the analytical solution
#' breaks down.
#'
#' @return Vector of primary event censored CDFs
#'
#' @family primary_censored_dist
#'
#' @export
primary_censored_cdf <- function(
    object, q, pwindow, use_numeric = FALSE) {
  UseMethod("primary_censored_cdf")
}

#' Default method for computing primary event censored CDF
#'
#' This method serves as a fallback for combinations of delay and primary
#' event distributions that don't have specific implementations. It uses
#' the numeric integration method.
#'
#' @inheritParams primary_censored_cdf
#' @inheritParams pprimarycensoreddist
#'
#' @details
#' This method implements the numerical integration approach for computing
#' the primary event censored CDF. It uses the same mathematical formulation
#' as described in the details section of [pprimarycensoreddist()], but
#' applies numerical integration instead of analytical solutions.
#'
#' @seealso [pprimarycensoreddist()] for the mathematical details of the
#' primary event censored CDF computation.
#'
#' @family primary_censored_dist
#'
#' @export
primary_censored_cdf.default <- function(
    object, q, pwindow, use_numeric = FALSE) {
  result <- vapply(q, function(d) {
    if (d <= 0) {
      return(0) # Return 0 for non-positive delays
    } else {
      integrand <- function(p) {
        d_adj <- d - p
        do.call(object$pdist, c(list(q = d_adj), object$args)) *
          do.call(
            object$dprimary,
            c(list(x = p, min = 0, max = pwindow), object$dprimary_args)
          )
      }

      stats::integrate(integrand, lower = 0, upper = pwindow)$value
    }
  }, numeric(1))

  # Ensure the result is not greater than 1 (accounts for numerical errors)
  result <- pmin(1, result)

  return(result)
}

#' Method for Gamma delay with uniform primary
#'
#' @inheritParams primary_censored_cdf
#'
#' @family primary_censored_dist
#'
#' @export
primary_censored_cdf.pcens_pgamma_dunif <- function(
    object, q, pwindow, use_numeric = FALSE) {
  if (isTRUE(use_numeric)) {
    return(
      primary_censored_cdf.default(object, q, pwindow, use_numeric)
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
    stop("shape parameter is required for Gamma distribution")
  }
  if (is.null(scale)) {
    stop("scale or rate parameter is required for Gamma distribution")
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
#' @inheritParams primary_censored_cdf
#'
#' @family primary_censored_dist
#'
#' @export
primary_censored_cdf.pcens_plnorm_dunif <- function(
    object, q, pwindow, use_numeric = FALSE) {
  if (isTRUE(use_numeric)) {
    return(
      primary_censored_cdf.default(object, q, pwindow, use_numeric)
    )
  }

  # Extract Log-Normal distribution parameters
  mu <- object$args$meanlog
  sigma <- object$args$sdlog
  if (is.null(mu)) {
    stop("meanlog parameter is required for Log-Normal distribution")
  }
  if (is.null(sigma)) {
    stop("sdlog parameter is required for Log-Normal distribution")
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
#' @inheritParams primary_censored_cdf
#'
#' @family primary_censored_dist
#'
#' @export
primary_censored_cdf.pcens_pweibull_dunif <- function(
    object, q, pwindow, use_numeric = FALSE) {
  if (isTRUE(use_numeric)) {
    return(
      primary_censored_cdf.default(object, q, pwindow, use_numeric)
    )
  }

  # Extract Weibull distribution parameters
  shape <- object$args$shape
  scale <- object$args$scale
  if (is.null(shape)) {
    stop("shape parameter is required for Weibull distribution")
  }
  if (is.null(scale)) {
    stop("scale parameter is required for Weibull distribution")
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
    vapply(scaled_t, function(x) {
      if (x <= 0) {
        return(0)
      }
      pracma::gammainc(x, 1 + inv_shape)["lowinc"]
    }, numeric(1))
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
    pweibull_q <- partial_pweibull(non_zero_q)
    pweibull_q_pwindow <- partial_pweibull(non_zero_q + pwindow)
    g_q <- g(non_zero_q)
    g_q_pwindow <- g(non_zero_q + pwindow)

    Q_T <- 1 - pweibull_q_pwindow
    Delta_g <- (g_q_pwindow - g_q)
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
