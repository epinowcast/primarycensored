#' S3 class for primary event censored distribution computation
#'
#' @inheritParams pprimarycensoreddist
#'
#' @return An object of class primary_censored_cdf
#'
#' @family primary_censored_dist
#'
#' @export
new_primary_censored_dist <- function(pdist, dprimary, dprimary_args, ...) {
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
        deparse(substitute(pdist)), "_",
        deparse(substitute(dprimary))
      ),
      "pcens_numeric"
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

#' Numeric method for computing primary event censored CDF
#'
#' This method uses numerical integration to compute the primary event censored
#' CDF for any combination of delay distribution and primary event distribution.
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
#' @export
primary_censored_cdf.pcens_numeric <- function(
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
  use_numeric <- TRUE
  if (isTRUE(use_numeric)) {
    NextMethod()
  }

  result <- vapply(q, function(n) {
    # Implement analytical solution here
  }, numeric(1))

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
  use_numeric <- TRUE
  if (isTRUE(use_numeric)) {
    NextMethod()
  }

  result <- vapply(q, function(n) {
    # Implement analytical solution here
  }, numeric(1))

  return(result)
}
