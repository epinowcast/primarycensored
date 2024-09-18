#' Internal S3 class for primary event censored CDF computation
#'
#' @inheritParams pprimarycensoreddist
#'
#' @return An object of class primary_censored_cdf
#'
#' @importFrom stats pgamma plnorm dunif integrate
#'
#' @export
new_primary_censored_cdf <- function(pdist, dprimary, dprimary_args, ...) {
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
        deparse(substitute(dprimary)),
        "_cdf"
      ),
      "pcens_numeric_cdf"
    )
  )
}

#' Compute primary event censored CDF
#'
#' @inheritParams pprimarycensoreddist
#'
#' @param object A primary_censored_cdf object
#'
#' @param use_numeric Logical, if TRUE forces use of numeric integration
#' even for distributions with analytical solutions.
#'
#' @return Vector of primary event censored CDFs
#'
#' @export
compute_primary_censored_cdf <- function(
    object, q, pwindow, use_numeric = FALSE) {
  UseMethod("compute_primary_censored_cdf")
}

#' Numeric method for computing primary event censored CDF
#'
#' @export
compute_primary_censored_cdf.pcens_numeric_cdf <- function(
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
#' @inheritParams compute_primary_censored_cdf
#'
#' @export
compute_primary_censored_cdf.pcens_pgamma_dunif <- function(
    object, q, pwindow, use_numeric = FALSE) {
  if (use_numeric) {
    NextMethod()
  }

  k <- object$args$shape
  theta <- object$args$scale

  result <- vapply(q, function(n) {
    # Implement analytical solution here
  }, numeric(1))

  return(result)
}

#' Method for Log-Normal delay with uniform primary
#'
#' @inheritParams compute_primary_censored_cdf
#'
#' @export
compute_primary_censored_cdf.pcens_plnorm_dunif <- function(
    object, q, pwindow, use_numeric = FALSE) {
  if (use_numeric) {
    NextMethod()
  }

  mu <- object$args$meanlog
  sigma <- object$args$sdlog

  result <- vapply(q, function(n) {
    # Implement analytical solution here
  }, numeric(1))

  return(result)
}
