#' Check if a function is a valid cumulative distribution function (CDF)
#'
#' This function tests whether a given function behaves like a valid CDF by
#' checking if it's monotonically increasing and bounded between 0 and 1.
#'
#' @inheritParams pprimarycensored
#' @return NULL. The function will stop execution with an error message if
#'         pdist is not a valid CDF.
#' @export
#'
#' @family check
#'
#' @examples
#' check_pdist(pnorm, D = 10)
check_pdist <- function(pdist, D, ...) {
  if (is.infinite(D)) {
    D <- 1000
  }
  test_values <- sort(runif(4, 0, D))
  test_results <- pdist(test_values, ...)

  if (
    !all(diff(test_results) >= 0) ||
      !all(test_results >= 0 & test_results <= 1)
  ) {
    stop(
      "pdist is not a valid cumulative distribution function (CDF). ",
      "Please ensure you're using a p-function (e.g., pnorm, punif) and not ",
      "a d-function (e.g., dnorm, dbinom).",
      "For values ",
      toString(test_values),
      " the results were ",
      toString(round(test_results, 3)),
      ". ",
      "You can use the `check_pdist` function to check if your p-function ",
      "is correct.",
      call. = FALSE
    )
  }
  return(invisible(NULL))
}

#' Check if a function is a valid bounded probability density function (PDF)
#'
#' This function tests whether a given function behaves like a valid PDF by
#' checking if it integrates to approximately 1 over the specified range
#' and if it takes the arguments min and max.
#'
#' @inheritParams pprimarycensored
#' @param tolerance The tolerance for the integral to be considered close to 1
#'
#' @return NULL. The function will stop execution with an error message if
#'         dprimary is not a valid PDF.
#' @export
#'
#' @family check
#'
#' @examples
#' check_dprimary(dunif, pwindow = 1)
check_dprimary <- function(
    dprimary,
    pwindow,
    dprimary_args = list(),
    tolerance = 1e-3) {
  # check if dprimary takes min and max as arguments
  if (!all(c("min", "max") %in% names(formals(dprimary)))) {
    stop("dprimary must take min and max as arguments", call. = FALSE)
  }

  integrand <- function(x) {
    do.call(dprimary, c(list(x = x, min = 0, max = pwindow), dprimary_args))
  }
  integral <- stats::integrate(integrand, lower = 0, upper = pwindow)$value

  if (abs(integral - 1) > tolerance) {
    stop(
      "dprimary is not a valid probability density function (PDF). ",
      "It should integrate to approximately 1 over the range [0, pwindow]. ",
      "Calculated integral: ",
      round(integral, 4),
      "You can use the `check_dprimary` function to check if your d-function ",
      "is correct.",
      call. = FALSE
    )
  }
  return(invisible(NULL))
}

#' Validate truncation bounds L and D
#'
#' Internal function to validate that L (lower truncation) and D (upper
#' truncation) parameters are valid: L must be non-negative and less than D.
#'
#' @param L Lower truncation bound
#' @param D Upper truncation bound
#'
#' @return Invisible NULL if valid, otherwise stops with an error message.
#'
#' @keywords internal
.check_truncation_bounds <- function(L, D) {
  if (L < 0) {
    stop("L must be non-negative.", call. = FALSE)
  }
  if (L >= D) {
    stop("L must be less than D.", call. = FALSE)
  }
  invisible(NULL)
}

#' Check if truncation time is appropriate relative to the maximum delay
#'
#' This function checks if the truncation time D is appropriate relative to the
#' maximum delay. If D is much larger than necessary, it suggests
#' considering setting it to `Inf` for better efficiency with minimal accuracy
#' cost.
#'
#' @param delays A numeric vector of delay times
#'
#' @param D The truncation time
#'
#' @param multiplier The multiplier for the maximum delay to compare with D.
#'   Default is 2.
#'
#' @return Invisible NULL. Prints a message if the condition is met.
#'
#' @export
#' @family check
#'
#' @examples
#' check_truncation(delays = c(1, 2, 3, 4), D = 10, multiplier = 2)
check_truncation <- function(delays, D, multiplier = 2) {
  if (!is.numeric(delays) || !is.numeric(D) || !is.numeric(multiplier)) {
    stop("All arguments must be numeric.", call. = FALSE)
  }

  if (length(D) > 1) {
    stop("D must be a single value for check_truncation", call. = FALSE)
  }

  if (D <= 0 || multiplier <= 1) {
    stop(
      "Invalid argument values. D must be positive and multiplier must be ",
      "greater than 1.",
      call. = FALSE
    )
  }

  if (is.infinite(D)) {
    return(invisible(NULL))
  }

  # Remove NA
  delays <- delays[!is.na(delays)]

  if (length(delays) == 0) {
    warning("No finite observed delays to check.", call. = FALSE)
    return(invisible(NULL))
  }

  max_delay <- max(delays)
  expected_D <- max_delay * multiplier

  # Check if D is much larger than expected
  if (D > expected_D) {
    message(
      sprintf(
        paste0(
          "The truncation time D (%g) is larger than %g times the maximum ",
          "observed delay (%g). Consider setting D to Inf for better ",
          "efficiency with minimal accuracy cost for this case."
        ),
        D,
        multiplier,
        max_delay
      )
    )
  }
  invisible(NULL)
}
