#' Check if a function is a valid cumulative distribution function (CDF)
#'
#' This function tests whether a given function behaves like a valid CDF by
#' checking if it's monotonically increasing and bounded between 0 and 1.
#'
#' @inheritParams pprimarycensoreddist
#' @return NULL. The function will stop execution with an error message if
#'         pdist is not a valid CDF.
#' @export
#' @examples
#' check_pdist(pnorm, D = 10)
check_pdist <- function(pdist, D, ...) {
  if (is.infinite(D)) {
    D <- 1000
  }
  test_values <- sort(sample(seq(0, D, by = 1), 4))
  test_results <- pdist(test_values, ...)

  if (!all(diff(test_results) >= 0) ||
    !all(test_results >= 0 & test_results <= 1)) {
    stop(
      "pdist is not a valid cumulative distribution function (CDF). ",
      "Please ensure you're using a p-function (e.g., pnorm, punif) and not ",
      "a d-function (e.g., dnorm, dbinom).",
      "For values ", toString(test_values),
      " the results were ", toString(round(test_results, 3)), ". ",
      "You can use the `check_pdist` function to check if your p-function ",
      "is correct."
    )
  }
  return(invisible(NULL))
}

#' Check if a function is a valid probability density function (PDF)
#'
#' This function tests whether a given function behaves like a valid PDF by
#' checking if it integrates to approximately 1 over the specified range.
#'
#' @inheritParams pprimarycensoreddist
#' @param tolerance The tolerance for the integral to be considered close to 1
#'
#' @return NULL. The function will stop execution with an error message if
#'         dprimary is not a valid PDF.
#' @export
#'
#' @examples
#' check_dprimary(dunif, pwindow = 1)
check_dprimary <- function(dprimary, pwindow, dprimary_args = list(),
                           tolerance = 1e-3) {
  integrand <- function(x) {
    do.call(dprimary, c(list(x = x, min = 0, max = pwindow), dprimary_args))
  }
  integral <- stats::integrate(integrand, lower = 0, upper = pwindow)$value

  if (abs(integral - 1) > tolerance) {
    stop(
      "dprimary is not a valid probability density function (PDF). ",
      "It should integrate to approximately 1 over the range [0, pwindow]. ",
      "Calculated integral: ", round(integral, 4),
      "You can use the `check_dprimary` function to check if your d-function ",
      "is correct."
    )
  }
  return(invisible(NULL))
}
