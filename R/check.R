#' Check if a function is a valid cumulative distribution function (CDF)
#'
#' This function tests whether a given function behaves like a valid CDF by
#' checking if it's monotonically increasing and bounded between 0 and 1.
#'
#' @param dist_func The distribution function to be checked
#'
#' @param D The maximum delay (truncation point)
#'
#' @param swindow The secondary event window
#'
#' @param ... Additional arguments to be passed to dist_func
#'
#' @return NULL. The function will stop execution with an error message if
#'         dist_func is not a valid CDF.
check_dist_func <- function(dist_func, D, swindow, ...) {
  test_values <- sample(seq(0, D, by = swindow), 3)
  test_results <- dist_func(test_values, ...)

  if (!all(diff(test_results) >= 0) ||
    !all(test_results >= 0 & test_results <= 1)) {
    stop(
      "dist_func is not a valid cumulative distribution function (CDF). ",
      "Please ensure you're using a p-function (e.g., pnorm, punif) and not ",
      "a d-function (e.g., dnorm, dbinom)."
    )
  }
  return(invisible(NULL))
}
