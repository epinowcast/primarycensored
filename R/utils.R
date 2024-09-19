#' Extract Base Function Name
#'
#' This helper function extracts the base name of a function, removing an
#' namespace prefixes.
#'
#' @param func The output of `substitute` on a function.
#'
#' @return A character string representing the base name of the function.
#'
#' @examples
#' .extract_function_name(substitute(plnorm))
#' .extract_function_name(substitute(stats::dunif))
#'
#' @keywords internal
.extract_function_name <- function(func) {
  func_name <- deparse(func)
  base_name <- sub("^.*::", "", func_name)
  return(base_name)
}
