#' Extract Base Function Name
#'
#' This helper function extracts the base name of a function, removing an
#' namespace prefixes.
#'
#' @param func The output of `substitute` on a function.
#'
#' @return A character string representing the base name of the function.
#'
#' @keywords internal
.extract_function_name <- function(func) {
  bd <- grep(".Call", deparse(body(func)), value = TRUE, fixed = TRUE)
  return(sub("^.*\\.Call\\(C_(\\w+),.+$", "\\1", x = bd))
}
