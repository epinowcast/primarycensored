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
  browser()
  my_call <- quote(substitute(func))
  func_name <- eval(my_call)
  
  for(i in rev(head(sys.frames(), 1L))) {
    my_call[[2]] <- func_name
    func_name <- eval(my_call, i)
  }
  
  func_name <- deparse(func_name)
  base_name <- sub("^.*::", "", func_name)
  return(base_name)
}
