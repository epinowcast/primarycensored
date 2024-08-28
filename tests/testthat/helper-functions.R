on_ci <- function() {
  isTRUE(as.logical(Sys.getenv("CI")))
}

not_on_cran <- function() {
  identical(Sys.getenv("NOT_CRAN"), "true")
}

skip_on_local <- function() {
  if (on_ci()) {
    return(invisible(TRUE))
  }
  testthat::skip("Not on CI")
}
