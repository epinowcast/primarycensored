#' Extract Base Function Name
#'
#' This helper function extracts the base name of a function, removing an
#' namespace prefixes.
#'
#' @param func a density or distribution function
#'
#' @return A character string representing the base name of the function.
#'
#' @keywords internal
.extract_function_name <- function(func) {
  bd <- grep(".Call", deparse(body(func)), value = TRUE, fixed = TRUE)
  if (length(bd) == 1) {
    return(sub("^.*\\.Call\\(C_(\\w+),.+$", "\\1", x = bd))
  } else {
    return("unknown")
  }
}

#' @title Helper Method for Custom Distributions
#'
#' @description
#' [pprimarycensored()] and related functions can identify which distributions
#' are provided via the `pdist` and `dprimary` arguments when those are base R
#' functions (e.g. `punif`, `dexp`).
#'
#' If you need to use a non-base R implementation, but know the distribution
#' name, you can use this helper function to set it in a way that will be
#' detected.
#'
#' @param func a function, the `p`- or `d`- form of a distribution function
#'
#' @param name a string, starting with "p" or "d" indicating the underlying
#' distribution
#'
#' @return the function, with a "name" attribute added
#'
attach_distribution_name <- function(func, name) {
  attr(func, "name") <- name
  func
}

#' Extract and Combine Distribution Names
#'
#' This helper function attempts to determine distribution names and uses those
#' to establish a class name for potential analytical solutions.
#'
#' @inheritParams pprimarycensored
#'
#' @return a character string representing the combined distribution class
#'
#' @keywords internal
.format_class <- function(pdist, dprimary) {
  pdist_name <- attr(pdist, "name")
  if (is.null(pdist_name)) {
    pdist_name <- .extract_function_name(pdist)
  }
  dprim_name <- attr(dprimary, "name")
  if (is.null(dprim_name)) {
    dprim_name <- .extract_function_name(dprimary)
  }
  sprintf("pcens_%s_%s", pdist_name, dprim_name)
}

#' @title Deprecation helper
#'
#' @description
#' Provides lifecycle management consistently across several functions.
#' Currently uses [lifecycle::deprecate_soft()] - i.e. warns only when used
#' directly. In future versions, this will switch to warning unconditionally
#' with [lifecycle::deprecate_warn()], then throwing via
#' [lifecycle::deprecate_warn()], and finally be deleted along with the subject
#' arguments.
#'
#' @param pdist_name the deprecated variable to check
#' @param dprimary_name the deprecated variable to check
#' @inheritParams lifecycle::deprecate_soft
#'
#' @keywords internal
.name_deprecation <- function(
  pdist_name, dprimary_name,
  env = rlang::caller_env(),
  user_env = rlang::caller_env(2)
) {
  test_use <- c(
    lifecycle::is_present(pdist_name),
    lifecycle::is_present(dprimary_name)
  )
  res <- list(pdist = NULL, dprimary = NULL)
  if (any(test_use)) {
    lifecycle::deprecate_soft(
      when = "1.0.0",
      what = I("`pdist_name` and `dprimary_name` are deprecated across all
        functions and will be ignored in future versions;"
      ),
      details = "Use `attach_distribution_name()` instead.",
      env = env, user_env = user_env
    )
    if (test_use[1]) res$pdist <- pdist_name
    if (test_use[2]) res$dprimary <- dprimary_name
  }
  return(res)
}
