#' Extract base function name
#'
#' This helper function extracts the base name of a function, removing
#' namespace prefixes.
#'
#' @inheritParams add_name_attribute
#'
#' @return Character string representing the base name of the function.
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

#' Helper method for custom distributions
#'
#' [pprimarycensored()] and related functions can identify which distributions
#' are provided via the `pdist` and `dprimary` arguments when those are base R
#' functions (e.g. `punif`, `dexp`) via the `name` attribute.
#'
#' If you need to use a non-base R implementation, but know the distribution
#' name, you can use this helper function to set it in a way that will be
#' detected by [pprimarycensored()] and related functions.
#'
#' This is useful as it enables the automatic use of analytical solutions for
#' distributions where they exist.
#'
#' @param func Function, for example the `p`- or `d`- form of a distribution
#' function.
#'
#' @param name Character string, starting with "p" or "d" indicating the
#' underlying distribution.
#'
#' @return Function, with a "name" attribute added
#' @family utils
#' @export
#' @examples
#' dist <- add_name_attribute(pnorm, "hello")
#' attr(dist, "name")
add_name_attribute <- function(func, name) {
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

#' Deprecation name helper
#'
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
      details = "Use `add_name_attribute()` instead.",
      env = env, user_env = user_env
    )
    if (test_use[1]) res$pdist <- pdist_name
    if (test_use[2]) res$dprimary <- dprimary_name
  }
  return(res)
}
