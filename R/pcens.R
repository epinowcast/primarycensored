#' S3 class for primary event censored distribution computation
#'
#' @inheritParams pprimarycensored
#'
#' @param pprimary CDF of the primary event distribution. May be a function
#'   or a character string naming a primary distribution in
#'   \code{pcd_primary_distributions}. When \code{NULL} (the default), it is
#'   looked up automatically from the registry using the \code{"name"}
#'   attribute of \code{dprimary}. When both \code{dprimary} and
#'   \code{pprimary} carry a name, the two must agree on everything other
#'   than the leading \code{d}/\code{p} prefix; mismatches such as
#'   \code{dunif} + \code{pexpgrowth} raise an error.
#'
#' @param primary_args List of additional arguments to be passed to
#'   \code{dprimary} (and the looked-up \code{pprimary}). Replaces the
#'   deprecated \code{dprimary_args}.
#'
#' @param dprimary_args \[Deprecated\] Use \code{primary_args} instead.
#'
#' @return An object with class hierarchy
#'   \code{c("pcens_{pdist_name}_{dprimary_name}", "pcens_{pdist_name}",
#'   "pcens")}. This contains the primary event distribution, the delay
#'   distribution, the delay distribution arguments, the primary event CDF
#'   (if available), and any additional arguments. It can be used with the
#'   \code{pcens_cdf()} function to compute the primary event censored CDF.
#'
#' @family pcens
#'
#' @export
#' @examples
#' new_pcens(
#'   pdist = pgamma, dprimary = dunif,
#'   primary_args = list(min = 0, max = 1),
#'   shape = 1, scale = 1
#' )
new_pcens <- function(
    pdist,
    dprimary,
    primary_args = NULL,
    pprimary = NULL,
    dprimary_args = NULL,
    ...) {
  primary_args <- .resolve_primary_args(
    primary_args, dprimary_args, "new_pcens"
  )
  pprimary <- .resolve_pprimary(
    dprimary, pprimary
  )
  obj <- list(
    pdist = pdist,
    dprimary = dprimary,
    primary_args = primary_args,
    # Keep dprimary_args alias for backward compatibility with downstream
    # consumers that read object$dprimary_args directly.
    dprimary_args = primary_args,
    pprimary = pprimary,
    args = list(...)
  )
  class(obj) <- .format_class(pdist, dprimary)
  obj
}
