#' S3 class for primary event censored distribution computation
#'
#' @inheritParams pprimarycensored
#'
#' @return An object of class `pcens_{pdist_name}_{dprimary_name}`. This
#'  contains the primary event distribution, the delay distribution, the
#'  delay distribution arguments, and any additional arguments. It can be
#'  used with the `pcens_cdf()` function to compute the primary event censored
#'  cdf.
#'
#' @family pcens
#'
#' @export
#' @examples
#' new_pcens(
#'   pdist = pgamma, dprimary = dunif, dprimary_args = list(min = 0, max = 1),
#'   shape = 1, scale = 1
#' )
new_pcens <- function(
    pdist,
    dprimary,
    dprimary_args,
    ...) {
  structure(
    list(
      pdist = pdist,
      dprimary = dprimary,
      dprimary_args = dprimary_args,
      args = list(...)
    ),
    class = .format_class(pdist, dprimary)
  )
}
