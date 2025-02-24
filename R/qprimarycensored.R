#' Compute quantiles corresponding to target probabilities for primary event
#' censored delays
#'
#' This function computes the quantiles (delay values) that correspond to
#' specified probabilities in the primary event censored distribution. For a
#' given probability p, it computes the delay value q such that P(X ≤ q) = p,
#' where X follows the primary event censored distribution. The distribution
#' accounts for both the delay distribution and the primary event timing
#' distribution.
#'
#' @param p Vector of probabilities between 0 and 1 for which to compute
#'  corresponding quantiles
#'
#' @inheritParams pprimarycensored
#'
#' @return Vector of delay values (quantiles) corresponding to the input
#'  probabilities
#'
#' @aliases qpcens
#'
#' @importFrom stats dunif
#'
#' @export
#'
#' @details
#' For each probability p, the function computes the delay value q such that
#' P(X ≤ q) = p, where X follows the primary event censored distribution.
#' This is done by inverting the primary event censored CDF.
#'
#' The function creates a `primarycensored` object using [new_pcens()] and then
#' computes the quantiles using [pcens_quantile()]. This approach allows for
#' analytical solutions when available, falling back to numerical methods when
#' necessary.
#'
#' For example, if p = 0.5, the function returns the median delay - the value
#' where 50% of censored events occur by this time and 50% occur after.
#'
#' See `methods(pcens_quantile)` for which combinations have analytical
#' solutions implemented.
#'
#' @family primarycensored
#' @seealso [new_pcens()] and [pcens_quantile()]
#'
#' @examples
#' # Compute delays where 25%, 50%, and 75% of events occur by (quartiles)
#' # Using lognormal delays with uniform primary events
#' qprimarycensored(c(0.25, 0.5, 0.75), plnorm, meanlog = 0, sdlog = 1)
#'
#' # Same quartiles but with exponential growth in primary events
#' qprimarycensored(
#'   c(0.25, 0.5, 0.75), plnorm,
#'   dprimary = dexpgrowth,
#'   dprimary_args = list(r = 0.2), meanlog = 0, sdlog = 1
#' )
#'
#' # Same quartiles but with truncation at 10
#' qprimarycensored(
#'   c(0.25, 0.5, 0.75), plnorm,
#'   dprimary = dexpgrowth,
#'   dprimary_args = list(r = 0.2), meanlog = 0, sdlog = 1, D = 10
#' )
qprimarycensored <- function(
    p,
    pdist,
    pwindow = 1,
    D = Inf,
    dprimary = stats::dunif,
    dprimary_args = list(),
    ...) {
  check_pdist(pdist, Inf, ...)
  check_dprimary(dprimary, pwindow, dprimary_args)

  # Create a new primarycensored object
  pcens_obj <- new_pcens(
    pdist,
    dprimary,
    dprimary_args,
    ...
  )

  # Compute the quantiles using the S3 method
  pcens_quantile(pcens_obj, p, pwindow, D)
}

#' @rdname qprimarycensored
#' @export
qpcens <- qprimarycensored
