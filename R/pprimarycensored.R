#' Compute the primary event censored CDF for delays
#'
#' This function computes the primary event censored cumulative distribution
#' function (CDF) for a given set of quantiles. It adjusts the CDF of the
#' primary event distribution by accounting for the delay distribution and
#' potential truncation at a maximum delay (D). The function allows for
#' custom primary event distributions and delay distributions.
#'
#' @param q Vector of quantiles
#'
#' @param pdist Distribution function (CDF). The package can identify base R
#'  distributions for potential analytical solutions. For non-base R functions,
#'  users can apply [add_name_attribute()] to yield properly tagged
#'  functions if they wish to leverage the analytical solutions.
#'
#' @param pwindow Primary event window
#'
#' @param D Maximum delay (truncation point). If finite, the distribution is
#'  truncated at D. If set to Inf, no truncation is applied. Defaults to Inf.
#'
#' @param dprimary Function to generate the probability density function
#'  (PDF) of primary event times. This function should take a value `x` and a
#'  `pwindow` parameter, and return a probability density. It should be
#'  normalized to integrate to 1 over \[0, pwindow\]. Defaults to a uniform
#'  distribution over \[0, pwindow\]. Users can provide custom functions or use
#'  helper functions like `dexpgrowth` for an exponential growth distribution.
#'  See [pcd_primary_distributions()] for examples. The package can identify
#'  base R distributions for potential analytical solutions. For non-base R
#'  functions, users can apply [add_name_attribute()] to yield properly tagged
#'  functions if they wish to leverage analytical solutions.
#'
#' @param dprimary_args List of additional arguments to be passed to
#'  dprimary. For example, when using `dexpgrowth`, you would
#'  pass `list(min = 0, max = pwindow, r = 0.2)` to set the minimum, maximum,
#'  and rate parameters
#'
#' @param pdist_name `r lifecycle::badge("deprecated")` this argument will be
#'  ignored in future versions; use [add_name_attribute()] on `pdist`
#'  instead
#'
#' @param dprimary_name `r lifecycle::badge("deprecated")` this argument will be
#'  ignored in future versions; use [add_name_attribute()] on `dprimary`
#'  instead
#'
#' @param ... Additional arguments to be passed to pdist
#'
#' @return Vector of primary event censored CDFs, normalized by D if finite
#'  (truncation adjustment)
#'
#' @aliases ppcens
#'
#' @importFrom stats dunif
#'
#' @export
#'
#' @details
#' The primary event censored CDF is computed by integrating the product of
#' the delay distribution function (CDF) and the primary event distribution
#' function (PDF) over the primary event window. The integration is adjusted
#' for truncation if a finite maximum delay (D) is specified.
#'
#' The primary event censored CDF, \eqn{F_{\text{cens}}(q)}, is given by:
#' \deqn{
#' F_{\text{cens}}(q) = \int_{0}^{pwindow} F(q - p) \cdot f_{\text{primary}}(p)
#' \, dp
#' }
#' where \eqn{F} is the CDF of the delay distribution,
#' \eqn{f_{\text{primary}}} is the PDF of the primary event times, and
#' \eqn{pwindow} is the primary event window.
#'
#' If the maximum delay \eqn{D} is finite, the CDF is normalized by dividing
#' by \eqn{F_{\text{cens}}(D)}:
#' \deqn{
#' F_{\text{cens,norm}}(q) = \frac{F_{\text{cens}}(q)}{F_{\text{cens}}(D)}
#' }
#' where \eqn{F_{\text{cens,norm}}(q)} is the normalized CDF.
#'
#' This function creates a `primarycensored` object using
#' [new_pcens()] and then computes the primary event
#' censored CDF using [pcens_cdf()]. This abstraction allows
#' for automatic use of analytical solutions when available, while
#' seamlessly falling back to numerical integration when necessary.
#'
#' See `methods(pcens_cdf)` for which combinations have analytical
#' solutions implemented.
#'
#' @family primarycensored
#' @seealso [new_pcens()] and [pcens_cdf()]
#'
#' @examples
#' # Example: Lognormal distribution with uniform primary events
#' pprimarycensored(c(0.1, 0.5, 1), plnorm, meanlog = 0, sdlog = 1)
#'
#' # Example: Lognormal distribution with exponential growth primary events
#' pprimarycensored(
#'   c(0.1, 0.5, 1), plnorm,
#'   dprimary = dexpgrowth,
#'   dprimary_args = list(r = 0.2), meanlog = 0, sdlog = 1
#' )
pprimarycensored <- function(
    q,
    pdist,
    pwindow = 1,
    D = Inf,
    dprimary = stats::dunif,
    dprimary_args = list(),
    pdist_name = lifecycle::deprecated(),
    dprimary_name = lifecycle::deprecated(),
    ...) {
  nms <- .name_deprecation(pdist_name, dprimary_name)
  if (!is.null(nms$pdist)) {
    pdist <- add_name_attribute(pdist, nms$pdist)
  }
  if (!is.null(nms$dprimary)) {
    dprimary <- add_name_attribute(dprimary, nms$dprimary)
  }

  check_pdist(pdist, D, ...)
  check_dprimary(dprimary, pwindow, dprimary_args)

  # Create a new primarycensored object
  pcens_obj <- new_pcens(
    pdist,
    dprimary,
    dprimary_args,
    ...
  )

  # Compute the CDF using the S3 method
  result <- pcens_cdf(pcens_obj, q, pwindow)

  if (!is.infinite(D)) {
    result <- .normalise_cdf(result, q, D, pcens_obj, pwindow)
  }

  return(result)
}

#' Normalise a primary event censored CDF
#'
#'
#' Internal function to normalise a primary event censored CDF when truncation
#' is applied. The CDF is normalised by dividing by its value at the truncation
#' point D and setting all values beyond D to 1.
#'
#' @param result Numeric vector of CDF values to normalise.
#'
#' @param q Numeric vector of quantiles at which CDF was evaluated.
#'
#' @param D Numeric truncation point
#'
#' @param pcens_obj A primarycensored object as created by [new_pcens()].
#'
#' @param pwindow Secondary event window
#'
#' @return Normalised CDF values as a numeric vector
#'
#' @keywords internal
.normalise_cdf <- function(result, q, D, pcens_obj, pwindow) {
  if (max(q) == D) {
    adj <- result[which.max(q == D)]
  } else {
    adj <- pcens_cdf(pcens_obj, D, pwindow)
  }
  result <- result / adj
  result <- ifelse(q > D, 1, result)
  return(result)
}

#' @rdname pprimarycensored
#' @export
ppcens <- pprimarycensored
