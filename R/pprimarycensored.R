#' Compute the primary event censored CDF for delays
#'
#' This function computes the primary event censored cumulative distribution
#' function (CDF) for a given set of quantiles. It adjusts the CDF of the
#' primary event distribution by accounting for the delay distribution and
#' potential truncation at a maximum delay (D) and minimum delay (L). The
#' function allows for custom primary event distributions and delay
#' distributions.
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
#' @param D Maximum delay (upper truncation point). If finite, the distribution
#'  is truncated at D. If set to Inf, no upper truncation is applied. Defaults
#'  to Inf.
#'
#' @param L Minimum delay (lower truncation point). If greater than 0, the
#'  distribution is left-truncated at L. This is useful for modelling
#'  generation intervals where day 0 is excluded. Defaults to 0 (no left
#'  truncation).
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
#' @param ... Additional arguments to be passed to pdist
#'
#' @return Vector of primary event censored CDFs, normalized over \[L, D\] if
#'  truncation is applied
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
#' for truncation if specified.
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
#' If truncation is applied (finite D or L > 0), the CDF is normalized:
#' \deqn{
#' F_{\text{cens,norm}}(q) = \frac{F_{\text{cens}}(q) - F_{\text{cens}}(L)}{
#' F_{\text{cens}}(D) - F_{\text{cens}}(L)}
#' }
#' where \eqn{F_{\text{cens,norm}}(q)} is the normalized CDF. For values
#' \eqn{q \leq L}, the function returns 0; for values \eqn{q \geq D}, it
#' returns 1.
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
#'
#' # Example: Left-truncated distribution (e.g., for generation intervals)
#' pprimarycensored(
#'   c(1, 2, 3), plnorm,
#'   D = 10, L = 1,
#'   meanlog = 0, sdlog = 1
#' )
pprimarycensored <- function(
    q,
    pdist,
    pwindow = 1,
    D = Inf,
    L = 0,
    dprimary = stats::dunif,
    dprimary_args = list(),
    ...) {
  .check_truncation_bounds(L, D)

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

  # Apply truncation normalization if needed
  if (!is.infinite(D) || L > 0) {
    result <- .normalise_cdf(result, q, D, L, pcens_obj, pwindow)
  }

  return(result)
}

#' Normalise a primary event censored CDF
#'
#' Internal function to normalise a primary event censored CDF when truncation
#' is applied. The CDF is normalised using (F(q) - F(L)) / (F(D) - F(L)) and
#' values outside \[L, D\] are clamped to 0 or 1.
#'
#' @param result Numeric vector of CDF values to normalise.
#'
#' @param q Numeric vector of quantiles at which CDF was evaluated.
#'
#' @param D Numeric upper truncation point
#'
#' @param L Numeric lower truncation point
#'
#' @param pcens_obj A primarycensored object as created by [new_pcens()].
#'
#' @param pwindow Secondary event window
#'
#' @return Normalised CDF values as a numeric vector
#'
#' @keywords internal
.normalise_cdf <- function(result, q, D, L, pcens_obj, pwindow) {
  # Get CDF at upper truncation point D
  if (any(q == D)) {
    cdf_D <- result[which.max(q == D)]
  } else if (is.infinite(D)) {
    cdf_D <- 1
  } else {
    cdf_D <- pcens_cdf(pcens_obj, D, pwindow)
  }

  # Get CDF at lower truncation point L
  if (L == 0) {
    cdf_L <- 0
  } else if (any(q == L)) {
    cdf_L <- result[which.max(q == L)]
  } else {
    cdf_L <- pcens_cdf(pcens_obj, L, pwindow)
  }

  # Normalise: (F(q) - F(L)) / (F(D) - F(L))
  normaliser <- cdf_D - cdf_L
  result <- (result - cdf_L) / normaliser

  # Clamp values outside truncation range
  result <- ifelse(q <= L, 0, result)
  result <- ifelse(q >= D, 1, result)

  return(result)
}

#' @rdname pprimarycensored
#' @export
ppcens <- pprimarycensored
