#' Compute the primary event censored PMF for delays
#'
#'
#' This function computes the primary event censored probability mass function
#' (PMF) for a given set of quantiles. It adjusts the PMF of the primary event
#' distribution by accounting for the delay distribution and potential
#' truncation at a maximum delay (D) and minimum delay (L). The function allows
#' for custom primary event distributions and delay distributions.
#'
#' @inheritParams pprimarycensored
#'
#' @param x Vector of quantiles
#'
#' @param swindow Secondary event window (default: 1)
#'
#' @param log Logical; if TRUE, probabilities p are given as log(p)
#'
#' @param ... Additional arguments to be passed to the distribution function
#'
#' @return Vector of primary event censored PMFs, normalized over \[L, D\] if
#' truncation is applied
#'
#' @aliases dpcens
#'
#' @importFrom stats dunif
#'
#' @export
#'
#' @details
#' The primary event censored PMF is computed by taking the difference of the
#' primary event censored cumulative distribution function (CDF) at two points,
#' \eqn{d + \text{swindow}} and \eqn{d}. The primary event censored PMF,
#' \eqn{f_{\text{cens}}(d)}, is given by:
#' \deqn{
#' f_{\text{cens}}(d) = F_{\text{cens}}(d + \text{swindow}) - F_{\text{cens}}(d)
#' }
#' where \eqn{F_{\text{cens}}} is the primary event censored CDF.
#'
#' The function creates a `pcens` object with [new_pcens()] and computes
#' the PMF with [pcens_pmf()]. This evaluates the CDF once for all unique
#' points (including both \eqn{d} and \eqn{d + \text{swindow}}) and
#' reuses these values to calculate the PMF for each input value.
#'
#' When the secondary censoring interval extends past the upper truncation
#' point (\eqn{d + \text{swindow} > D}) but the lower endpoint satisfies
#' \eqn{d < D}, the upper endpoint is internally clipped to \eqn{D} before
#' evaluating the CDF. The likelihood for such an observation is
#' \eqn{P(X \in [d, \min(d + \text{swindow}, D)] \mid L \le X \le D)}, which
#' equals the usual interval probability when \eqn{d + \text{swindow} \le D}.
#' This avoids erroring when an observation's secondary window straddles the
#' truncation point (relevant for non-parametric delays such as
#' [pdiscretestep()]).
#'
#' Observations with \eqn{d \ge D} are rejected with an error: under the
#' truncation \eqn{X \le D}, no event with latent value \eqn{d \ge D} is
#' observable, and accepting such inputs would otherwise yield a 0/0
#' likelihood.
#'
#' The PMF is normalised to
#' ensure it sums to 1 over the range \[L, D\). This normalization uses:
#' \deqn{
#' f_{\text{cens,norm}}(d) = \frac{f_{\text{cens}}(d)}{
#'   F_{\text{cens}}(D) - F_{\text{cens}}(L)}
#' }
#' where \eqn{f_{\text{cens,norm}}(d)} is the normalized PMF. For the
#' explanation and mathematical details of the CDF, refer to the documentation
#' of [pprimarycensored()].
#'
#' @family primarycensored
#'
#' @importFrom stats setNames
#'
#' @examples
#' # Example: Weibull distribution with uniform primary events
#' dprimarycensored(c(0.1, 0.5, 1), pweibull, shape = 1.5, scale = 2.0)
#'
#' # Example: Weibull distribution with exponential growth primary events
#' dprimarycensored(
#'   c(0.1, 0.5, 1), pweibull,
#'   dprimary = dexpgrowth,
#'   primary_args = list(r = 0.2), shape = 1.5, scale = 2.0
#' )
#'
#' # Example: Left-truncated distribution (e.g., for generation intervals)
#' dprimarycensored(1:9, pweibull, L = 1, D = 10, shape = 1.5, scale = 2.0)
dprimarycensored <- function(
    x,
    pdist,
    pwindow = 1,
    swindow = 1,
    L = -Inf,
    D = Inf,
    dprimary = dunif,
    primary_args = NULL,
    pprimary = NULL,
    dprimary_args = NULL,
    log = FALSE,
    ...,
    check = TRUE) {
  .check_truncation_bounds(L, D)

  primary_args <- .resolve_primary_args(
    primary_args, dprimary_args, "dprimarycensored"
  )
  pcens_obj <- .build_pcens(
    pdist, dprimary, primary_args, pprimary, list(...),
    pwindow = pwindow, D = D, check = check
  )

  pcens_pmf(
    pcens_obj, x, pwindow,
    swindow = swindow, L = L, D = D, log = log
  )
}

#' @rdname dprimarycensored
#' @export
dpcens <- dprimarycensored
