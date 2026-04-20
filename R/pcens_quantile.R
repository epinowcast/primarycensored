#' Compute primary event censored quantiles
#'
#' This function inverts the primary event censored CDF to compute quantiles.
#' It uses numerical optimisation via optim to find the value q such that
#' [pcens_cdf()] is close to the specified probability. Currently, only the
#' default numerical inversion method is implemented. Future analytical
#' solutions may be added.
#'
#' @inheritParams pprimarycensored
#'
#' @param object A `primarycensored` object as created by [new_pcens()].
#'
#' @param p A vector of probabilities at which to compute the quantiles.
#'
#' @param use_numeric Logical; if TRUE forces the use of numeric inversion even
#' if an analytical solution is available (not yet implemented).
#'
#' @return Vector of primary event censored quantiles.
#'
#' @family pcens
#'
#' @export
pcens_quantile <- function(
    object,
    p,
    pwindow,
    L = -Inf,
    D = Inf,
    use_numeric = FALSE,
    ...) {
  UseMethod("pcens_quantile")
}

#' Default method for computing primary event censored quantiles
#'
#' This method inverts the primary event censored CDF by root-finding via
#' [stats::uniroot()]. The censored CDF is monotone in `q`, so a bracketed
#' root-finder is both simpler and more robust than the minimisation approach
#' used previously, and it works without special casing for `L = -Inf` or
#' `D = Inf` because [stats::uniroot()] can extend its search interval
#' automatically.
#'
#' @param init Half-width of the initial search interval used when one or
#'   both truncation bounds are infinite. The starting interval is taken as
#'   `[-init, init]` (or `[L, L + 2 * init]` / `[D - 2 * init, D]` when only
#'   one bound is finite) and then extended outward by [stats::uniroot()]
#'   as needed. Defaults to 5, which brackets the bulk of most commonly used
#'   delay distributions.
#'
#' @param tol Numeric tolerance passed to [stats::uniroot()].
#'
#' @param max_iter Maximum number of [stats::uniroot()] iterations.
#'
#' @param ... Additional arguments passed to underlying functions.
#'
#' @inheritParams pcens_quantile
#'
#' @family pcens
#'
#' @return A numeric vector containing the computed primary event censored
#'   quantiles.
#'
#' @export
#' @examples
#' # Create a primarycensored object with gamma delay and uniform primary
#' pcens_obj <- new_pcens(
#'   pdist = pgamma,
#'   dprimary = dunif,
#'   dprimary_args = list(min = 0, max = 1),
#'   shape = 3,
#'   scale = 2
#' )
#'
#' # Compute quantile for a single probability
#' pcens_quantile(pcens_obj, p = 0.8, pwindow = 1)
#'
#' # Compute quantiles for multiple probabilities
#' pcens_quantile(pcens_obj, p = c(0.25, 0.5, 0.75), pwindow = 1)
#'
#' # Compute quantiles for multiple probabilities with truncation
#' pcens_quantile(pcens_obj, p = c(0.25, 0.5, 0.75), pwindow = 1, D = 10)
#'
#' # Compute quantiles with left truncation
#' pcens_quantile(pcens_obj, p = c(0.25, 0.5, 0.75), pwindow = 1, L = 1, D = 10)
pcens_quantile.default <- function(
    object,
    p,
    pwindow,
    L = -Inf,
    D = Inf,
    use_numeric = FALSE,
    init = 5,
    tol = 1e-8,
    max_iter = 10000,
    ...) {
  .check_truncation_bounds(L, D)

  # Build an initial bracket for `stats::uniroot`. When a bound is infinite
  # we use a finite stand-in and rely on `extendInt = "yes"` to expand
  # outward until the censored CDF brackets `prob`. This avoids any need for
  # a "starting point" in the flat tail of the CDF (the failure mode that
  # would previously strand L-BFGS-B at a zero gradient).
  if (is.finite(L) && is.finite(D)) {
    lower <- L
    upper <- D
  } else if (is.finite(L)) {
    lower <- L
    upper <- L + 2 * init
  } else if (is.finite(D)) {
    lower <- D - 2 * init
    upper <- D
  } else {
    lower <- -init
    upper <- init
  }

  sapply(p, function(prob) {
    if (prob <= 0) {
      return(L)
    }
    if (prob >= 1) {
      return(NA_real_)
    }

    objective <- function(q) {
      cdf_val <- pcens_cdf(object, q, pwindow, use_numeric)
      cdf_val <- .normalise_cdf(cdf_val, q, L, D, object, pwindow)
      cdf_val - prob
    }

    stats::uniroot(
      objective,
      lower = lower,
      upper = upper,
      extendInt = "upX",
      tol = tol,
      maxiter = max_iter
    )$root
  })
}
