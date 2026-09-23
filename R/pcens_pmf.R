#' Compute the primary event censored PMF from a pcens object
#'
#' This is the PMF counterpart of [pcens_cdf()]. It computes the primary
#' event censored PMF from a `pcens` object as created by [new_pcens()] or
#' [update.pcens()]. It handles secondary event windows and truncation in the
#' same way as [dprimarycensored()], so callers do not need to difference
#' CDFs themselves. Unlike [dprimarycensored()], it does not look up
#' distributions by name or validate `pdist` and `dprimary`, which makes it
#' cheaper when one distribution is evaluated many times.
#'
#' @inheritParams pcens_cdf
#' @inheritParams dprimarycensored
#'
#' @param object A `pcens` object as created by [new_pcens()].
#'
#' @details
#' The PMF at `x` is the difference of the primary event censored CDF at
#' `min(x + swindow, D)` and at `x`, normalised by
#' \eqn{F_{\text{cens}}(D) - F_{\text{cens}}(L)}. Values of `x` below `L`
#' and, for finite `D`, values of `x` at or above `D` raise an error. See
#' [dprimarycensored()] for the details.
#'
#' @inherit dprimarycensored return
#'
#' @family pcens
#' @seealso [dprimarycensored()], [pcens_cdf()] and [update.pcens()]
#'
#' @export
#' @examples
#' obj <- new_pcens(
#'   pdist = pgamma, dprimary = dunif,
#'   primary_args = list(min = 0, max = 1),
#'   shape = 3, scale = 2
#' )
#' pcens_pmf(obj, x = 0:9, pwindow = 1, D = 10)
#'
#' # Evaluate the same distribution for a new parameter set
#' pcens_pmf(update(obj, shape = 2), x = 0:9, pwindow = 1, D = 10)
pcens_pmf <- function(
    object,
    x,
    pwindow,
    swindow = 1,
    L = -Inf,
    D = Inf,
    log = FALSE) {
  if (!inherits(object, "pcens")) {
    stop(
      "object must be a pcens object as created by new_pcens().",
      call. = FALSE
    )
  }
  .check_truncation_bounds(L, D)

  if (min(x) < L) {
    stop(
      "Some values of x are below L. Minimum x is ",
      min(x),
      " and L is ",
      L,
      ". Resolve this by filtering x to only include values >= L.",
      call. = FALSE
    )
  }

  if (is.finite(D) && max(x) >= D) {
    stop(
      "Upper truncation point is greater than D. Maximum x is ",
      max(x),
      " and D is ",
      D,
      ". Under truncation at D no event with latent value >= D is ",
      "observable; resolve this by filtering x to values strictly less than D.",
      call. = FALSE
    )
  }

  # Clip the upper end of each secondary interval at D so observations with
  # `x + swindow > D` (legitimate when the secondary censoring interval
  # straddles D) are still valid. The likelihood becomes
  # `P(X in [x, min(x + swindow, D)] | L <= X <= D)`, which equals the usual
  # interval probability when `x + swindow <= D` (the parametric default) and
  # captures the residual mass between `x` and `D` otherwise.
  upper_raw <- x + swindow
  upper <- pmin(upper_raw, D)
  if (is.finite(D) && any(upper_raw > D)) {
    message(
      "Upper truncation point is greater than D. It is ",
      max(upper_raw),
      " and D is ",
      D,
      "; clipping the upper end of secondary intervals at D."
    )
  }

  # Compute raw (unnormalised) CDFs for all unique points so PMF differences
  # below can be normalised with the truncation-aware F_cens(L) and F_cens(D).
  unique_points <- sort(unique(c(x, upper)))
  if (length(unique_points) == 0) {
    return(rep(0, length(x)))
  }
  cdfs <- pcens_cdf(object, unique_points, pwindow)
  # Match `pprimarycensored(L = -Inf, D = Inf)` at infinite points.
  cdfs[unique_points == -Inf] <- 0
  cdfs[unique_points == Inf] <- 1

  result <- cdfs[match(upper, unique_points)] -
    cdfs[match(x, unique_points)]

  # Fast path: with no truncation on either side the raw PMF needs no
  # renormalisation, so skip the two extra CDF lookups below.
  if (!(is.infinite(L) && is.infinite(D))) {
    cdf_D <- .pcens_cdf_at(object, D, pwindow, unique_points, cdfs, 1)
    cdf_L <- .pcens_cdf_at(object, L, pwindow, unique_points, cdfs, 0)

    # Divide by (F(D) - F(L)). Skip the division when the normaliser is 1
    # (e.g. a finite `L` that sits below the support of the delay, so
    # `F_cens(L) = 0`, paired with `D = Inf` where `F_cens(D) = 1`).
    normaliser <- cdf_D - cdf_L
    if (normaliser != 1) {
      result <- result / normaliser
    }
  }

  # Ensure non-negative values (can become slightly negative due to
  # floating-point precision when computing CDF differences)
  result <- pmax(0, result)

  if (log) {
    return(log(result))
  } else {
    return(result)
  }
}

#' Primary event censored CDF at a truncation point
#'
#' Returns the primary event censored CDF, before truncation, at a single
#' truncation point, reusing an already computed value where possible.
#'
#' @param object A `pcens` object.
#'
#' @param bound Numeric truncation point (`L` or `D`).
#'
#' @param pwindow Primary event window.
#'
#' @param points Numeric vector of points at which `cdfs` was computed.
#'
#' @param cdfs Numeric vector of CDF values at `points`.
#'
#' @param inf_value CDF value to return when `bound` is infinite.
#'
#' @return A single numeric CDF value.
#'
#' @keywords internal
.pcens_cdf_at <- function(object, bound, pwindow, points, cdfs, inf_value) {
  if (is.infinite(bound)) {
    return(inf_value)
  }
  idx <- match(bound, points)
  if (!is.na(idx)) {
    return(cdfs[[idx]])
  }
  pcens_cdf(object, bound, pwindow)
}
