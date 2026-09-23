#' Compute primary event censored PMF
#'
#' Computes the primary event censored PMF for a `pcens` object as created
#' by [new_pcens()]. Secondary event windows and truncation are handled as in
#' [dprimarycensored()].
#'
#' @inheritParams dprimarycensored
#'
#' @param object A `pcens` object as created by [new_pcens()].
#'
#' @param ... Additional arguments passed to methods.
#'
#' @inherit dprimarycensored return
#'
#' @family pcens
#'
#' @export
pcens_pmf <- function(
    object,
    x,
    pwindow,
    swindow = 1,
    L = -Inf,
    D = Inf,
    log = FALSE,
    ...) {
  UseMethod("pcens_pmf")
}

#' Default method for computing primary event censored PMF
#'
#' Computes the PMF by differencing [pcens_cdf()] at `x` and
#' `min(x + swindow, D)`, normalised over \[L, D\]. See [dprimarycensored()]
#' for the details.
#'
#' @inheritParams pcens_pmf
#' @inheritParams dprimarycensored
#'
#' @inherit dprimarycensored return
#'
#' @family pcens
#'
#' @export
#' @examples
#' obj <- new_pcens(
#'   pdist = pgamma, dprimary = dunif,
#'   primary_args = list(min = 0, max = 1),
#'   shape = 3, scale = 2
#' )
#' pcens_pmf(obj, x = 0:9, pwindow = 1, D = 10)
pcens_pmf.default <- function(
    object,
    x,
    pwindow,
    swindow = 1,
    L = -Inf,
    D = Inf,
    log = FALSE,
    ...) {
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

  # Clip the upper end of each secondary interval at D
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

  # Compute CDFs for all unique points
  unique_points <- sort(unique(c(x, upper)))
  if (length(unique_points) == 0) {
    return(rep(0, length(x)))
  }
  cdfs <- pcens_cdf(object, unique_points, pwindow)
  # Some analytical methods return NaN at Inf
  cdfs[unique_points == -Inf] <- 0
  cdfs[unique_points == Inf] <- 1

  result <- cdfs[match(upper, unique_points)] -
    cdfs[match(x, unique_points)]

  # Normalise by F(D) - F(L) when truncated
  if (!(is.infinite(L) && is.infinite(D))) {
    cdf_D <- .pcens_cdf_at(object, D, pwindow, unique_points, cdfs, 1)
    cdf_L <- .pcens_cdf_at(object, L, pwindow, unique_points, cdfs, 0)
    normaliser <- cdf_D - cdf_L
    if (normaliser != 1) {
      result <- result / normaliser
    }
  }

  # Ensure non-negative values
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
