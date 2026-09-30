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
#' `min(x + swindow, D)`, normalised over \[L, D\]. Where `swindow = 0` the
#' primary event censored density at `x` is returned instead. See
#' [dprimarycensored()] for the details.
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
  if (length(x) == 0) {
    return(numeric(0))
  }

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
  upper <- x + swindow
  if (is.finite(D) && any(upper > D)) {
    upper_raw <- upper
    upper <- pmin(upper_raw, D)
    message(
      "Upper truncation point is greater than D. It is ",
      max(upper_raw),
      " and D is ",
      D,
      "; clipping the upper end of secondary intervals at D."
    )
  }

  # Rows with a zero-width secondary window contribute a density
  exact <- rep_len(swindow == 0, length(x))
  result <- numeric(length(x))

  # Compute CDFs for all unique points
  unique_points <- unique(c(x[!exact], upper[!exact]))
  # Skip the sort when the points are already in order (e.g. x = 0:n)
  if (anyNA(unique_points) || is.unsorted(unique_points)) {
    unique_points <- sort(unique_points)
  }
  cdfs <- numeric(0)
  if (length(unique_points) > 0) {
    cdfs <- pcens_cdf(object, unique_points, pwindow)
    # Some analytical methods return NaN at Inf
    cdfs[unique_points == -Inf] <- 0
    cdfs[unique_points == Inf] <- 1

    result[!exact] <- cdfs[match(upper[!exact], unique_points)] -
      cdfs[match(x[!exact], unique_points)]
  }
  if (any(exact)) {
    result[exact] <- .pcens_density(object, x[exact], pwindow)
  }

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

#' Primary event censored density
#'
#' Computes the density of the primary event censored delay, the derivative
#' of [pcens_cdf()] in `x`. This is the contribution of an observation with a
#' zero-width secondary window.
#'
#' With `pwindow = 0` this is the delay density. With a uniform primary
#' event distribution it is \eqn{(F(x) - F(x - pwindow)) / pwindow}, which
#' only needs the delay CDF. Otherwise the delay density is integrated
#' against the primary event density over \[0, pwindow\].
#'
#' @param object A `pcens` object as created by [new_pcens()].
#'
#' @param x Vector of points at which to evaluate the density.
#'
#' @param pwindow Primary event window.
#'
#' @return Vector of densities, not normalised for truncation.
#'
#' @keywords internal
.pcens_density <- function(object, x, pwindow) {
  exact_primary <- .is_exact_window(pwindow)
  # The class ends in "_dunif" for a uniform primary, as set by
  # .format_class() and used to dispatch the analytical pcens_cdf() methods
  if (!exact_primary && endsWith(class(object)[1], "_dunif")) {
    return(
      (.delay_cdf(object, x) - .delay_cdf(object, x - pwindow)) / pwindow
    )
  }
  ddist <- .lookup_ddist(object$pdist)
  if (exact_primary) {
    return(do.call(ddist, c(list(x), object$args)))
  }
  vapply(
    x,
    function(d) {
      integrand <- function(p) {
        do.call(ddist, c(list(d - p), object$args)) *
          do.call(
            object$dprimary,
            c(list(x = p, min = 0, max = pwindow), object$primary_args)
          )
      }
      stats::integrate(integrand, lower = 0, upper = pwindow)$value
    },
    numeric(1)
  )
}

#' Look up the density matching a delay CDF
#'
#' Finds the `d` function paired with a `p` function by name, for example
#' `dgamma()` for `pgamma()`. The name is taken from the `"name"` attribute
#' of `pdist` or inferred with [.dist_name()]. The density is searched
#' for from the environment of `pdist`, then in `stats` and
#' `primarycensored`.
#'
#' @param pdist Delay distribution CDF.
#'
#' @return The density function. An error is raised if none is found.
#'
#' @keywords internal
.lookup_ddist <- function(pdist) {
  pdist_name <- .dist_name(pdist)
  ddist_name <- sub("^p", "d", pdist_name)
  if (startsWith(pdist_name, "p")) {
    envs <- list(
      environment(pdist), asNamespace("stats"),
      asNamespace("primarycensored")
    )
    # Primitives have no environment
    for (env in Filter(Negate(is.null), envs)) {
      ddist <- get0(ddist_name, envir = env, mode = "function")
      if (!is.null(ddist)) {
        return(ddist)
      }
    }
  }
  stop(
    "A zero-width secondary window (swindow = 0) needs the delay density, ",
    "but no density '", ddist_name, "' matching pdist ('", pdist_name,
    "') was found. Define it where pdist is defined, or name pdist with ",
    "add_name_attribute().",
    call. = FALSE
  )
}
