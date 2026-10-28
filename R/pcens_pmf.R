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
  .check_truncation_bounds(L, D, vector = TRUE)
  if (length(x) == 0) {
    return(numeric(0))
  }
  n <- max(lengths(list(x, pwindow, swindow, L, D)))
  if (length(x) < n) {
    x <- rep_len(x, n)
  }

  if (any(x < L)) {
    below <- which(x < L)
    i <- below[which.min(x[below])]
    stop(
      "Some values of x are below L. Minimum x is ", x[i],
      " and L is ", rep_len(L, n)[i],
      if (length(L) > 1) c(" (row ", i, ")"),
      ". Resolve this by filtering x to only include values >= L.",
      call. = FALSE
    )
  }

  if (any(is.finite(D) & x >= D)) {
    above <- which(is.finite(D) & x >= D)
    i <- above[which.max(x[above])]
    stop(
      "Upper truncation point is greater than D. Maximum x is ", x[i],
      " and D is ", rep_len(D, n)[i],
      if (length(D) > 1) c(" (row ", i, ")"),
      ". Under truncation at D no event with latent value >= D is ",
      "observable; resolve this by filtering x to values strictly less than D.",
      call. = FALSE
    )
  }

  # Clip the upper end of each secondary interval at D
  upper <- x + swindow
  if (any(upper > D)) {
    clipped <- which(upper > D)
    i <- clipped[which.max(upper[clipped])]
    message(
      "Upper truncation point is greater than D. It is ", upper[i],
      " and D is ", rep_len(D, n)[i],
      if (length(D) > 1) c(" (row ", i, ")"),
      "; clipping the upper end of secondary intervals at D."
    )
    upper <- pmin(upper, D)
  }

  # Rows with a zero-width secondary window contribute a density, so their
  # interval ends are set to -Inf, which needs no CDF evaluation
  exact <- swindow == 0
  lower <- x
  if (any(exact)) {
    exact <- rep_len(exact, n)
    lower[exact] <- -Inf
    upper[exact] <- -Inf
  }
  cdf <- .pcens_cdf_pooled(object, list(lower, upper, L, D), pwindow)
  result <- cdf[[2]] - cdf[[1]]
  if (any(exact)) {
    pwindow <- rep_len(pwindow, n)
    for (rows in .row_groups(pwindow[exact])) {
      rows <- which(exact)[rows]
      result[rows] <- .pcens_density(object, x[rows], pwindow[rows[1]])
    }
  }

  # Normalise by F(D) - F(L), skipped when it is 1 without truncation
  normaliser <- cdf[[4]] - cdf[[3]]
  if (any(normaliser != 1)) {
    result <- result / normaliser
  }
  result <- pmax(0, result)

  if (log) {
    return(log(result))
  } else {
    return(result)
  }
}

#' Primary event censored CDF at pooled points
#'
#' Evaluates [pcens_cdf()] once at the unique finite points of each set of
#' rows sharing a `pwindow`, as the CDF depends on no other setting.
#'
#' @param object A `pcens` object.
#'
#' @param points List of numeric vectors of points, recycled with `pwindow`.
#'
#' @param pwindow Primary event window.
#'
#' @inheritParams pcens_cdf
#'
#' @return A list like `points` of CDF values, which are 0 at `-Inf` and 1
#'   at `Inf`.
#'
#' @keywords internal
.pcens_cdf_pooled <- function(object, points, pwindow,
                              use_numeric = FALSE) {
  if (length(pwindow) > 1L) {
    n <- max(lengths(c(points, list(pwindow))))
    pwindow <- rep_len(pwindow, n)
    # Single points, such as L and D, are kept as they are in every batch
    vary <- lengths(points) > 1L
    points[vary] <- lapply(points[vary], rep_len, n)
    cdf <- rep(list(numeric(n)), length(points))
    for (rows in .row_groups(pwindow)) {
      batch <- points
      batch[vary] <- lapply(points[vary], `[`, rows)
      batch <- .pcens_cdf_pooled(
        object, batch, pwindow[[rows[[1L]]]], use_numeric
      )
      for (k in seq_along(cdf)) {
        cdf[[k]][rows] <- batch[[k]]
      }
    }
    return(cdf)
  }
  all_points <- unlist(points, use.names = FALSE)
  pts <- unique(all_points)
  # Infinite points are not evaluated, as some analytical methods give NaN
  finite <- !is.infinite(pts)
  cdf <- as.numeric(pts == Inf)
  if (any(finite)) {
    cdf[finite] <- pcens_cdf(object, pts[finite], pwindow, use_numeric)
  }
  if (length(pts) < length(all_points)) {
    cdf <- cdf[match(all_points, pts)]
  }
  if (length(all_points) == length(points)) {
    return(as.list(cdf))
  }
  # Split back into the elements of `points`
  done <- 0
  for (k in seq_along(points)) {
    len <- length(points[[k]])
    points[[k]] <- cdf[done + seq_len(len)]
    done <- done + len
  }
  points
}

#' Group rows that share values
#'
#' @param ... Vectors of the same length, one element per row.
#'
#' @return A list of integer row indices, one element per unique combination
#'   of values, in order of first appearance.
#'
#' @keywords internal
.row_groups <- function(...) {
  keys <- list(...)
  id <- match(keys[[1L]], unique(keys[[1L]]))
  for (key in keys[-1L]) {
    values <- unique(key)
    # Renumber the combinations so far, so ids stay at most the row count
    id <- (id - 1L) * length(values) + match(key, values)
    id <- match(id, unique(id))
  }
  unname(split(seq_along(id), id))
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
