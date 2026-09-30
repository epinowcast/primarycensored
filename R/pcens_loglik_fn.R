#' Build a reusable log-likelihood function for primary event censored data
#'
#' `r lifecycle::badge("experimental")`
#'
#' Sets up the log-likelihood of a set of observations once and returns a
#' function of the delay distribution parameters. It is intended for use
#' inside an optimiser or sampler, where the same data are evaluated many
#' times with different parameters.
#' The setup that [dprimarycensored()] repeats on every call is done once
#' here.
#'
#' @inheritParams dprimarycensored
#'
#' @param x Vector of observed delays, the lower bounds of the secondary
#'   event intervals.
#'
#' @param pwindow,swindow,L,D Primary event window, secondary event window,
#'   lower truncation point and upper truncation point. Each is either a
#'   single value, used for every observation, or a vector with one value per
#'   element of `x`. See [dprimarycensored()] for their meaning, including
#'   zero-width windows and secondary intervals that extend past `D`.
#'
#' @param ... Delay distribution parameters held fixed, for example
#'   `sdlog = 0.5` for a log-normal delay. They are merged with, and can be
#'   overridden by, the parameters given to the returned function.
#'
#' @param check Logical; if `TRUE` (the default), validate `dprimary` and
#'   `pdist`. Set to `FALSE` to skip the checks.
#'
#' @details
#' The returned function takes the delay distribution parameters as named
#' scalars, other than vector parameters such as `boundaries` and `pmf`.
#' Errors from `pdist` are not caught.
#' The message about secondary intervals that extend past `D` is given only
#' when the function is built.
#'
#' @return A function of the delay distribution parameters that returns a
#'   numeric vector of log-likelihood contributions, one per element of `x`
#'   in the order of `x`. Sum them for the total log-likelihood.
#'
#' @family modelhelpers
#' @seealso [fitdistdoublecens()] which uses the same machinery,
#'   [dprimarycensored()] for a single evaluation, and [new_pcens()] and
#'   [update()][update.pcens()] for the underlying objects.
#'
#' @export
#' @examples
#' set.seed(1)
#' pw <- 2
#' x <- rprimarycensored(
#'   200, rlnorm,
#'   meanlog = 1.3, sdlog = 0.5,
#'   pwindow = pw, swindow = 1, D = 20
#' )
#' ll <- pcens_loglik_fn(
#'   floor(x), plnorm,
#'   pwindow = pw, swindow = 1, D = 20
#' )
#' # Log-likelihood contributions for one set of parameters
#' head(ll(meanlog = 1.3, sdlog = 0.5))
#'
#' # Maximise the log-likelihood directly
#' fit <- optim(
#'   c(meanlog = 1, log_sdlog = 0),
#'   function(par) -sum(ll(meanlog = par[[1]], sdlog = exp(par[[2]])))
#' )
#' c(fit$par[[1]], exp(fit$par[[2]]))
#'
#' # Observations can have their own windows and truncation points
#' ll_mixed <- pcens_loglik_fn(
#'   c(1, 2, 3, 4), pgamma,
#'   pwindow = c(1, 1, 2, 2), swindow = c(1, 0, 1, 1), D = c(Inf, Inf, 10, 10)
#' )
#' ll_mixed(shape = 2, rate = 1)
pcens_loglik_fn <- function(
  x,
  pdist,
  pwindow = 1,
  swindow = 1,
  L = -Inf,
  D = Inf,
  dprimary = dunif,
  primary_args = NULL,
  pprimary = NULL,
  ...,
  check = TRUE
) {
  .check_row_inputs(x, pwindow, swindow, L, D)
  n <- length(x)
  pwindow <- rep_len(pwindow, n)
  swindow <- rep_len(swindow, n)
  L <- rep_len(L, n)
  D <- rep_len(D, n)

  primary_args <- .resolve_primary_args(
    primary_args, NULL, "pcens_loglik_fn"
  )
  base <- .build_pcens(
    pdist, dprimary, primary_args, pprimary, list(),
    pwindow = NULL, D = NULL, check = FALSE
  )
  if (...length() > 0L) {
    # Checks the names of the fixed parameters against pdist
    base <- update(base, ...)
  }
  if (isTRUE(check)) {
    for (pw in unique(pwindow)) {
      check_dprimary(dprimary, pw, primary_args)
    }
  }

  groups <- .pcens_row_groups(x, pwindow, swindow, L, D)
  .message_if_clipped(x + swindow, D)
  max_D <- if (n > 0L) max(D) else Inf

  # Check `pdist` on the first call without missing values, so invalid
  # starting values give `NaN` and not an error
  state <- new.env(parent = emptyenv())
  state$names <- NA_character_
  state$pdist_pending <- isTRUE(check)
  function(...) {
    nms <- names(list(...))
    if (identical(nms, state$names)) {
      obj <- update(base, ..., check = FALSE)
    } else {
      obj <- update(base, ...)
      state$names <- nms
      state$pdist_pending <- isTRUE(check)
    }
    out <- log(.pcens_pmf_groups(obj, groups, n))
    if (state$pdist_pending && !anyNA(out)) {
      do.call(check_pdist, c(list(obj$pdist, D = max_D), obj$args))
      state$pdist_pending <- FALSE
    }
    out
  }
}

#' Validate the per-observation inputs of `pcens_loglik_fn()`
#'
#' @inheritParams pcens_loglik_fn
#'
#' @return `NULL` invisibly. Called for its errors.
#'
#' @noRd
.check_row_inputs <- function(x, pwindow, swindow, L, D) {
  if (!is.numeric(x)) {
    stop("x must be numeric.", call. = FALSE)
  }
  if (anyNA(x)) {
    stop("x must not contain missing values.", call. = FALSE)
  }
  n <- length(x)
  inputs <- list(pwindow = pwindow, swindow = swindow, L = L, D = D)
  for (nm in names(inputs)) {
    if (!is.numeric(inputs[[nm]])) {
      stop(nm, " must be numeric.", call. = FALSE)
    }
    len <- length(inputs[[nm]])
    if (len != 1L && len != n) {
      stop(
        nm, " must have length 1 or the length of x (", n, ").",
        call. = FALSE
      )
    }
    if (anyNA(inputs[[nm]])) {
      stop(nm, " must not contain missing values.", call. = FALSE)
    }
  }
  for (nm in c("pwindow", "swindow")) {
    if (any(inputs[[nm]] < 0)) {
      stop(nm, " must be non-negative.", call. = FALSE)
    }
  }
  if (n == 0L) {
    return(invisible(NULL))
  }
  .check_truncation_bounds_df(
    data.frame(L = rep_len(L, n), D = rep_len(D, n)), "L", "D"
  )
  .check_x_bounds(x, L, D)
}

#' Group observations that share settings and collapse repeated delays
#'
#' Rows are grouped by `pwindow`, `swindow`, `L` and `D`, and each group keeps
#' its unique delays with a map back to its rows.
#' The CDF points of groups that share a `pwindow` are pooled, so each point
#' is evaluated once per `pwindow`.
#'
#' @param x Numeric vector of delays.
#'
#' @param pwindow,swindow,L,D Per-observation settings, each of length
#'   `length(x)`.
#'
#' @return A list with `groups`, one element per group, and `sets`, one
#'   element per unique `pwindow` with its pooled CDF `points`.
#'
#' @noRd
.pcens_row_groups <- function(x, pwindow, swindow, L, D) {
  n <- length(x)
  if (n == 0L) {
    return(list(groups = list(), sets = list()))
  }
  settings <- list(pwindow = pwindow, swindow = swindow, L = L, D = D)
  groups <- lapply(.param_groups(settings, names(settings)), function(g) {
    idx <- g$idx
    xs <- if (is.null(idx)) x else x[idx]
    ux <- unique(xs)
    map <- NULL
    if (length(ux) < length(xs)) {
      map <- match(xs, ux)
    }
    c(g, list(x = ux, map = map))
  })
  .pcens_share_points(groups)
}

#' Pool the CDF points of groups that share a primary event window
#'
#' @param groups Groups as made by `.pcens_row_groups()`.
#'
#' @return As `.pcens_row_groups()`, with the pooled points and the position
#'   of each group's points in them.
#'
#' @noRd
.pcens_share_points <- function(groups) {
  pwindows <- vapply(groups, function(g) g$pwindow, numeric(1))
  unique_pw <- unique(pwindows)
  set_id <- match(pwindows, unique_pw)
  ends <- lapply(groups, function(g) {
    .pcens_group_ends(g$x, g$swindow, g$L, g$D)
  })
  sets <- lapply(seq_along(unique_pw), function(k) {
    members <- which(set_id == k)
    needed <- unique(unlist(
      lapply(ends[members], function(e) e$needed),
      use.names = FALSE
    ))
    needed <- sort(as.numeric(needed))
    list(
      pwindow = unique_pw[[k]],
      points = needed,
      at_minf = which(needed == -Inf),
      at_inf = which(needed == Inf)
    )
  })
  for (i in seq_along(groups)) {
    groups[[i]] <- c(
      groups[[i]],
      list(set = set_id[[i]]),
      .pcens_group_positions(
        ends[[i]], groups[[i]], sets[[set_id[[i]]]]$points
      )
    )
  }
  list(groups = groups, sets = sets)
}

#' Points at which the CDF is needed for one group of observations
#'
#' These are the delays, the upper ends of their secondary intervals clipped
#' at `D`, and any finite truncation points.
#'
#' @param x Unique delays of a group.
#'
#' @param swindow,L,D Settings of the group, each a single value.
#'
#' @return A list with `exact` (zero-width secondary window), `upper` and
#'   `needed`, the points, which may repeat.
#'
#' @noRd
.pcens_group_ends <- function(x, swindow, L, D) {
  exact <- swindow == 0
  upper <- NULL
  if (!exact) {
    upper <- x + swindow
    if (is.finite(D)) {
      upper <- pmin(upper, D)
    }
  }
  bounds <- c(L, D)
  list(
    exact = exact,
    upper = upper,
    needed = c(if (!exact) c(x, upper), bounds[is.finite(bounds)])
  )
}

#' Positions of a group's CDF points in a pooled set of points
#'
#' @param ends Output of `.pcens_group_ends()` for the group.
#'
#' @param group The group.
#'
#' @param points Sorted unique points of the group's set.
#'
#' @return A list with `exact`, `lower` and `upper`, the positions of the ends
#'   of each secondary interval, and `pos_L` and `pos_D`, the positions of
#'   `L` and `D` (`NA` if infinite).
#'
#' @noRd
.pcens_group_positions <- function(ends, group, points) {
  position <- function(bound) {
    if (is.finite(bound)) match(bound, points) else NA_integer_
  }
  list(
    exact = ends$exact,
    lower = if (ends$exact) NULL else match(group$x, points),
    upper = if (ends$exact) NULL else match(ends$upper, points),
    pos_L = position(group$L),
    pos_D = position(group$D)
  )
}

#' Evaluate the primary event censored PMF for grouped observations
#'
#' @param object A `pcens` object.
#'
#' @param grouped Output of `.pcens_row_groups()`.
#'
#' @param n Number of observations.
#'
#' @return Numeric vector of length `n` with the PMF, or the density where
#'   `swindow = 0`, of each observation in row order.
#'
#' @noRd
.pcens_pmf_groups <- function(object, grouped, n) {
  cdfs <- lapply(grouped$sets, function(set) {
    cdf <- numeric(0)
    if (length(set$points) > 0L) {
      cdf <- pcens_cdf(object, set$points, set$pwindow)
      # Some analytical methods return NaN at Inf
      cdf[set$at_minf] <- 0
      cdf[set$at_inf] <- 1
    }
    cdf
  })
  groups <- grouped$groups
  if (length(groups) == 1L && is.null(groups[[1L]]$idx)) {
    return(.pcens_pmf_group(object, groups[[1L]], cdfs[[1L]]))
  }
  result <- numeric(n)
  for (g in groups) {
    result[g$idx] <- .pcens_pmf_group(object, g, cdfs[[g$set]])
  }
  result
}

#' Evaluate the primary event censored PMF for one group of observations
#'
#' Gives the values of `pcens_pmf()` for the unique delays of the group,
#' copied to its rows, from the CDF values at the points of the group's set.
#'
#' @param object A `pcens` object.
#'
#' @param group One element of `groups` from `.pcens_row_groups()`.
#'
#' @param cdfs CDF at the points of the group's set.
#'
#' @return Numeric vector with one value per row of the group.
#'
#' @noRd
.pcens_pmf_group <- function(object, group, cdfs) {
  if (group$exact) {
    pmf <- .pcens_density(object, group$x, group$pwindow)
  } else {
    pmf <- cdfs[group$upper] - cdfs[group$lower]
  }
  cdf_D <- if (is.na(group$pos_D)) 1 else cdfs[[group$pos_D]]
  cdf_L <- if (is.na(group$pos_L)) 0 else cdfs[[group$pos_L]]
  pmf <- .pcens_normalise(pmf, cdf_D, cdf_L)
  if (is.null(group$map)) pmf else pmf[group$map]
}
