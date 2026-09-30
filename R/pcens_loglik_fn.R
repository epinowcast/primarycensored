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
#' @param check Logical; if `TRUE` (the default) `dprimary` is validated with
#'   [check_dprimary()] when the function is built, and `pdist` is validated
#'   with [check_pdist()] on the first valid call for each set of parameter
#'   names. Set to `FALSE` to skip both.
#'
#' @details
#' The returned function is called as `ll(...)` with the delay distribution
#' parameters as named arguments, for example `ll(shape = 2, rate = 1)`.
#' Each call starts from the fixed parameters given in `...` here, so
#' parameters from an earlier call are not carried over.
#' A misspelt parameter name raises an error.
#' Repeated delays are evaluated once.
#' The result matches `log(dprimarycensored())` called on each observation.
#'
#' Parameters that make `pdist` return `NaN` give `NaN`.
#' Errors raised by `pdist` are not caught, so wrap the call in `tryCatch()`
#' if an optimiser should see a missing value instead.
#' Zero probabilities are returned as `-Inf`.
#'
#' The construction checks that every observation satisfies `L <= x < D`.
#' A message is given at construction when any secondary interval extends
#' past `D`. None is given when the returned function is called.
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
  .message_if_clipped(x, swindow, D)
  max_D <- if (n > 0L) max(D) else Inf

  # Parameter names are checked when they change, as `update()` with
  # `check = FALSE` would otherwise hide a misspelt name. `pdist` is checked
  # once a call gives a result without missing values, so invalid starting
  # parameters give `NaN` rather than an error from `check_pdist()`.
  state <- new.env(parent = emptyenv())
  state$checked <- FALSE
  state$names <- NULL
  state$pdist_pending <- isTRUE(check)
  function(...) {
    nms <- names(list(...))
    if (!state$checked || !identical(nms, state$names)) {
      obj <- update(base, ...)
      state$names <- nms
      state$checked <- TRUE
      state$pdist_pending <- isTRUE(check)
    } else {
      obj <- update(base, ..., check = FALSE)
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
#' @keywords internal
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
  L <- rep_len(L, n)
  D <- rep_len(D, n)
  .check_truncation_bounds_df(data.frame(L = L, D = D), "L", "D")
  below <- which(x < L)
  if (length(below) > 0L) {
    stop(
      "Some values of x are below L. First invalid row: ", below[1],
      " (x = ", x[below[1]], ", L = ", L[below[1]], "). Resolve this by ",
      "filtering x to only include values >= L.",
      call. = FALSE
    )
  }
  above <- which(is.finite(D) & x >= D)
  if (length(above) > 0L) {
    stop(
      "Upper truncation point is greater than D. First invalid row: ",
      above[1], " (x = ", x[above[1]], ", D = ", D[above[1]], "). Under ",
      "truncation at D no event with latent value >= D is observable; ",
      "resolve this by filtering x to values strictly less than D.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Message when secondary intervals extend past the upper truncation point
#'
#' @inheritParams pcens_loglik_fn
#'
#' @return `NULL` invisibly. Called for its message.
#'
#' @keywords internal
.message_if_clipped <- function(x, swindow, D) {
  upper <- x + swindow
  over <- is.finite(D) & upper > D
  if (any(over)) {
    message(
      "Upper truncation point is greater than D for ", sum(over),
      " observation(s); clipping the upper end of their secondary ",
      "intervals at D."
    )
  }
  invisible(NULL)
}

#' Group observations that share settings and collapse repeated delays
#'
#' Groups rows by their `pwindow`, `swindow`, `L` and `D`. Within each group
#' the unique values of `x` are kept with a map back to the rows.
#' The primary event censored CDF depends on the primary event window and the
#' point only, so the points at which it is needed are pooled over every
#' group that shares a `pwindow`. Each group stores the positions of its own
#' points in the pooled set. This means each endpoint is evaluated once per
#' `pwindow` and reused by every group that needs it, as in the vectorised
#' Stan PMF, and [.pcens_pmf_groups()] only has to evaluate the CDF once per
#' set.
#'
#' @param x Numeric vector of delays.
#'
#' @param pwindow,swindow,L,D Per-observation settings, each of length 1
#'   or `length(x)`.
#'
#' @return A list with `groups` and `sets`, both empty if `x` is empty.
#'   `groups` has one element per group. Each is a list with `idx`, the rows
#'   of the group or `NULL` if the group is every row, `x`, the unique delays
#'   of the group, `map`, the position of each row of the group in `x` or
#'   `NULL` if `x` has no repeats, and the group's `pwindow`, `swindow`, `L`
#'   and `D`. It also has `set`, the position of the group's set in `sets`,
#'   and the positions within that set's points as made by
#'   [.pcens_group_positions()].
#'   `sets` has one element per unique `pwindow`. Each is a list with
#'   `pwindow`, `points`, the sorted unique points at which the CDF is needed
#'   (empty if there is nothing to evaluate), and `at_minf` and `at_inf`, the
#'   positions of `-Inf` and `Inf` in `points`, where the CDF is known to be 0
#'   and 1.
#'
#' @keywords internal
.pcens_row_groups <- function(x, pwindow, swindow, L, D) {
  n <- length(x)
  if (n == 0L) {
    return(list(groups = list(), sets = list()))
  }
  settings <- list(pwindow = pwindow, swindow = swindow, L = L, D = D)
  # Integer code of each row's combination of settings. Re-matching after
  # each column keeps the codes small.
  id <- rep.int(1L, n)
  for (s in settings) {
    if (length(s) == 1L) {
      next
    }
    code <- match(s, unique(s))
    k <- max(code)
    if (k > 1L) {
      key <- as.numeric(id - 1L) * k + code
      id <- match(key, unique(key))
    }
  }
  rows <- list(NULL)
  if (max(id) > 1L) {
    rows <- unname(split(seq_len(n), id))
  }
  groups <- lapply(rows, function(idx) {
    first <- if (is.null(idx)) 1L else idx[[1L]]
    setting <- function(s) s[[if (length(s) == 1L) 1L else first]]
    xs <- if (is.null(idx)) x else x[idx]
    ux <- unique(xs)
    map <- NULL
    if (length(ux) < length(xs)) {
      map <- match(xs, ux)
    }
    list(
      idx = idx,
      x = ux,
      map = map,
      pwindow = setting(pwindow),
      swindow = setting(swindow),
      L = setting(L),
      D = setting(D)
    )
  })
  .pcens_share_points(groups)
}

#' Pool the CDF points of groups that share a primary event window
#'
#' @param groups List of groups as made in [.pcens_row_groups()], before
#'   their points are pooled.
#'
#' @return A list with `groups` and `sets`, as described in
#'   [.pcens_row_groups()].
#'
#' @keywords internal
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
    needed <- as.numeric(needed)
    # Skip the sort when the points are already in order (e.g. x = 0:n)
    if (is.unsorted(needed)) {
      needed <- sort(needed)
    }
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
#' Works out what [pcens_pmf()] would repeat on every call. These are the
#' delays and the clipped upper ends of their secondary intervals, and any
#' finite truncation points, which are in the same set so that one CDF
#' evaluation serves them all.
#'
#' @param x Numeric vector of unique delays of a group.
#'
#' @param swindow,L,D Secondary window and truncation points of the group,
#'   each a single value.
#'
#' @return A list with `exact`, whether the group has a zero-width secondary
#'   window and so contributes densities, `upper`, the upper end of each
#'   secondary interval clipped at `D` (`NULL` if `exact`), `truncated`,
#'   whether the PMF is normalised, and `needed`, the points, which may
#'   repeat, at which the CDF is needed.
#'
#' @keywords internal
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
    truncated = !(is.infinite(L) && is.infinite(D)),
    needed = c(if (!exact) c(x, upper), bounds[is.finite(bounds)])
  )
}

#' Positions of a group's CDF points in a pooled set of points
#'
#' @param ends Output of [.pcens_group_ends()] for the group.
#'
#' @param group The group, with its delays `x` and truncation points `L`
#'   and `D`.
#'
#' @param points Sorted unique points of the set the group belongs to.
#'
#' @return A list with `exact` and `truncated` as in [.pcens_group_ends()],
#'   `lower` and `upper`, the positions in `points` of the ends of each
#'   secondary interval (`NULL` if `exact`), and `pos_L` and `pos_D`, the
#'   positions of `L` and `D` in `points` (`NA` if infinite).
#'
#' @keywords internal
.pcens_group_positions <- function(ends, group, points) {
  position <- function(bound) {
    if (is.finite(bound)) match(bound, points) else NA_integer_
  }
  list(
    exact = ends$exact,
    truncated = ends$truncated,
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
#' @param grouped Groups of observations and their shared CDF points as
#'   made by [.pcens_row_groups()].
#'
#' @param n Number of observations.
#'
#' @return Numeric vector of length `n` with the PMF (or density, where
#'   `swindow = 0`) of each observation, in row order. Not on the log scale.
#'
#' @keywords internal
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
#' Gives the same values as [pcens_pmf()] for the unique delays of the group
#' and copies them to the rows of the group. The CDF values are taken from
#' the set the group belongs to, and the positions needed to difference and
#' normalise them from the group, rather than worked out on each call.
#' A message about clipping at `D` is not given.
#' The clipping at `D`, the normalisation by `F(D) - F(L)` and the
#' non-negative clamp repeat [pcens_pmf()], and the input checks of
#' [.check_row_inputs()] repeat [.check_truncation_bounds_df()] and the
#' checks of [dprimarycensored()]. Keep them in step with those.
#'
#' @inheritParams .pcens_pmf_groups
#'
#' @param group One element of `groups` in the list made by
#'   [.pcens_row_groups()].
#'
#' @param cdfs Numeric vector of the CDF at the points of the group's set.
#'
#' @return Numeric vector with one value per row of the group.
#'
#' @keywords internal
.pcens_pmf_group <- function(object, group, cdfs) {
  if (group$exact) {
    pmf <- .pcens_density(object, group$x, group$pwindow)
  } else {
    pmf <- cdfs[group$upper] - cdfs[group$lower]
  }
  if (group$truncated) {
    cdf_D <- if (is.na(group$pos_D)) 1 else cdfs[[group$pos_D]]
    cdf_L <- if (is.na(group$pos_L)) 0 else cdfs[[group$pos_L]]
    normaliser <- cdf_D - cdf_L
    # A missing normaliser, from invalid parameters, gives missing values
    if (!isTRUE(normaliser == 1)) {
      pmf <- pmf / normaliser
    }
  }
  # Ensure non-negative values, keeping missing values
  pmf[which(pmf < 0)] <- 0
  if (is.null(group$map)) pmf else pmf[group$map]
}
