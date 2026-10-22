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
#' @param check Logical; if `TRUE` (the default), validate `dprimary` when
#'   the function is built and `pdist` on the first call that gives no
#'   missing values, so invalid starting values give `NaN` rather than an
#'   error. Set to `FALSE` to skip the checks.
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
#' @seealso [dprimarycensored()] for a single evaluation, and [new_pcens()]
#'   and [update()][update.pcens()] for the underlying objects.
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

  if (is.null(primary_args)) {
    primary_args <- list()
  }
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

  # `pdist` is checked once per set of parameter names, on the first call
  # without missing values
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
#' Each group keeps the positions of its points in the pool.
#'
#' @inheritParams .check_row_inputs
#'
#' @param pwindow,swindow,L,D Per-observation settings, each of length
#'   `length(x)`.
#'
#' @return A list with `groups`, one element per group, and `sets`, one
#'   element per unique `pwindow` with its pooled CDF `points`.
#'
#' @noRd
.pcens_row_groups <- function(x, pwindow, swindow, L, D) {
  settings <- list(pwindow = pwindow, swindow = swindow, L = L, D = D)
  groups <- lapply(.param_groups(settings, names(settings)), function(g) {
    xs <- if (is.null(g$idx)) x else x[g$idx]
    g$x <- unique(xs)
    if (length(g$x) < length(xs)) {
      g$map <- match(xs, g$x)
    }
    g$exact <- g$swindow == 0
    if (!g$exact) {
      g$upper <- pmin(g$x + g$swindow, g$D)
    }
    g
  })
  pwindows <- vapply(groups, function(g) g$pwindow, numeric(1))
  unique_pw <- unique(pwindows)
  sets <- lapply(unique_pw, function(pw) {
    pts <- unlist(lapply(groups[pwindows == pw], function(g) {
      bounds <- c(g$L, g$D)
      c(g$x[!g$exact], g$upper, bounds[is.finite(bounds)])
    }))
    list(pwindow = pw, points = sort(unique(pts)))
  })
  for (i in seq_along(groups)) {
    g <- groups[[i]]
    g$set <- match(g$pwindow, unique_pw)
    pts <- sets[[g$set]]$points
    # An infinite bound is only among the points if a delay reaches it
    g$pos_L <- match(g$L, pts)
    g$pos_D <- match(g$D, pts)
    if (!g$exact) {
      g$lower <- match(g$x, pts)
      g$upper <- match(g$upper, pts)
    }
    groups[[i]] <- g
  }
  list(groups = groups, sets = sets)
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
    if (length(set$points) == 0L) {
      return(numeric(0))
    }
    cdf <- pcens_cdf(object, set$points, set$pwindow)
    # Some analytical methods return NaN at Inf
    cdf[set$points == -Inf] <- 0
    cdf[set$points == Inf] <- 1
    cdf
  })
  result <- numeric(n)
  for (g in grouped$groups) {
    rows <- if (is.null(g$idx)) seq_len(n) else g$idx
    result[rows] <- .pcens_pmf_group(object, g, cdfs[[g$set]])
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
