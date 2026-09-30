#' Build a reusable log-likelihood function for primary event censored data
#'
#' `r lifecycle::badge("experimental")`
#'
#' Sets up the log-likelihood of a set of observations once and returns a
#' function of the delay distribution parameters. It is intended for use
#' inside an optimiser or sampler, where the same data are evaluated many
#' times with different parameters.
#' The setup that [dprimarycensored()] repeats on every call is done once
#' here. Observations are grouped by their censoring and truncation
#' settings, `pdist` and `dprimary` are resolved and checked, and a `pcens`
#' object is built.
#' Each call then only updates the delay parameters and evaluates
#' [pcens_pmf()] once per group.
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
#'   with [check_pdist()] the first time the function is called, and again
#'   whenever it is called with a different set of parameter names. Set to
#'   `FALSE` to skip both.
#'
#' @details
#' The returned function is called as `ll(...)` with the delay distribution
#' parameters as named arguments, for example `ll(shape = 2, rate = 1)`.
#' Each call starts from the fixed parameters given in `...` here, so
#' parameters from an earlier call are not carried over.
#' Parameter names are checked against `pdist` whenever they differ from the
#' previous call, so a misspelt parameter raises an error and is never
#' silently ignored.
#'
#' Within each group of observations that share settings, the likelihood is
#' evaluated once at each unique value of `x` and copied to the matching
#' rows. This makes the cost depend on the number of unique delays rather
#' than the number of observations, which helps most with daily data.
#' Values are identical to those from `log(dprimarycensored())` called on
#' each observation.
#'
#' Errors from `pdist`, for example for parameters outside its support, are
#' not caught. Wrap the call in `tryCatch()` if an optimiser should see a
#' missing value instead. Zero probabilities are returned as `-Inf`.
#'
#' The construction checks that every observation satisfies `L <= x < D`.
#' A single message is given when any secondary interval extends past `D`
#' and is clipped at `D`. No message is given when the returned function is
#' called.
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
#'   200, rlnorm, meanlog = 1.3, sdlog = 0.5,
#'   pwindow = pw, swindow = 1, D = 20
#' )
#' ll <- pcens_loglik_fn(
#'   floor(x), plnorm, pwindow = pw, swindow = 1, D = 20
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
    check = TRUE) {
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
  # `.check = FALSE` would otherwise hide a misspelt name.
  state <- new.env(parent = emptyenv())
  state$checked <- FALSE
  state$names <- NULL
  function(...) {
    nms <- names(list(...))
    if (!state$checked || !identical(nms, state$names)) {
      obj <- update(base, ...)
      if (isTRUE(check)) {
        do.call(check_pdist, c(list(obj$pdist, D = max_D), obj$args))
      }
      state$names <- nms
      state$checked <- TRUE
    } else {
      obj <- update(base, ..., .check = FALSE)
    }
    log(.pcens_pmf_groups(obj, groups, n))
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
#' the unique values of `x` are kept with a map back to the rows, so that
#' [pcens_pmf()] is evaluated once per unique delay.
#'
#' @param x Numeric vector of delays.
#'
#' @param pwindow,swindow,L,D Per-observation settings, each of length 1
#'   or `length(x)`.
#'
#' @return A list with one element per group. Each element is a list with
#'   `idx`, the rows of the group or `NULL` if the group is every row, `x`,
#'   the unique delays of the group, `map`, the position of each row of the
#'   group in `x` or `NULL` if `x` has no repeats, the group's `pwindow`,
#'   `swindow`, `L` and `D`, and `clip`, whether any secondary interval of
#'   the group extends past `D`. Empty if `x` is empty.
#'
#' @keywords internal
.pcens_row_groups <- function(x, pwindow, swindow, L, D) {
  n <- length(x)
  if (n == 0L) {
    return(list())
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
  lapply(rows, function(idx) {
    first <- if (is.null(idx)) 1L else idx[[1L]]
    setting <- function(s) s[[if (length(s) == 1L) 1L else first]]
    xs <- if (is.null(idx)) x else x[idx]
    ux <- unique(xs)
    map <- NULL
    if (length(ux) < length(xs)) {
      map <- match(xs, ux)
    }
    g_swindow <- setting(swindow)
    g_D <- setting(D)
    list(
      idx = idx,
      x = ux,
      map = map,
      pwindow = setting(pwindow),
      swindow = g_swindow,
      L = setting(L),
      D = g_D,
      clip = is.finite(g_D) && any(ux + g_swindow > g_D)
    )
  })
}

#' Evaluate the primary event censored PMF for grouped observations
#'
#' @param object A `pcens` object.
#'
#' @param groups Groups of observations as made by [.pcens_row_groups()].
#'
#' @param n Number of observations.
#'
#' @return Numeric vector of length `n` with the PMF (or density, where
#'   `swindow = 0`) of each observation, in row order. Not on the log scale.
#'
#' @keywords internal
.pcens_pmf_groups <- function(object, groups, n) {
  if (length(groups) == 1L && is.null(groups[[1L]]$idx)) {
    return(.pcens_pmf_group(object, groups[[1L]]))
  }
  result <- numeric(n)
  for (g in groups) {
    result[g$idx] <- .pcens_pmf_group(object, g)
  }
  result
}

#' Evaluate the primary event censored PMF for one group of observations
#'
#' @inheritParams .pcens_pmf_groups
#'
#' @param group One element of the list made by [.pcens_row_groups()].
#'
#' @return Numeric vector with one value per row of the group.
#'
#' @keywords internal
.pcens_pmf_group <- function(object, group) {
  # The clipping notice from pcens_pmf() would be repeated on every call
  if (group$clip) {
    pmf <- suppressMessages(pcens_pmf(
      object, group$x, group$pwindow,
      swindow = group$swindow, L = group$L, D = group$D
    ))
  } else {
    pmf <- pcens_pmf(
      object, group$x, group$pwindow,
      swindow = group$swindow, L = group$L, D = group$D
    )
  }
  if (is.null(group$map)) pmf else pmf[group$map]
}
