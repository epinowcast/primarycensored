#' Methods for delays with a truncated logistic primary
#'
#' Analytical primary event censored CDFs for exponential, gamma and normal
#' delay distributions with a truncated logistic primary event window, the
#' [dtlogis()] primary distribution with `location` and `scale`. They honour
#' `use_numeric`, and use the numerical method of [pcens_cdf.default()] when
#' no closed form applies.
#'
#' @inheritParams pcens_cdf
#'
#' @details
#' With the window \eqn{[0, w]} and the logistic CDF \eqn{L} with location
#' \eqn{m} and scale \eqn{s}, the window density is \eqn{L'(z) / D_L} with
#' \eqn{D_L = L(w) - L(0)}. For a delay with density \eqn{f} and CDF
#' \eqn{F}, the primary event censored CDF at \eqn{q} is
#' \deqn{F_L(q) = F(q - w) + \frac{1}{D_L}
#'   \int_{q - w}^{q} \{L(q - u) - L(0)\} f(u) du.}
#' Expanding \eqn{L} as a geometric series writes the integral as sums of the
#' transforms \eqn{T_f(\xi; \tau)} of [tilt_transform] at the tilts
#' \eqn{\xi = \pm n / s}. The positive tilts are for the part of the window
#' above the location and the negative tilts for the part below it, so the
#' window is split at the location when it is inside it. Each transform
#' depends on one endpoint, so each endpoint is evaluated once and reused
#' when several `q` share it, as for the integer delays of [pcens_pmf()].
#' Terms that hold \eqn{L(0)} are written as differences so nothing cancels
#' when the window is far from the location.
#'
#' **Truncation.** Near the location the series converge slowly. Each series
#' is summed with binomial taper weights (the Euler transform after \eqn{n_0}
#' terms), which bounds the error for each term by
#' \eqn{x^{n_0} ((1 - x) / 2)^M / (1 + x)}. The fewest terms, at most 64,
#' with a total error below \eqn{10^{-10} D_L} are used. The numerical method
#' is used where there is no such rule.
#'
#' **Admissibility.** The exponential and gamma forms need the largest
#' positive tilt of the series to be below the rate, which holds for a
#' location after the window or a large rate. Otherwise, and for a window that
#' is not a single positive finite number, [pcens_cdf.default()] is used. The
#' normal form has no restriction.
#'
#' **Precision.** The CDF agrees with a reference integral to a relative
#' difference of about 1e-8 in general.
#' In the far tails the error grows with the cancellation in the integral,
#' to about 1e-7 for a scale much larger than the window.
#'
#' A new delay distribution is supported by the methods of [tilt_transform]
#' and a `pcens_cdf` method for its class that calls `.pcens_cdf_tlogis()`.
#'
#' @inherit pcens_cdf return
#'
#' @concept pcens
#'
#' @name pcens_cdf_tlogis
#'
#' @examples
#' # Normal delay, for example a difference of event times, with a logistic
#' # window that saturates inside the primary event window
#' pnorm_obj <- new_pcens(
#'   pdist = pnorm, dprimary = dtlogis,
#'   primary_args = list(location = 0.5, scale = 0.2), mean = 3, sd = 2
#' )
#' pcens_cdf(pnorm_obj, q = c(-1, 3, 8), pwindow = 2)
#'
#' # Gamma delay with the midpoint after the window, where the analytic
#' # solution applies for any rate
#' pgamma_obj <- new_pcens(
#'   pdist = pgamma, dprimary = dtlogis,
#'   primary_args = list(location = 3, scale = 0.5), shape = 2, rate = 1
#' )
#' pcens_cdf(pgamma_obj, q = c(1, 4, 8), pwindow = 2)
NULL

#' @rdname pcens_cdf_tlogis
#' @export
pcens_cdf.pcens_pexp_dtlogis <- function(
  object,
  q,
  pwindow,
  use_numeric = FALSE
) {
  .pcens_cdf_tlogis(object, q, pwindow, use_numeric)
}

#' @rdname pcens_cdf_tlogis
#' @export
pcens_cdf.pcens_pgamma_dtlogis <- function(
  object,
  q,
  pwindow,
  use_numeric = FALSE
) {
  .pcens_cdf_tlogis(object, q, pwindow, use_numeric)
}

#' @rdname pcens_cdf_tlogis
#' @export
pcens_cdf.pcens_pnorm_dtlogis <- function(
  object,
  q,
  pwindow,
  use_numeric = FALSE
) {
  .pcens_cdf_tlogis(object, q, pwindow, use_numeric)
}

# Target for the truncation error, relative to the mass of the window
.tlogis_tol <- 1e-10

# Most terms of a series before the numerical method is used
.tlogis_max_terms <- 64L

# Fraction of the terms of a series kept with weight 1 by the taper
.tlogis_taper <- 0.35

# Delays on the non-negative reals with q below this multiple of the scale
# use the small delay form
.tlogis_small <- 1e-4

#' Number of terms and taper of a series for the truncated logistic window
#'
#' Searches the plain partial sum and two tapers for the fewest terms whose
#' bound on the error, see [pcens_cdf_tlogis], is below the tolerance for
#' every x in the range.
#'
#' @param log_lo,log_hi Log of the smallest and largest x of the series.
#'
#' @param log_tol Log of the tolerance.
#'
#' @param max_terms The most terms allowed.
#'
#' @return An integer vector `c(n0, M)`, or `NULL` if no rule meets the
#'   tolerance.
#'
#' @noRd
.tlogis_series_terms <- function(log_lo, log_hi, log_tol,
                                 max_terms = .tlogis_max_terms) {
  if (!is.finite(log_tol)) {
    return(NULL)
  }
  # Log of the largest error over the range for n0 direct terms and M
  # tapered, for each of the candidates
  log_bound <- function(n0, m_taper) {
    out <- rep(0, length(n0))
    plain <- m_taper == 0
    out[plain] <- n0[plain] * log_hi
    none <- !plain & n0 == 0
    out[none] <- m_taper[none] * (.log1m_exp(log_lo) - log(2))
    both <- !plain & !none
    log_x <- pmin(
      pmax(log(n0[both] / (n0[both] + m_taper[both])), log_lo), log_hi
    )
    out[both] <- m_taper[both] * (.log1m_exp(log_x) - log(2)) +
      n0[both] * log_x
    out
  }
  k <- seq_len(max_terms)
  # The plain partial sum, the taper and the taper from the first term
  candidates <- rbind(k, floor(.tlogis_taper * k), 0, deparse.level = 0)
  ok <- rbind(
    log_bound(candidates[1L, ], k - candidates[1L, ]),
    log_bound(candidates[2L, ], k - candidates[2L, ]),
    log_bound(candidates[3L, ], k - candidates[3L, ])
  ) <= log_tol
  first <- which(ok, arr.ind = TRUE)
  if (nrow(first) == 0L) {
    return(NULL)
  }
  # The fewest terms, and for those the first candidate that meets the bound
  first <- first[order(first[, 2L], first[, 1L])[1L], ]
  n0 <- unname(candidates[first[[1L]], first[[2L]]])
  c(n0 = n0, M = first[[2L]] - n0)
}

#' Weights of the terms of a series of the truncated logistic window
#'
#' 1 for the first `n0` terms, then the upper tail of a Binomial(`M`, 1/2).
#'
#' @param n0,M Number of terms with weight 1 and number that taper to 0.
#'
#' @return A numeric vector of length `n0 + M`.
#'
#' @noRd
.tlogis_weights <- function(n0, M) {
  c(
    rep(1, n0),
    stats::pbinom(seq_len(M) - 1L, M, 0.5, lower.tail = FALSE)
  )
}

#' Plan the series for a truncated logistic primary
#'
#' @param object A `pcens` object with a [dtlogis()] primary.
#'
#' @param pwindow Primary event window.
#'
#' @return `NULL` if the analytical solution does not apply, otherwise a list
#'   with the primary parameters, the log mass of the window, the split
#'   point of the window and the series `pos` and `neg` of positive and
#'   negative tilts, each with its `form`, first tilt index `first`, `n0`, `M`
#'   and `weights`.
#'
#' @noRd
.tlogis_plan <- function(object, pwindow) {
  if (length(pwindow) != 1L || !is.finite(pwindow) || pwindow <= 0) {
    return(NULL)
  }
  primary <- .tlogis_primary_args(object)
  log_mass <- .tlogis_log_diff(
    0, pwindow, primary$location, primary$scale
  )
  series <- .tlogis_plan_series(
    primary$location, primary$scale, pwindow, log_mass
  )
  if (is.null(series) ||
    !.tlogis_tilts_available(object, series, primary$scale)) {
    return(NULL)
  }
  list(
    location = primary$location, scale = primary$scale, pwindow = pwindow,
    log_mass = log_mass, split = min(max(primary$location, 0), pwindow),
    pos = series$pos, neg = series$neg
  )
}

#' Series of the truncated logistic primary for a location, scale and window
#'
#' @param location,scale Location and scale of the primary.
#'
#' @param pwindow Primary event window, positive.
#'
#' @param log_mass Log of the mass of the window.
#'
#' @return `NULL` if a series that is needed has no truncation rule, otherwise
#'   a list with the series `pos` and `neg`, each `NULL` if not needed.
#'
#' @noRd
.tlogis_plan_series <- function(location, scale, pwindow, log_mass) {
  log_tol <- log(.tlogis_tol) + log_mass
  series <- function(form, first, log_lo, log_hi) {
    rule <- .tlogis_series_terms(log_lo, log_hi, log_tol)
    if (is.null(rule)) {
      return(NULL)
    }
    list(
      form = form, first = first, n0 = rule[["n0"]], M = rule[["M"]],
      weights = .tlogis_weights(rule[["n0"]], rule[["M"]])
    )
  }
  pos <- NULL
  neg <- NULL
  if (location < 0) {
    # The whole window is above the location. x = exp(-(p - m) / s) is at
    # most exp(m / s) and at least exp((m - w) / s).
    pos <- series(
      "A", 1L, (location - pwindow) / scale, location / scale
    )
  } else {
    if (location < pwindow) {
      pos <- series("P", 0L, -(pwindow - location) / scale, 0)
    }
    if (location > 0) {
      # y = exp((p - m) / s) for p in [0, min(m, w)]
      neg <- series(
        "N", 1L, -location / scale, min(0, (pwindow - location) / scale)
      )
    }
  }
  if ((is.null(pos) && location < pwindow) || (is.null(neg) && location > 0)) {
    return(NULL)
  }
  list(pos = pos, neg = neg)
}

#' Check the delay has the tilts of the truncated logistic series
#'
#' @param object A `pcens` object.
#'
#' @param series List with the series `pos` and `neg`.
#'
#' @param scale Scale of the primary.
#'
#' @return `TRUE` if the largest tilt of each series is available.
#'
#' @noRd
.tlogis_tilts_available <- function(object, series, scale) {
  largest <- function(x) (x$first + x$n0 + x$M - 1L) / scale
  (is.null(series$pos) ||
    .pcens_tilt_available(object, largest(series$pos))) &&
    (is.null(series$neg) ||
      .pcens_tilt_available(object, -largest(series$neg)))
}

#' Primary event censored CDF for a truncated logistic primary
#'
#' Shared implementation of the [pcens_cdf_tlogis] methods. It dispatches
#' on the delay class of `object` through the generics of [tilt_transform].
#'
#' @inheritParams pcens_cdf
#'
#' @return Vector of computed primary event censored CDFs.
#'
#' @noRd
.pcens_cdf_tlogis <- function(object, q, pwindow, use_numeric = FALSE) {
  if (isTRUE(use_numeric)) {
    return(pcens_cdf.default(object, q, pwindow, use_numeric))
  }
  plan <- .tlogis_plan(object, pwindow)
  if (is.null(plan)) {
    return(pcens_cdf.default(object, q, pwindow, use_numeric))
  }
  result <- rep(NA_real_, length(q))
  result[!is.na(q) & q == Inf] <- 1
  result[!is.na(q) & q == -Inf] <- 0
  finite <- which(is.finite(q))
  if (length(finite) > 0L) {
    result[finite] <- .tlogis_cdf_finite(object, q[finite], plan)
  }
  result
}

#' Truncated logistic CDF at finite points
#'
#' @param object A `pcens` object.
#'
#' @param q Numeric vector of finite quantiles.
#'
#' @param plan The plan from `.tlogis_plan()`.
#'
#' @return Vector of CDFs, clamped to \[0, 1\].
#'
#' @noRd
.tlogis_cdf_finite <- function(object, q, plan) {
  lower <- .pcens_tilt_lower(object)
  positive <- is.finite(lower)
  # Below the support of the delay no mass has arrived
  active <- !positive | q > lower
  # The direct form cancels for a delay on the non-negative reals with q
  # small relative to the scale, where the window is close to its derivative
  # at 0
  small <- positive & active & q < plan$pwindow &
    q < .tlogis_small * plan$scale
  direct <- active & !small
  log_cdf <- rep(-Inf, length(q))
  if (any(small)) {
    log_cdf[small] <- .tlogis_lcdf_small_delay(object, q[small], plan)
  }
  if (any(direct)) {
    log_cdf[direct] <- .tlogis_lcdf(object, q[direct], plan)
  }
  pmin(1, exp(log_cdf))
}

#' Small delay form of the truncated logistic log CDF
#'
#' For delays on the non-negative reals with `q / scale` below 1e-4 the
#' direct form cancels, so `L(p) - L(0)` is expanded to second order in `p`.
#' The truncation error is about `(q / scale)^2 / 6`.
#'
#' @inheritParams .tlogis_cdf_finite
#'
#' @return Vector of log CDFs.
#'
#' @noRd
.tlogis_lcdf_small_delay <- function(object, q, plan) {
  s <- plan$scale
  moments <- .pcens_tilt_moments(object, q)
  log_l0 <- stats::plogis(-plan$location / s, log.p = TRUE)
  log_l0c <- stats::plogis(-plan$location / s, lower.tail = FALSE, log.p = TRUE)
  log_l1 <- log_l0 + log_l0c - log(s)
  ratio <- (1 - 2 * exp(log_l0)) / s
  out <- moments[, 1L] + log_l1 - plan$log_mass +
    log1p(0.5 * ratio * exp(moments[, 2L] - moments[, 1L]))
  out[!is.finite(moments[, 1L])] <- -Inf
  out
}

#' Transforms of a delay at the endpoints of a truncated logistic series
#'
#' Evaluates the transforms at the unique endpoints, so each is computed once
#' however many `q` use it.
#'
#' @param object A `pcens` object.
#'
#' @param t Numeric vector of finite endpoints, which may repeat.
#'
#' @param xi Numeric vector of tilts.
#'
#' @return A list with the sorted unique endpoints `t` and the matrices
#'   `lower` and `upper` of the log transforms, a row for each endpoint and a
#'   column for each tilt.
#'
#' @noRd
.tlogis_endpoint_terms <- function(object, t, xi) {
  lower <- .pcens_tilt_lower(object)
  if (is.finite(lower)) {
    t <- pmax(t, lower)
  }
  t <- sort.int(unique(t))
  # One call for all tilts, each at every endpoint
  grid_t <- rep(t, times = length(xi))
  grid_xi <- rep(xi, each = length(t))
  list(
    t = t,
    lower = matrix(
      .pcens_tilt_transform(object, grid_t, grid_xi),
      nrow = length(t)
    ),
    upper = matrix(
      .pcens_tilt_transform(object, grid_t, grid_xi, upper = TRUE),
      nrow = length(t)
    )
  )
}

#' Log of the difference of a transform between two endpoints
#'
#' Uses the lower or the upper transforms, whichever loses less precision.
#'
#' @param terms Output of `.tlogis_endpoint_terms()`.
#'
#' @param lo,hi Numeric vectors of interval ends among the endpoints of
#'   `terms`.
#'
#' @param lower Lower end of the support of the delay.
#'
#' @return A matrix with a row for each interval and a column for each tilt.
#'
#' @noRd
.tlogis_endpoint_diff <- function(terms, lo, hi, lower) {
  if (is.finite(lower)) {
    lo <- pmax(lo, lower)
    hi <- pmax(hi, lower)
  }
  i_lo <- match(lo, terms$t)
  i_hi <- match(hi, terms$t)
  n_tilts <- ncol(terms$lower)
  out <- .exptilt_tail_diff(
    as.vector(terms$lower[i_hi, , drop = FALSE]),
    as.vector(terms$lower[i_lo, , drop = FALSE]),
    as.vector(terms$upper[i_hi, , drop = FALSE]),
    as.vector(terms$upper[i_lo, , drop = FALSE])
  )
  matrix(out, nrow = length(lo), ncol = n_tilts)
}

#' Weighted alternating sum of the terms of a series
#'
#' @param log_terms Matrix of the log of the terms, a row for each `q`.
#'
#' @param weights Numeric vector of the weights, one for each column.
#'
#' @param scale Numeric vector of the log scale of each row.
#'
#' @return Numeric vector of the sums divided by `exp(scale)`.
#'
#' @noRd
.tlogis_weighted_sum <- function(log_terms, weights, scale) {
  signs <- rep_len(c(1, -1), length(weights))
  drop(exp(log_terms - scale) %*% (signs * weights))
}

#' Truncated logistic log CDF, direct form
#'
#' @inheritParams .tlogis_cdf_finite
#'
#' @return Vector of log CDFs.
#'
#' @noRd
.tlogis_lcdf <- function(object, q, plan) {
  lower <- .pcens_tilt_lower(object)
  a <- q - plan$pwindow
  u <- q - plan$split
  # F at the ends of the window and the split point
  f_terms <- .tlogis_endpoint_terms(object, c(a, u, q), 0)
  diff_f <- function(lo, hi) {
    .tlogis_endpoint_diff(f_terms, lo, hi, lower)[, 1L]
  }
  parts <- list()
  if (!is.null(plan$pos) && plan$pos$form == "A") {
    parts$pos <- .tlogis_part_before_window(object, q, a, plan, diff_f)
  } else if (!is.null(plan$pos)) {
    parts$pos <- .tlogis_part_above(object, q, a, u, plan, diff_f)
  }
  if (!is.null(plan$neg)) {
    parts$neg <- .tlogis_part_below(object, q, u, plan, diff_f)
  }
  log_phi <- rep(-Inf, length(q))
  for (part in parts) {
    log_phi <- .log_sum_exp(log_phi, part$scale + log(pmax(part$sum, 0)))
  }
  # F(a) is in `f_terms`, at the same ends as they were evaluated
  log_f_a <- f_terms$lower[
    match(if (is.finite(lower)) pmax(a, lower) else a, f_terms$t), 1L
  ]
  .log_sum_exp(log_f_a, log_phi - plan$log_mass)
}

#' Part of the truncated logistic integral for a location before the window
#'
#' @inheritParams .tlogis_cdf_finite
#'
#' @param a Numeric vector of the lower ends of the integral, `q - pwindow`.
#'
#' @param diff_f Function of the ends of an interval giving the log of the
#'   difference of the delay CDF between them.
#'
#' @return A list with the `sum` of the part and the log `scale` it is in.
#'
#' @noRd
.tlogis_part_before_window <- function(object, q, a, plan, diff_f) {
  lower <- .pcens_tilt_lower(object)
  s <- plan$scale
  k <- plan$pos$first + seq_along(plan$pos$weights) - 1L
  endpoint_terms <- .tlogis_endpoint_terms(object, c(a, q), k / s)
  log_dt <- .tlogis_endpoint_diff(endpoint_terms, a, q, lower)
  log_c <- .log_diff_exp(
    rep(diff_f(a, q), times = length(k)),
    as.vector(-outer(q, k) / s + log_dt)
  )
  log_c <- matrix(log_c, nrow = length(q)) +
    rep(k * plan$location / s, each = length(q))
  .tlogis_weighted_part(log_c, plan$pos$weights)
}

#' Part of the truncated logistic integral above the location
#'
#' @inheritParams .tlogis_part_before_window
#'
#' @param u Numeric vector of the split points `q - location`.
#'
#' @inherit .tlogis_part_before_window return
#'
#' @noRd
.tlogis_part_above <- function(object, q, a, u, plan, diff_f) {
  lower <- .pcens_tilt_lower(object)
  m <- plan$location
  s <- plan$scale
  k <- plan$pos$first + seq_along(plan$pos$weights) - 1L
  endpoint_terms <- .tlogis_endpoint_terms(object, c(a, u), k / s)
  log_c <- .tlogis_endpoint_diff(endpoint_terms, a, u, lower) -
    outer(q - m, k) / s
  # The constant term of the series, L(0) times the mass of the delay
  log_constant <- stats::plogis(-m / s, log.p = TRUE) + diff_f(a, u)
  scale <- pmax(.tlogis_row_max(log_c), log_constant)
  list(
    sum = .tlogis_weighted_sum(log_c, plan$pos$weights, scale) -
      exp(log_constant - scale),
    scale = scale
  )
}

#' Part of the truncated logistic integral below the location
#'
#' @inheritParams .tlogis_part_above
#'
#' @inherit .tlogis_part_before_window return
#'
#' @noRd
.tlogis_part_below <- function(object, q, u, plan, diff_f) {
  lower <- .pcens_tilt_lower(object)
  s <- plan$scale
  k <- plan$neg$first + seq_along(plan$neg$weights) - 1L
  endpoint_terms <- .tlogis_endpoint_terms(object, c(u, q), -k / s)
  log_dt <- .tlogis_endpoint_diff(endpoint_terms, u, q, lower)
  log_c <- .log_diff_exp(
    as.vector(outer(q, k) / s + log_dt),
    rep(diff_f(u, q), times = length(k))
  )
  log_c <- matrix(log_c, nrow = length(q)) -
    rep(k * plan$location / s, each = length(q))
  .tlogis_weighted_part(log_c, plan$neg$weights)
}

#' Weighted sum of the terms of a series with its scale
#'
#' @param log_c Matrix of the log of the terms, a row for each `q`.
#'
#' @param weights Numeric vector of the weights, one for each column.
#'
#' @inherit .tlogis_part_before_window return
#'
#' @noRd
.tlogis_weighted_part <- function(log_c, weights) {
  scale <- .tlogis_row_max(log_c)
  list(sum = .tlogis_weighted_sum(log_c, weights, scale), scale = scale)
}

#' Row maxima of a matrix of log terms
#'
#' @param x Numeric matrix.
#'
#' @return Numeric vector of the maximum of each row, with rows of only
#'   `-Inf` given a finite value so they scale to zero.
#'
#' @noRd
.tlogis_row_max <- function(x) {
  out <- apply(x, 1L, max)
  out[!is.finite(out)] <- 0
  out
}
