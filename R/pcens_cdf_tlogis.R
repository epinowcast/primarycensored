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
#' With the window \eqn{[0, w]} and the logistic CDF
#' \eqn{L(z) = 1 / (1 + e^{-(z - m) / s)})}, location \eqn{m} and scale
#' \eqn{s}, the window density is \eqn{L'(z) / D_L} with
#' \eqn{D_L = L(w) - L(0)}. The primary event censored CDF at \eqn{q} is
#' \deqn{
#' F_L(q) = F(q - w) + \frac{1}{D_L}
#'   \int_{q - w}^{q} \{L(q - u) - L(0)\} f(u) du,
#' }
#' for a delay with density \eqn{f} and CDF \eqn{F}.
#' The integral is a sum of the finite horizon transforms \eqn{T_f(\xi; \tau)}
#' of [tilt_transform] at the tilts \eqn{\xi = \pm n / s}, found by expanding
#' \eqn{L} as a geometric series.
#' With \eqn{p = q - u} the primary event time, \eqn{L(p) = 1 / (1 + x)} with
#' \eqn{x = e^{-(p - m) / s}}.
#' For \eqn{p > m}, \eqn{x < 1} and \eqn{L(p) = \sum_{n \ge 0} (-x)^n}, which
#' needs the positive tilts \eqn{\xi = n / s}.
#' For \eqn{p < m}, \eqn{L(p) = \sum_{n \ge 1} (-1)^{n - 1} y^n} with
#' \eqn{y = 1 / x}, which needs the negative tilts \eqn{\xi = -n / s}.
#' The window is split at \eqn{p = m}, at \eqn{u^\star = q - m}, when the
#' location is inside it. With \eqn{\Delta T(\xi; a, b) = T_f(\xi; b) -
#' T_f(\xi; a)}, \eqn{\Delta F(a, b) = F(b) - F(a)}, \eqn{a = q - w} and
#' \eqn{b = q} the integral \eqn{\Phi} is
#' * for \eqn{m < 0}, where the whole window is above the location,
#'   \deqn{\Phi = \sum_{n \ge 1} (-1)^{n - 1} e^{n m / s}
#'   \{\Delta F(a, b) - e^{-n q / s} \Delta T(n / s; a, b)\}.}
#' * for \eqn{0 \le m < w}, the sum of
#'   \deqn{\Phi_P = \sum_{n \ge 0} (-1)^n e^{-n (q - m) / s}
#'   \Delta T(n / s; a, u^\star) - L(0) \Delta F(a, u^\star)}
#'   over the part of the window above the location and
#'   \deqn{\Phi_N = \sum_{n \ge 1} (-1)^{n - 1} e^{-n m / s}
#'   \{e^{n q / s} \Delta T(-n / s; u^\star, b) - \Delta F(u^\star, b)\}}
#'   over the part below it. For \eqn{m \ge w} only \eqn{\Phi_N} is needed
#'   with \eqn{u^\star = a}, and for \eqn{m = 0} only \eqn{\Phi_P}.
#'
#' and \eqn{F_L(q) = F(q - w) + \Phi / D_L}.
#' The terms that hold \eqn{L(0)} in place of the constant \eqn{1} of the
#' series are written as differences, so that nothing cancels when the window
#' is far from the location and \eqn{D_L} is small.
#' For delays on the non-negative reals \eqn{F(x) = T_f(\xi; x) = 0} for
#' \eqn{x \le 0}, and the normal delay has full support.
#' Each transform depends on one endpoint, \eqn{q} or \eqn{q - w} or
#' \eqn{u^\star}, so each endpoint is evaluated once and reused when several
#' `q` share one, as for the integer delays of [pcens_pmf()].
#' The split points \eqn{q - m} are shared when the location and window are
#' integers and are separate for each `q` otherwise.
#'
#' **Truncation and its error bound.** The geometric series have a ratio
#' that approaches 1 near the location, so the partial sums converge slowly
#' there. Each series \eqn{\sum (-x)^n c_n} is instead summed with weights
#' \eqn{W_n}, which are 1 for \eqn{n < n_0} and taper to 0 at
#' \eqn{n_0 + M} as the upper tail of a Binomial(\eqn{M}, 1/2). This is the
#' average of the partial sums \eqn{S_{n_0}, \ldots, S_{n_0 + M}} with
#' binomial weights, the Euler transform applied after \eqn{n_0} terms. For
#' each \eqn{x} the error is exactly
#' \eqn{x^{n_0} ((1 - x) / 2)^M / (1 + x)}. The terms of every series are
#' integrals against a positive measure of powers of \eqn{x} with \eqn{x} in
#' a known range, so the total error is at most the largest value of that
#' error over the range times the mass of the delay in the window. The rule
#' picks the smallest \eqn{n_0 + M}, among plain partial sums
#' (\eqn{M = 0}) and two tapers, for which that error is below
#' \eqn{10^{-10} D_L}. The error of \eqn{F_L(q)} is then below
#' \eqn{10^{-10}} of the mass of the delay in the window. This takes about
#' 6 terms for \eqn{x \le 0.01} and about 23 for ranges that reach 1.
#' If no rule with at most 64 terms exists the numerical method is used.
#'
#' **Admissibility.** The series need the positive tilts up to
#' \eqn{(n_0 + M) / s} and the negative tilts down to \eqn{-(n_0 + M) / s}.
#' The exponential and gamma forms need \eqn{\lambda - \xi > 0}, so that the
#' tilted delay distribution exists, which fails for large positive tilts
#' unless the rate \eqn{\lambda} exceeds the largest tilt. The negative tilts
#' always exist. The method falls back to [pcens_cdf.default()] for every `q`
#' where the largest positive tilt is not available or there is no
#' truncation rule. For gamma and exponential delays this is the case for a
#' location at or before the window end unless the rate is large. The normal
#' form has no restriction. The fallback is also used for a window that is not
#' a single positive finite number.
#'
#' **Precision.** The CDF agrees with a reference integral to a relative
#' difference of about 1e-9 or better, including in the tails.
#' The terms carry a relative error of about the machine precision divided
#' by \eqn{D_L} in the amount of cancellation of \eqn{\Phi / D_L}, which is
#' large only for a scale much larger than the window where the window
#' approaches the uniform one.
#'
#' **Extending.** A new delay distribution is supported by the methods of
#' [tilt_transform] and a `pcens_cdf` method for its class that calls
#' `.pcens_cdf_tlogis()`. A new window reuses the same transforms with its
#' own tilts and coefficients, see [pcens_cdf_exptilt].
#'
#' @family pcens
#'
#' @inherit pcens_cdf return
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

# Multiples of the scale from the centre of a narrow window density at which
# the numerical method breaks the integral. The density falls by about e^-x
# at x scales, so the mass beyond 30 scales is below 1e-13.
.tlogis_break_multiples <- c(-30, -12, -6, -3, -1.5, 0, 1.5, 3, 6, 12, 30)

# Fraction of the terms of a series kept with weight 1 by the taper
.tlogis_taper <- 0.35

# Delays on the non-negative reals with q below this multiple of the scale
# use the small delay form
.tlogis_small <- 1e-4

#' Number of terms and taper of a series for the truncated logistic window
#'
#' The series \eqn{\sum (-x)^n c_n} is summed with the weights of
#' [.tlogis_weights()], 1 for the first `n0` terms and then tapering to 0
#' over `M` more. The error of the weighted sum of the geometric series for
#' one \eqn{x} is \eqn{x^{n_0} ((1 - x) / 2)^M / (1 + x)}, see
#' [pcens_cdf_tlogis]. Its largest value over \eqn{[x_{lo}, x_{hi}]} is at the
#' point \eqn{n_0 / (n_0 + M)} clipped to the range, and at \eqn{x_{lo}} for
#' \eqn{n_0 = 0}. This searches the plain partial sum and two tapers for the
#' fewest terms with that bound below the tolerance.
#'
#' @param log_lo,log_hi Log of the smallest and largest \eqn{x} of the series.
#'
#' @param log_tol Log of the tolerance on the bound of the error.
#'
#' @param max_terms The most terms allowed.
#'
#' @return An integer vector with elements `n0` and `M`, or `NULL` if no rule
#'   with at most `max_terms` terms meets the tolerance.
#'
#' @keywords internal
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
#' The weights of the average of the partial sums \eqn{S_{n_0}, \ldots,
#' S_{n_0 + M}} with Binomial(\eqn{M}, 1/2) weights. They are 1 for
#' \eqn{n < n_0} and the probability that the Binomial is more than
#' \eqn{n - n_0} for the rest.
#'
#' @param n0 Number of terms with weight 1, at least 0.
#'
#' @param M Number of terms that taper to 0, at least 0.
#'
#' @return A numeric vector of length `n0 + M`.
#'
#' @keywords internal
.tlogis_weights <- function(n0, M) {
  c(
    rep(1, n0),
    stats::pbinom(seq_len(M) - 1L, M, 0.5, lower.tail = FALSE)
  )
}

#' Plan the series for a truncated logistic primary
#'
#' Finds which series the primary event censored CDF needs for the window and
#' location and scale of `object`, with the terms and weights of each, and
#' whether the delay has the transforms they need.
#'
#' @param object A `pcens` object with a [dtlogis()] primary.
#'
#' @param pwindow Primary event window.
#'
#' @return `NULL` if the analytical solution does not apply, which is where
#'   the window is not a single positive finite number, no truncation rule
#'   exists or a tilt of the series is not available for the delay, see
#'   [.pcens_tilt_available()]. Otherwise a list with the `location`,
#'   `scale`, `pwindow`, the log mass `log_mass` of the window, the point
#'   `split` of the window where the expansion changes and up to two
#'   series, `pos` for the positive tilts and `neg` for the negative tilts.
#'   Each series is a list with its `form`, the first tilt index `first`, the
#'   number of terms `n0` and `M` and the `weights` of the terms.
#'
#' @keywords internal
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
#' @return `NULL` if a series that is needed has no truncation rule. Otherwise
#'   a list with the series `pos` and `neg`, each `NULL` where the window does
#'   not need it, see [.tlogis_plan()].
#'
#' @keywords internal
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
#' The series run to the tilts \eqn{\pm (\text{first} + n_0 + M - 1) / s}. The
#' exponential and gamma delays need the positive one to be below their rate.
#'
#' @param object A `pcens` object.
#'
#' @param series List with the series `pos` and `neg`, see
#'   [.tlogis_plan_series()].
#'
#' @param scale Scale of the primary.
#'
#' @return `TRUE` if [.pcens_tilt_available()] is `TRUE` for the largest tilt
#'   of each series, otherwise `FALSE`.
#'
#' @keywords internal
.tlogis_tilts_available <- function(object, series, scale) {
  largest <- function(x) (x$first + x$n0 + x$M - 1L) / scale
  (is.null(series$pos) ||
    .pcens_tilt_available(object, largest(series$pos))) &&
    (is.null(series$neg) ||
      .pcens_tilt_available(object, -largest(series$neg)))
}

#' Location and scale of the truncated logistic primary of a pcens object
#'
#' @inheritParams .tlogis_plan
#'
#' @return A list with `location` and `scale`, defaulting as in [dtlogis()].
#'
#' @keywords internal
.tlogis_primary_args <- function(object) {
  m <- object$primary_args$location
  s <- object$primary_args$scale
  if (is.null(m)) {
    m <- 0
  }
  if (is.null(s)) {
    s <- 1
  }
  if (!is.numeric(m) || length(m) != 1L || !is.finite(m)) {
    stop(
      "location must be a single finite number for the truncated logistic ",
      "primary distribution",
      call. = FALSE
    )
  }
  if (!is.numeric(s) || length(s) != 1L || !is.finite(s) || s <= 0) {
    stop(
      "scale must be a single positive finite number for the truncated ",
      "logistic primary distribution",
      call. = FALSE
    )
  }
  list(location = m, scale = s)
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
#' @keywords internal
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
#' @param plan The plan from [.tlogis_plan()].
#'
#' @return Vector of CDFs, clamped to \[0, 1\].
#'
#' @keywords internal
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
#' For delays on the non-negative reals with `q` below the window and
#' \eqn{q / s} below 1e-4 only the primary event times \eqn{p \le q}
#' contribute and \eqn{L(p) - L(0)} is expanded in \eqn{p}. With
#' \eqn{G_k(q) = \int_0^q (q - u)^k f(u) du},
#' \deqn{F_L(q) = \{L'(0) G_1(q) + L''(0) G_2(q) / 2\} / D_L + O((q / s)^2),}
#' where \eqn{L'(0) = \sigma_0 / s}, \eqn{L''(0) = L'(0) (1 - 2 L(0)) / s}
#' and \eqn{\sigma_0 = L(0) (1 - L(0))}. The direct form cancels here, losing
#' about 1e-16 s / q of relative precision, while the truncation error of this
#' form is about \eqn{(q / s)^2 / 6}. Both are below 1e-8 at the threshold.
#'
#' @inheritParams .tlogis_cdf_finite
#'
#' @return Vector of log CDFs.
#'
#' @keywords internal
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
#' Evaluates the transforms at the unique endpoints, with all the endpoints
#' at or below the lower end of the support sharing one entry, so each is
#' computed once however many `q` use it.
#'
#' @param object A `pcens` object.
#'
#' @param t Numeric vector of finite endpoints, which may repeat.
#'
#' @param xi Numeric vector of tilts.
#'
#' @return A list with the sorted unique endpoints `t`, and the matrices
#'   `lower` and `upper` of the log transforms over the lower and the upper
#'   part of the support, with a row for each endpoint and a column for each
#'   tilt.
#'
#' @keywords internal
.tlogis_endpoint_terms <- function(object, t, xi) {
  lower <- .pcens_tilt_lower(object)
  if (is.finite(lower)) {
    t <- pmax(t, lower)
  }
  t <- sort.int(unique(t))
  if (.pcens_tilt_vectorised(object)) {
    # One call for all tilts, each at every endpoint
    grid_t <- rep(t, times = length(xi))
    grid_xi <- rep(xi, each = length(t))
    return(list(
      t = t,
      lower = matrix(
        .pcens_tilt_transform(object, grid_t, grid_xi), nrow = length(t)
      ),
      upper = matrix(
        .pcens_tilt_transform(object, grid_t, grid_xi, upper = TRUE),
        nrow = length(t)
      )
    ))
  }
  lower_terms <- matrix(0, length(t), length(xi))
  upper_terms <- matrix(0, length(t), length(xi))
  for (j in seq_along(xi)) {
    lower_terms[, j] <- .pcens_tilt_transform(object, t, xi[[j]])
    upper_terms[, j] <- .pcens_tilt_transform(object, t, xi[[j]], upper = TRUE)
  }
  list(t = t, lower = lower_terms, upper = upper_terms)
}

#' Log of the difference of a transform between two endpoints
#'
#' Evaluates \eqn{\log(T_f(\xi; hi) - T_f(\xi; lo))} for each tilt of
#' `terms` and each pair of `lo` and `hi`, from the lower or the upper
#' transforms, whichever loses less precision, see `.exptilt_tail_diff()`.
#'
#' @param terms Output of [.tlogis_endpoint_terms()].
#'
#' @param lo,hi Numeric vectors of the ends of the intervals, `lo <= hi`,
#'   that are among the endpoints of `terms`.
#'
#' @param lower Lower end of the support of the delay.
#'
#' @return A matrix with a row for each interval and a column for each tilt.
#'
#' @keywords internal
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
#' Evaluates \eqn{\sum_n (-1)^n W_n c_n} from the log of the terms \eqn{c_n}.
#' The terms are scaled by `scale` so nothing underflows or overflows.
#'
#' @param log_terms Matrix of the log of the terms, a row for each `q`.
#'
#' @param weights Numeric vector of the weights, one for each column.
#'
#' @param scale Numeric vector of the log scale of each row, finite.
#'
#' @return Numeric vector of the sums divided by `exp(scale)`.
#'
#' @keywords internal
.tlogis_weighted_sum <- function(log_terms, weights, scale) {
  signs <- rep_len(c(1, -1), length(weights))
  drop(exp(log_terms - scale) %*% (signs * weights))
}

#' Truncated logistic log CDF
#'
#' The direct form, from the series of the plan. Each series gives a part of
#' the integral \eqn{\Phi} on its own scale, see [pcens_cdf_tlogis].
#'
#' @inheritParams .tlogis_cdf_finite
#'
#' @return Vector of log CDFs.
#'
#' @keywords internal
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
  .log_sum_exp(
    .pcens_tilt_transform(object, a, 0), log_phi - plan$log_mass
  )
}

#' Part of the truncated logistic integral for a location before the window
#'
#' For \eqn{m < 0} the whole window is above the location and
#' \eqn{\Phi = \sum_{n \ge 1} (-1)^{n - 1} e^{n m / s}
#' \{\Delta F(a, b) - e^{-n q / s} \Delta T(n / s; a, b)\}}, with the
#' terms weighted by the plan.
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
#' @keywords internal
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
#' For \eqn{0 \le m < w}, the part of the window with primary event time
#' above the location, between `a` and \eqn{u^\star = q - m}, is
#' \eqn{\Phi_P = \sum_{n \ge 0} (-1)^n e^{-n (q - m) / s}
#' \Delta T(n / s; a, u^\star) - L(0) \Delta F(a, u^\star)}, with the terms
#' weighted by the plan.
#'
#' @inheritParams .tlogis_part_before_window
#'
#' @param u Numeric vector of the split points \eqn{u^\star = q - m}.
#'
#' @inherit .tlogis_part_before_window return
#'
#' @keywords internal
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
#' For \eqn{m > 0}, the part of the window with primary event time below the
#' location, between \eqn{u^\star = q - \min(m, w)} and `q`, is
#' \eqn{\Phi_N = \sum_{n \ge 1} (-1)^{n - 1} e^{-n m / s}
#' \{e^{n q / s} \Delta T(-n / s; u^\star, b) - \Delta F(u^\star, b)\}},
#' with the terms weighted by the plan.
#'
#' @inheritParams .tlogis_part_above
#'
#' @inherit .tlogis_part_before_window return
#'
#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
.tlogis_row_max <- function(x) {
  out <- apply(x, 1L, max)
  out[!is.finite(out)] <- 0
  out
}
