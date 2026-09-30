#' Methods for delays with a truncated Gumbel primary
#'
#' Analytical primary event censored CDFs for exponential, gamma and normal
#' delay distributions with a truncated Gumbel primary event window, the
#' [dtgumbel()] primary distribution. They honour `use_numeric`, and use the
#' numerical method of [pcens_cdf.default()] when no accurate closed form
#' applies.
#'
#' @inheritParams pcens_cdf
#'
#' @details
#' Write \eqn{G(z) = \exp[-\exp\{-(z - \mu) / \beta\}]} and
#' \eqn{D = G(w) - G(0)} for a window of width \eqn{w}, and let \eqn{F} be
#' the delay CDF and \eqn{T_f(\xi; \tau)} the transform of [tilt_transform].
#' With \eqn{c(q) = \exp\{-(q - \mu) / \beta\}} the reversed window
#' \eqn{G(q - u) = \exp\{-c(q) e^{u / \beta}\}} is a power series in
#' \eqn{e^{u / \beta}}, which gives
#' \deqn{
#' F_G(q) = F(q - w) + \frac{1}{D} \Big[
#'   \sum_{n = 0}^\infty \frac{\{-c(q)\}^n}{n!}
#'   \{T_f(n / \beta; q) - T_f(n / \beta; q - w)\}
#'   - G(0) \{F(q) - F(q - w)\} \Big].
#' }
#' This is the solution of the paper, and was checked numerically.
#' The \eqn{n = 0} term is \eqn{F(q) - F(q - w)}, so the bracket is computed
#' as
#' \deqn{B(q) = (1 - G(0)) \{F(q) - F(q - w)\} +
#'   \sum_{n \ge 1} \frac{\{-c(q)\}^n}{n!} \Delta_n(q),}
#' with \eqn{\Delta_n(q) = T_f(n / \beta; q) - T_f(n / \beta; q - w)}. Then
#' \eqn{1 - G(0)} is evaluated as `-expm1(-exp(mu / beta))`, which keeps
#' precision when \eqn{\mu} is far below the window and \eqn{D} is small.
#' The transforms at each endpoint, \eqn{q} and \eqn{q - w}, are evaluated
#' once at every tilt and reused when several `q` share an endpoint. Every
#' \eqn{\Delta_n} is taken between the lower tail transforms or the upper
#' tail transforms, whichever loses less precision, as for
#' [pcens_cdf_exptilt].
#'
#' **Sums of alternating terms.** With \eqn{s(z) = e^{(\mu - z) / \beta}},
#' the \eqn{n}th term is
#' \eqn{a_n = \int s(z)^n f(q - z) dz / n!}, at most
#' \eqn{s_0^n / n!} times the delay mass in the window, where
#' \eqn{s_0 = e^{\mu / \beta}} is the largest value of \eqn{s} on the window.
#' The terms do not depend on \eqn{q} beyond this bound because
#' \eqn{c(q)} grows as the transforms \eqn{T_f(n / \beta; q)} shrink, and
#' they are evaluated on the log scale so they do not overflow. The sum
#' alternates, so a term that is larger than the result loses precision.
#' The series is truncated at the first \eqn{N} for which the next term
#' relative to \eqn{\max(s_0, 1)} is below 1e-18, see `.gumbel_n_terms()`.
#' The odd and even terms are accumulated separately on the log scale, and
#' the difference is taken last.
#'
#' **Accuracy region.** The relative error is estimated for each `q` from
#' the terms, as the machine precision times the ratio of the sum of the
#' absolute terms to the result, times one plus the largest magnitude of a
#' log term, plus the truncation bound relative to the result. This is
#' `.gumbel_error_bound()`. Where it is above
#' \eqn{10^{-9}} the method uses the numerical method of
#' [pcens_cdf.default()] for that `q`, so no value is silently inaccurate.
#' Two things cause it. The first is a large \eqn{s_0 = e^{\mu / \beta}},
#' where the terms reach \eqn{e^{s_0}} and the result is of order 1, so about
#' \eqn{s_0 \gtrsim 10} loses precision as \eqn{e^{s_0}} times the rounding
#' error of the log terms. The second is a window that is narrow relative to
#' the scale, \eqn{w \ll \beta}, where the bracket is a small difference of
#' large terms and the loss is about \eqn{\beta / w}.
#' Also \eqn{s_0 > 15} is never used (see `.gumbel_max_s0`) as the number of
#' terms grows with it.
#' In the grid of `mu` in -0.5, 0, 0.5, 1, 1.5, `beta` in 0.1, 0.2, 1 and
#' `pwindow` in 1, 2, the series is used for the normal delay when
#' \eqn{\mu / \beta} is at most about 2, and agrees with numerical
#' integration to a relative difference of about 1e-9 or better.
#'
#' **Admissibility.** The terms need the transform at the tilts
#' \eqn{n / \beta} for \eqn{n = 1, \ldots, N}. The exponential and gamma forms
#' need the tilted delay to exist, \eqn{\lambda > N / \beta} for rate
#' \eqn{\lambda}. Otherwise the method falls back to [pcens_cdf.default()] for
#' every `q`. For delays of a day or more the rate is of order 1, and the
#' bound is only met for very short delays, so in practice these delays use
#' the numerical method unless \eqn{\mu / \beta} is very negative or the
#' scale is large. The normal delay has no restriction.
#' A window with \eqn{\mu / \beta} very negative is close to the exponentially
#' tilted window with \eqn{\rho = -1 / \beta}, see [pcens_cdf_exptilt], which
#' only needs \eqn{\lambda + \rho > 0}.
#'
#' **Extending.** A new delay distribution is supported by the methods of
#' [tilt_transform] and a `pcens_cdf` method for the class that calls
#' `.pcens_cdf_gumbel()`. The tilts are \eqn{n / \beta}, all positive, so a
#' delay needs a transform for positive tilts.
#'
#' @family pcens
#'
#' @inherit pcens_cdf return
#'
#' @name pcens_cdf_gumbel
#'
#' @examples
#' # Normal delay, for example a difference of event times
#' pnorm_obj <- new_pcens(
#'   pdist = pnorm, dprimary = dtgumbel,
#'   primary_args = list(mu = 0.5, beta = 0.5), mean = 3, sd = 2
#' )
#' pcens_cdf(pnorm_obj, q = c(-1, 3, 8), pwindow = 2)
#'
#' # Exponential delay with a short mean, rate 40
#' pexp_obj <- new_pcens(
#'   pdist = pexp, dprimary = dtgumbel,
#'   primary_args = list(mu = -1, beta = 0.5), rate = 40
#' )
#' pcens_cdf(pexp_obj, q = c(0.01, 0.1, 1), pwindow = 1)
NULL

#' @rdname pcens_cdf_gumbel
#' @export
pcens_cdf.pcens_pexp_dtgumbel <- function(
  object,
  q,
  pwindow,
  use_numeric = FALSE
) {
  .pcens_cdf_gumbel(object, q, pwindow, use_numeric)
}

#' @rdname pcens_cdf_gumbel
#' @export
pcens_cdf.pcens_pgamma_dtgumbel <- function(
  object,
  q,
  pwindow,
  use_numeric = FALSE
) {
  .pcens_cdf_gumbel(object, q, pwindow, use_numeric)
}

#' @rdname pcens_cdf_gumbel
#' @export
pcens_cdf.pcens_pnorm_dtgumbel <- function(
  object,
  q,
  pwindow,
  use_numeric = FALSE
) {
  .pcens_cdf_gumbel(object, q, pwindow, use_numeric)
}

# The series is not used above this value of exp(mu / beta). The terms are
# as large as exp(s0) and the number of terms grows with s0 (about 60 at 12),
# so beyond it the result needs more than 1e-9 of the rounding error.
.gumbel_max_s0 <- 15

# Largest estimated relative error for which the series is used
.gumbel_tol <- 1e-9

# Truncation threshold of the series, relative to max(s0, 1)
.gumbel_log_trunc <- log(1e-18)

#' Number of terms of the Gumbel series
#'
#' The smallest \eqn{N} such that the bound \eqn{s_0^n / n!} of the terms
#' after \eqn{N} is below `1e-18 * max(s0, 1)`, where
#' \eqn{s_0 = e^{\mu / \beta}}.
#'
#' @param log_s0 \eqn{\mu / \beta}, the log of the largest value of
#'   \eqn{s} on the window.
#'
#' @return The number of terms, at least 2.
#'
#' @keywords internal
.gumbel_n_terms <- function(log_s0) {
  n <- seq_len(200L)
  log_bound <- n * log_s0 - lgamma(n + 1) - max(log_s0, 0)
  which(log_bound < .gumbel_log_trunc)[[1L]]
}

#' Primary arguments of a truncated Gumbel pcens object
#'
#' @inheritParams .pcens_cdf_gumbel
#'
#' @return A list with the location `mu` and the scale `beta`, checked to be
#'   single valid numbers.
#'
#' @keywords internal
.gumbel_primary_args <- function(object) {
  primary <- object$primary_args[c("mu", "beta")]
  if (is.null(primary$mu) || is.null(primary$beta)) {
    stop(
      "mu and beta parameters are required for the truncated Gumbel ",
      "primary distribution",
      call. = FALSE
    )
  }
  .check_tgumbel(0, 1, primary$mu, primary$beta)
  if (length(primary$mu) != 1L || length(primary$beta) != 1L) {
    stop(
      "mu and beta must be single numbers for the truncated Gumbel ",
      "primary distribution",
      call. = FALSE
    )
  }
  primary
}

#' Test whether the Gumbel series applies for a delay and primary
#'
#' The series needs \eqn{\mu / \beta} below `.gumbel_max_s0` and the
#' transform of the delay at the largest tilt \eqn{N / \beta}, with
#' \eqn{N} from `.gumbel_n_terms()`. This depends on the parameters only,
#' whether it is accurate at a `q` is decided by `.gumbel_error_bound()`.
#'
#' @inheritParams .pcens_cdf_gumbel
#'
#' @param mu,beta Location and scale of the truncated Gumbel primary.
#'
#' @return `TRUE` if the series can be used.
#'
#' @keywords internal
.gumbel_available <- function(object, mu, beta) {
  mu / beta <= log(.gumbel_max_s0) &&
    .pcens_tilt_available(object, .gumbel_n_terms(mu / beta) / beta)
}

#' Primary event censored CDF for a truncated Gumbel primary
#'
#' Shared implementation of the [pcens_cdf_gumbel] methods. It dispatches
#' on the delay class of `object` through the generics of [tilt_transform].
#'
#' @inheritParams pcens_cdf
#'
#' @return Vector of computed primary event censored CDFs.
#'
#' @keywords internal
.pcens_cdf_gumbel <- function(object, q, pwindow, use_numeric = FALSE) {
  if (isTRUE(use_numeric)) {
    return(pcens_cdf.default(object, q, pwindow, use_numeric))
  }
  primary <- .gumbel_primary_args(object)
  mu <- primary$mu
  scale <- primary$beta
  # The closed forms are for a single window and need a bounded number of
  # terms, and the tilted delay at the largest tilt
  if (length(pwindow) != 1L || !is.finite(pwindow) || pwindow <= 0 ||
    !.gumbel_available(object, mu, scale)) {
    return(.gumbel_numeric(object, q, pwindow, mu))
  }
  n_terms <- .gumbel_n_terms(mu / scale)

  result <- rep(NA_real_, length(q))
  result[!is.na(q) & q == Inf] <- 1
  result[!is.na(q) & q == -Inf] <- 0
  finite <- which(is.finite(q))
  if (length(finite) > 0L) {
    result[finite] <- .gumbel_cdf_finite(
      object, q[finite], pwindow, mu, scale, n_terms
    )
  }
  result
}

#' Numerical primary event censored CDF for a truncated Gumbel primary
#'
#' This is [pcens_cdf.default()]. The window density can go from zero to
#' its peak within a fraction of the window, for example for a large
#' location with a small scale, and the default integration can then fail
#' with a roundoff error. For a point where it does, the integral is taken
#' again with more subdivisions and a break at the location `mu`, where the
#' window density changes most.
#'
#' @inheritParams pcens_cdf
#'
#' @param mu Location of the truncated Gumbel primary.
#'
#' @return Vector of CDFs, in \[0, 1\].
#'
#' @keywords internal
.gumbel_numeric <- function(object, q, pwindow, mu) {
  vapply(
    q,
    function(d) {
      tryCatch(
        pcens_cdf.default(object, d, pwindow, FALSE),
        error = function(e) {
          integrand <- function(p) {
            do.call(object$pdist, c(list(q = d - p), object$args)) *
              do.call(
                object$dprimary,
                c(list(x = p, min = 0, max = pwindow), object$dprimary_args)
              )
          }
          breaks <- sort(unique(c(
            0, pwindow,
            if (d > 0 && d < pwindow) d,
            if (mu > 0 && mu < pwindow) mu
          )))
          value <- sum(vapply(
            seq_len(length(breaks) - 1L),
            function(i) {
              stats::integrate(
                integrand, breaks[i], breaks[i + 1L],
                rel.tol = 1e-9, subdivisions = 1000L, stop.on.error = FALSE
              )$value
            },
            numeric(1)
          ))
          min(1, max(0, value))
        }
      )
    },
    numeric(1)
  )
}

#' Truncated Gumbel CDF at finite points
#'
#' Evaluates the series at the unique endpoints, and uses the numerical
#' method for the points where the series loses accuracy.
#'
#' @inheritParams pcens_cdf
#'
#' @param mu,beta Location and scale of the truncated Gumbel primary.
#'
#' @param n_terms Number of terms of the series, see `.gumbel_n_terms()`.
#'
#' @return Vector of CDFs, clamped to \[0, 1\].
#'
#' @keywords internal
.gumbel_cdf_finite <- function(object, q, pwindow, mu, beta, n_terms) {
  lower <- .pcens_tilt_lower(object)
  positive <- is.finite(lower)
  log_cdf <- rep(-Inf, length(q))
  # Below the support of the delay no mass has arrived
  active <- which(!positive | q > lower)
  numeric_needed <- logical(length(q))
  if (length(active) > 0L) {
    fit <- .gumbel_lcdf(
      object, q[active], pwindow, mu, beta, n_terms, lower
    )
    log_cdf[active] <- fit$log_cdf
    numeric_needed[active] <- fit$error > .gumbel_tol | is.na(fit$error)
  }
  result <- pmin(1, exp(log_cdf))
  if (any(numeric_needed)) {
    result[numeric_needed] <- .gumbel_numeric(
      object, q[numeric_needed], pwindow, mu
    )
  }
  result
}

#' Log CDF for a truncated Gumbel primary and its estimated error
#'
#' @inheritParams .gumbel_cdf_finite
#'
#' @param lower Lower end of the support of the delay.
#'
#' @return A list with `log_cdf`, and `error`, the estimated relative error,
#'   see `.gumbel_error_bound()`.
#'
#' @keywords internal
.gumbel_lcdf <- function(object, q, pwindow, mu, beta, n_terms, lower) {
  endpoints <- .exptilt_endpoints(q, pwindow, lower)
  index_q <- .exptilt_index(q, endpoints, lower)
  index_y <- .exptilt_index(q - pwindow, endpoints, lower)
  # Transforms at each endpoint for the tilts n / beta, n = 0, ..., n_terms,
  # as log T_f over the lower and upper parts of the support. Every term
  # depends on one endpoint only.
  tilts <- seq.int(0L, n_terms) / beta
  lower_terms <- vapply(
    tilts, function(xi) .pcens_tilt_transform(object, endpoints, xi),
    numeric(length(endpoints))
  )
  upper_terms <- vapply(
    tilts, function(xi) {
      .pcens_tilt_transform(object, endpoints, xi, upper = TRUE)
    },
    numeric(length(endpoints))
  )
  if (length(endpoints) == 1L) {
    lower_terms <- matrix(lower_terms, nrow = 1L)
    upper_terms <- matrix(upper_terms, nrow = 1L)
  }
  # Log of Delta_n(q) = T_f(n / beta; q) - T_f(n / beta; q - w), n = 0, ...
  log_delta <- vapply(
    seq_along(tilts),
    function(i) {
      .exptilt_tail_diff(
        lower_terms[index_q, i], lower_terms[index_y, i],
        upper_terms[index_q, i], upper_terms[index_y, i]
      )
    },
    numeric(length(q))
  )
  if (length(q) == 1L) {
    log_delta <- matrix(log_delta, nrow = 1L)
  }
  log_f_y <- lower_terms[index_y, 1L]
  # The scale of the log transforms at both ends, from either tail
  finite_abs <- function(x) {
    x <- abs(x)
    x[!is.finite(x)] <- 0
    x
  }
  scale_t <- pmax(
    finite_abs(lower_terms[index_q, , drop = FALSE]),
    finite_abs(upper_terms[index_q, , drop = FALSE]),
    finite_abs(lower_terms[index_y, , drop = FALSE]),
    finite_abs(upper_terms[index_y, , drop = FALSE])
  )
  .gumbel_combine(
    log_delta, log_f_y, scale_t, q, pwindow, mu, beta
  )
}

#' Combine the Gumbel series terms into a log CDF
#'
#' @param log_delta Matrix of the log of \eqn{\Delta_n(q)}, one row per `q`
#'   and a column for each tilt \eqn{n = 0, \ldots, N}.
#'
#' @param log_f_y Log of the delay CDF at `q - pwindow`.
#'
#' @param scale_t Matrix of the scale of the log transforms at `q` and
#'   `q - pwindow` for each tilt, for the rounding error of the terms.
#'
#' @inheritParams .gumbel_cdf_finite
#'
#' @return A list with `log_cdf` and `error`.
#'
#' @keywords internal
.gumbel_combine <- function(log_delta, log_f_y, scale_t, q, pwindow, mu,
                            beta) {
  n_terms <- ncol(log_delta) - 1L
  n <- seq_len(n_terms)
  log_s0 <- mu / beta
  # log c(q) and the log of the terms a_n = c^n / n! Delta_n, n >= 1
  log_c <- -(q - mu) / beta
  log_a <- sweep(log_delta[, n + 1L, drop = FALSE], 2L, lgamma(n + 1), "-") +
    outer(log_c, n)
  # Positive terms are the even n and (1 - G(0)) Delta_0, negative the odd
  is_odd <- n %% 2L == 1L
  log_pos <- .log_sum_exp_rows(cbind(
    .log1m_exp(-exp(log_s0)) + log_delta[, 1L],
    log_a[, !is_odd, drop = FALSE]
  ))
  log_neg <- .log_sum_exp_rows(log_a[, is_odd, drop = FALSE])
  log_bracket <- .log_diff_exp(log_pos, log_neg)
  # The normalisation is the difference of G at the ends of the window. On
  # the log scale it is minus s_w plus the log of one minus exp(-delta),
  # where delta is the difference of s across the window, which is s_w times
  # expm1 of the window over beta
  log_s_w <- -(pwindow - mu) / beta
  log_d <- -exp(log_s_w) +
    .log1m_exp(-exp(log_s_w + .log_expm1(pwindow / beta)))
  log_cdf <- .log_sum_exp(log_f_y, log_bracket - log_d)
  # Mass in the window is zero, so all terms vanish together
  no_mass <- is.infinite(log_delta[, 1L]) & log_delta[, 1L] < 0
  log_cdf[no_mass] <- log_f_y[no_mass]
  error <- .gumbel_error_bound(
    log_pos, log_neg, log_bracket, log_delta[, 1L], log_c, log_s0, n_terms,
    scale_t, n
  )
  error[no_mass] <- 0
  list(log_cdf = log_cdf, error = error)
}

#' Estimated relative error of the Gumbel series
#'
#' The rounding error of the sum is the machine precision times the sum of
#' the absolute terms, which is amplified by the rounding error in the log
#' terms, about the machine precision times the magnitude of the logs that
#' are added, and the truncation error is bounded by
#' \eqn{s_0^{N + 1} / (N + 1)!} times the delay mass in the window.
#' Both are relative to the bracket.
#'
#' @param log_pos,log_neg,log_bracket Log of the positive part, the negative
#'   part and their difference.
#'
#' @param log_mass Log of the delay mass in the window, \eqn{\Delta_0}.
#'
#' @param log_c Log of \eqn{c(q)}.
#'
#' @param log_s0 \eqn{\mu / \beta}.
#'
#' @param n_terms Number of terms.
#'
#' @param scale_t Matrix of the scale of the log transforms.
#'
#' @param n Integer vector `1:n_terms`.
#'
#' @return Vector of estimated relative errors, `Inf` where the bracket is
#'   not positive.
#'
#' @keywords internal
.gumbel_error_bound <- function(log_pos, log_neg, log_bracket, log_mass,
                                log_c, log_s0, n_terms, scale_t, n) {
  # Magnitude of the log arithmetic behind each term
  magnitude <- 1 + apply(
    abs(outer(log_c, n)) + scale_t[, n + 1L, drop = FALSE], 1L, max
  )
  log_rounding <- log(.Machine$double.eps) + log(magnitude) +
    .log_sum_exp(log_pos, log_neg) - log_bracket
  log_trunc <- log_mass + (n_terms + 1) * log_s0 - lgamma(n_terms + 2) -
    log1p(-min(exp(log_s0) / (n_terms + 2), 0.5)) - log_bracket
  error <- exp(log_rounding) + exp(log_trunc)
  error[!is.finite(log_bracket)] <- Inf
  error
}

#' Row-wise log of a sum of exponentials
#'
#' @param x Numeric matrix on the log scale.
#'
#' @return Vector with the log of the sum of each row, `-Inf` for a row of
#'   `-Inf`.
#'
#' @keywords internal
.log_sum_exp_rows <- function(x) {
  top <- apply(x, 1L, max)
  shifted <- x - ifelse(is.finite(top), top, 0)
  top + log(rowSums(exp(shifted)))
}
