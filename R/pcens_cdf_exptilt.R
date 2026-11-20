#' Methods for delays with an exponentially tilted primary
#'
#' Analytical primary event censored CDFs for exponential, gamma, normal,
#' lognormal, Weibull and generalised gamma delay distributions with an
#' exponentially tilted primary event window, the [dexpgrowth()] primary
#' distribution with `r` equal to the tilt \eqn{\rho}.
#' They honour `use_numeric`, and use the numerical method of
#' [pcens_cdf.default()] when no closed form applies.
#'
#' @inheritParams pcens_cdf
#'
#' @details
#' With window width \eqn{w}, delay CDF \eqn{F} and
#' \eqn{J(x) = \int e^{-\rho u} f(u) du} up to \eqn{x}, the CDF at \eqn{q} is
#' \deqn{
#' F_\rho(q) = F(q - w) + \frac{e^{\rho q} \{J(q) - J(q - w)\} -
#'   \{F(q) - F(q - w)\}}{e^{\rho w} - 1}.
#' }
#' Every term depends on one endpoint, \eqn{q} or \eqn{q - w}, so each
#' endpoint is evaluated once and shared between neighbouring `q`.
#' Terms are on the log scale, and each difference is taken between the lower
#' or the upper tail terms, whichever loses less precision.
#'
#' The exponential and gamma forms need \eqn{\lambda + \rho > 0} for rate
#' \eqn{\lambda}, so that the tilted delay distribution exists.
#' Otherwise the numerical method is used.
#' The numerical method is less accurate in the lower tail of a gamma with
#' shape below 1, with a relative error of about 1e-5 in R and 1e-4 in Stan.
#' The lognormal transform has no closed form.
#' It is evaluated by quadrature for \eqn{\rho > 0} and by a series of
#' partial moments for \eqn{\rho < 0}.
#' The numerical method is used where the series needs too many terms or is
#' slower than numerical integration.
#'
#' The Weibull and the generalised gamma with power \eqn{a}, scale
#' \eqn{\theta} and \eqn{k} (\eqn{k = 1} for the Weibull) have the transform
#' \eqn{J(x) = \sum_n (-\rho \theta)^n / n!\,
#' \Gamma(k + n / a) / \Gamma(k)\, P(k + n / a, (x / \theta)^a)},
#' for the regularised incomplete gamma function \eqn{P}.
#' It exists for every tilt, and needs \eqn{Q > 0} for
#' `flexsurv::pgengamma()`.
#' For \eqn{\rho > 0} its terms alternate in sign and cancel.
#' Where they cancel badly, or \eqn{\rho q \ge 400}, the numerical method is
#' used.
#'
#' The direct form cancels as \eqn{\rho \to 0}.
#' With \eqn{G_k(t) = \int (t - u)^k f(u) du} up to \eqn{t}, two series in
#' \eqn{\rho} truncated after the second order replace it.
#' * For \eqn{|\rho| w} below \eqn{10^{-2}}, or \eqn{10^{-3}} for delays on
#'   the reals, and \eqn{|\rho| (|q| + w) < 0.1}, the uniform window limit
#'   with its corrections,
#'   \eqn{\{G_1(q) - G_1(q - w)\} / w +
#'   \rho \{G_2(q) - w G_1(q) - G_2(q - w) - w G_1(q - w)\} / (2 w) +
#'   \rho^2 \{G_3(q) - G_3(q - w) - 3 w [G_2(q) + G_2(q - w)] / 2 +
#'   w^2 [G_1(q) - G_1(q - w)] / 2\} / (6 w)}.
#'   Its relative error is below \eqn{(|\rho| w)^3 / 160}.
#' * Otherwise, for delays on the non-negative reals with \eqn{q < w} and
#'   \eqn{|\rho| q < 5 \times 10^{-3}},
#'   \eqn{\rho \{G_1(q) + \rho G_2(q) / 2 + \rho^2 G_3(q) / 6\} /
#'   (e^{\rho w} - 1)}.
#'   Its relative error is below \eqn{(|\rho| q)^3 / 24}.
#'
#' Both errors are below \eqn{10^{-8}} at these limits.
#'
#' The CDF agrees with numerical integration to a relative 1e-6 in the
#' package tests.
#'
#' @inherit pcens_cdf return
#'
#' @name pcens_cdf_exptilt
#'
#' @family pcens
#'
#' @examples
#' # Exponential delay with a growing primary event process
#' pexp_obj <- new_pcens(
#'   pdist = pexp, dprimary = dexpgrowth,
#'   primary_args = list(r = 0.3), rate = 0.5
#' )
#' pcens_cdf(pexp_obj, q = c(1, 4, 8), pwindow = 2)
#'
#' # Gamma delay with a declining process
#' pgamma_obj <- new_pcens(
#'   pdist = pgamma, dprimary = dexpgrowth,
#'   primary_args = list(r = -0.3), shape = 2, rate = 1
#' )
#' pcens_cdf(pgamma_obj, q = c(1, 4, 8), pwindow = 2)
#'
#' # Normal delay, for example a difference of event times
#' pnorm_obj <- new_pcens(
#'   pdist = pnorm, dprimary = dexpgrowth,
#'   primary_args = list(r = 0.1), mean = 3, sd = 2
#' )
#' pcens_cdf(pnorm_obj, q = c(-1, 3, 8), pwindow = 2)
#'
#' # Weibull delay
#' pweibull_obj <- new_pcens(
#'   pdist = pweibull, dprimary = dexpgrowth,
#'   primary_args = list(r = 0.2), shape = 1.5, scale = 5
#' )
#' pcens_cdf(pweibull_obj, q = c(1, 4, 8), pwindow = 2)
NULL

#' @rdname pcens_cdf_exptilt
#' @export
pcens_cdf.pcens_pexp_dexpgrowth <- function(
  object,
  q,
  pwindow,
  use_numeric = FALSE
) {
  .pcens_cdf_exptilt(object, q, pwindow, use_numeric)
}

#' @rdname pcens_cdf_exptilt
#' @export
pcens_cdf.pcens_pgamma_dexpgrowth <- function(
  object,
  q,
  pwindow,
  use_numeric = FALSE
) {
  .pcens_cdf_exptilt(object, q, pwindow, use_numeric)
}

#' @rdname pcens_cdf_exptilt
#' @export
pcens_cdf.pcens_pnorm_dexpgrowth <- function(
  object,
  q,
  pwindow,
  use_numeric = FALSE
) {
  .pcens_cdf_exptilt(object, q, pwindow, use_numeric)
}

#' @rdname pcens_cdf_exptilt
#' @export
pcens_cdf.pcens_plnorm_dexpgrowth <- function(
  object,
  q,
  pwindow,
  use_numeric = FALSE
) {
  .pcens_cdf_exptilt(
    object, q, pwindow, use_numeric,
    min_q = .lnorm_exptilt_min_q, min_xw = .lnorm_exptilt_min_xw
  )
}

#' @rdname pcens_cdf_exptilt
#' @export
pcens_cdf.pcens_pweibull_dexpgrowth <- function(
  object,
  q,
  pwindow,
  use_numeric = FALSE
) {
  .pcens_cdf_exptilt(object, q, pwindow, use_numeric)
}

#' @rdname pcens_cdf_exptilt
#' @export
# nolint start: object_length_linter.
pcens_cdf.pcens_pgengamma.orig_dexpgrowth <- function(
  object,
  q,
  pwindow,
  use_numeric = FALSE
) {
  .pcens_cdf_exptilt(object, q, pwindow, use_numeric)
}
# nolint end

#' @rdname pcens_cdf_exptilt
#' @export
pcens_cdf.pcens_pgengamma_dexpgrowth <- function(
  object,
  q,
  pwindow,
  use_numeric = FALSE
) {
  .pcens_cdf_exptilt(object, q, pwindow, use_numeric)
}

# Largest ratio of the terms of a series or of the direct form to its result
.exptilt_max_loss <- 1e5

#' Tilt of the exponentially tilted primary of a pcens object
#'
#' @inheritParams pcens_cdf
#'
#' @return The tilt `r`, a single finite number.
#'
#' @keywords internal
.exptilt_rho <- function(object) {
  rho <- object$primary_args$r
  if (is.null(rho)) {
    stop(
      "r parameter is required for the exponential growth primary ",
      "distribution",
      call. = FALSE
    )
  }
  if (!is.numeric(rho) || length(rho) != 1L || !is.finite(rho)) {
    stop(
      "r must be a single finite number for the exponential growth primary ",
      "distribution",
      call. = FALSE
    )
  }
  rho
}

#' Primary event censored CDF for an exponentially tilted primary
#'
#' Shared implementation of the [pcens_cdf_exptilt] methods.
#'
#' @inheritParams pcens_cdf
#'
#' @param min_q Number of quantiles below which the numerical method is
#'   used, for delays whose transform has a fixed cost. The default 0 always
#'   uses the closed forms.
#'
#' @param min_xw Largest \eqn{|\rho| w} at which `min_q` applies. The
#'   numerical method loses accuracy for a larger tilt times the window, so
#'   the closed forms are used there for any number of quantiles.
#'
#' @return Vector of CDFs.
#'
#' @keywords internal
.pcens_cdf_exptilt <- function(
  object,
  q,
  pwindow,
  use_numeric = FALSE,
  min_q = 0L,
  min_xw = Inf
) {
  if (isTRUE(use_numeric)) {
    return(pcens_cdf.default(object, q, pwindow, use_numeric))
  }
  rho <- .exptilt_rho(object)
  if (length(pwindow) != 1L || !is.finite(pwindow) || pwindow <= 0 ||
    !.pcens_tilt_available(object, -rho)) {
    return(pcens_cdf.default(object, q, pwindow, use_numeric))
  }

  result <- rep(NA_real_, length(q))
  result[!is.na(q) & q == Inf] <- 1
  result[!is.na(q) & q == -Inf] <- 0
  finite <- which(is.finite(q))
  if (length(finite) > 0L) {
    # Use the numerical method for the q where the transform does not fit.
    # With fewer than `min_q` quantiles and a small |rho| w, `fits` stays
    # FALSE and every point uses the numerical method.
    fits <- rep(FALSE, length(finite))
    if (length(finite) >= min_q || abs(rho) * pwindow > min_xw) {
      fits <- .pcens_tilt_fits(object, q[finite], -rho, pwindow)
    }
    result[finite[fits]] <- .exptilt_cdf_finite(
      object, q[finite[fits]], pwindow, rho
    )
    if (!all(fits)) {
      result[finite[!fits]] <- pcens_cdf.default(
        object, q[finite[!fits]], pwindow, use_numeric
      )
    }
  }
  result
}

#' Exponentially tilted CDF at finite points
#'
#' Chooses the form at each point as in the Stan function
#' `primarycensored_exptilt_form()`, see [pcens_cdf_exptilt].
#'
#' @inheritParams pcens_cdf
#'
#' @param rho The tilt, the `r` of the exponential growth primary.
#'
#' @return Vector of CDFs, clamped to \[0, 1\].
#'
#' @keywords internal
.exptilt_cdf_finite <- function(object, q, pwindow, rho) {
  if (length(q) == 0L) {
    return(numeric(0))
  }
  lower <- .pcens_tilt_lower(object)
  positive <- is.finite(lower)
  log_cdf <- rep(-Inf, length(q))
  active <- !positive | q > lower
  small_window <- rep(FALSE, length(q))
  small_delay <- small_window
  if (abs(rho) * pwindow < if (positive) 1e-2 else 1e-3) {
    small_window <- active & abs(rho) * (abs(q) + pwindow) < 0.1
  } else if (positive) {
    small_delay <- active & q < pwindow & abs(rho) * q < 5e-3
  }
  direct <- active & !small_window & !small_delay

  if (any(small_window)) {
    log_cdf[small_window] <- .exptilt_lcdf_small_window(
      object, q[small_window], pwindow, rho, lower
    )
  }
  if (any(small_delay)) {
    log_cdf[small_delay] <- .exptilt_lcdf_small_delay(
      object, q[small_delay], pwindow, rho
    )
  }
  if (any(direct)) {
    log_cdf[direct] <- .exptilt_lcdf_direct(
      object, q[direct], pwindow, rho, lower
    )
  }
  result <- pmin(1, exp(log_cdf))
  # NaN is where the closed form is not reliable
  numerical <- which(is.nan(result))
  if (length(numerical) > 0L) {
    result[numerical] <- pcens_cdf.default(object, q[numerical], pwindow)
  }
  result
}

#' Endpoints and the positions of `q` and `q - pwindow` among them
#'
#' Each term of the exponentially tilted CDF depends on one endpoint, so it
#' is computed once per unique endpoint. Endpoints at or below a finite
#' `lower` share one entry.
#'
#' @inheritParams .exptilt_cdf_finite
#'
#' @param lower Lower end of the support of the delay, 0 or `-Inf`.
#'
#' @return List with the unique `endpoints` and the positions `q` and `y` of
#'   `q` and `q - pwindow` among them.
#'
#' @keywords internal
.exptilt_layout <- function(q, pwindow, lower) {
  at_q <- pmax(q, lower)
  at_y <- pmax(q - pwindow, lower)
  endpoints <- unique(c(at_q, at_y))
  list(
    endpoints = endpoints,
    q = match(at_q, endpoints),
    y = match(at_y, endpoints)
  )
}

#' Log of `abs(exp(rho * pwindow) - 1)`
#'
#' @inheritParams .exptilt_cdf_finite
#'
#' @return The log of the normalising constant of the tilted window.
#'
#' @keywords internal
.exptilt_log_den <- function(rho, pwindow) {
  if (rho > 0) .log_diff_exp(rho * pwindow, 0) else .log1m_exp(rho * pwindow)
}

#' Direct form of the exponentially tilted log CDF
#'
#' The expression of [pcens_cdf_exptilt] on the log scale.
#' The result is `NaN` where the transform at either endpoint is, and where
#' the transform has no upper tail and the terms of the numerator exceed the
#' smaller of the CDF and its survival function by more than
#' `.exptilt_max_loss`.
#'
#' @inheritParams .exptilt_layout
#'
#' @param object A `pcens` object.
#'
#' @return Vector of log CDFs.
#'
#' @keywords internal
.exptilt_lcdf_direct <- function(object, q, pwindow, rho, lower) {
  pos <- .exptilt_layout(q, pwindow, lower)
  endpoints <- pos$endpoints
  j_pair <- .pcens_tilt_pair(object, endpoints, -rho)
  endpoint_terms <- cbind(
    .pcens_tilt_pair(object, endpoints, 0),
    j_pair
  )
  at_q <- endpoint_terms[pos$q, , drop = FALSE]
  at_y <- endpoint_terms[pos$y, , drop = FALSE]
  log_diff_f <- .exptilt_tail_diff(
    at_q[, 1], at_y[, 1], at_q[, 2], at_y[, 2]
  )
  log_diff_j <- .exptilt_tail_diff(
    at_q[, 3], at_y[, 3], at_q[, 4], at_y[, 4]
  )
  log_num <- if (rho > 0) {
    .log_diff_exp(rho * q + log_diff_j, log_diff_f)
  } else {
    .log_diff_exp(log_diff_f, rho * q + log_diff_j)
  }
  log_den <- .exptilt_log_den(rho, pwindow)
  out <- .log_sum_exp(at_y[, 1], log_num - log_den)
  j_loss <- attr(j_pair, "log_loss")
  if (is.null(j_loss)) {
    j_loss <- rep(0, length(endpoints))
  }
  log_size <- rho * q + at_q[, 3] + pmax(j_loss[pos$q], j_loss[pos$y])
  no_tail <- is.na(at_q[, 4]) | is.na(at_y[, 4])
  # The delay survival function is a lower bound on that of the CDF
  log_survival <- pmax(at_q[, 2], .log1m_exp(pmin(out, 0)))
  unreliable <- is.na(at_q[, 3]) | is.na(at_y[, 3]) |
    (no_tail & is.finite(out) &
      log_size - log_den - pmin(out, log_survival) > log(.exptilt_max_loss))
  out[which(unreliable)] <- NaN
  out
}

#' Log of a difference between two points from either tail
#'
#' With `L(t) + U(t)` constant, `L(q) - L(y)` equals `U(y) - U(q)`.
#' This uses the form with the smaller ratio of the terms.
#'
#' @param lower_q,lower_y,upper_q,upper_y Log lower and upper tail
#'   quantities at `q` and `y`.
#'
#' @return Vector of the log differences.
#'
#' @keywords internal
.exptilt_tail_diff <- function(lower_q, lower_y, upper_q, upper_y) {
  use_upper <- lower_y - lower_q > upper_q - upper_y
  # NaN from terms that underflow on both sides is a zero difference, and
  # an infinite upper term has no upper tail to use
  use_upper[is.na(use_upper) | is.infinite(upper_q)] <- FALSE
  out <- .log_diff_exp(lower_q, lower_y)
  if (any(use_upper)) {
    out[use_upper] <- .log_diff_exp(upper_y[use_upper], upper_q[use_upper])
  }
  out
}

#' Small window form of the exponentially tilted log CDF
#'
#' See [pcens_cdf_exptilt].
#'
#' @inheritParams .exptilt_lcdf_direct
#'
#' @return Vector of log CDFs.
#'
#' @keywords internal
.exptilt_lcdf_small_window <- function(object, q, pwindow, rho, lower) {
  pos <- .exptilt_layout(q, pwindow, lower)
  moments <- .pcens_tilt_moments(object, pos$endpoints)
  at_q <- moments[pos$q, , drop = FALSE]
  at_y <- moments[pos$y, , drop = FALSE]
  # Scaled by G_1(q) to avoid underflow for small CDFs
  scale <- at_q[, 1]
  g1_y <- exp(at_y[, 1] - scale)
  g2_q <- exp(at_q[, 2] - scale)
  g2_y <- exp(at_y[, 2] - scale)
  g3_q <- exp(at_q[, 3] - scale)
  g3_y <- exp(at_y[, 3] - scale)
  relative <- (1 - g1_y) +
    0.5 * rho * (g2_q - pwindow - g2_y - pwindow * g1_y) +
    rho^2 * ((g3_q - g3_y) / 6 - pwindow * (g2_q + g2_y) / 4 +
      pwindow^2 * (1 - g1_y) / 12)
  out <- scale + log(pmax(relative, 0)) - log(pwindow)
  out[!is.finite(scale)] <- -Inf
  out
}

#' Small delay form of the exponentially tilted log CDF
#'
#' See [pcens_cdf_exptilt].
#'
#' @inheritParams .exptilt_lcdf_direct
#'
#' @return Vector of log CDFs.
#'
#' @keywords internal
.exptilt_lcdf_small_delay <- function(object, q, pwindow, rho) {
  moments <- .pcens_tilt_moments(object, q)
  moments[, 1] + log(abs(rho)) - .exptilt_log_den(rho, pwindow) +
    log1p(rho * exp(moments[, 2] - moments[, 1]) / 2 +
      rho^2 * exp(moments[, 3] - moments[, 1]) / 6)
}
