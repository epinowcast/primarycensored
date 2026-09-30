#' Methods for delays with an exponentially tilted primary
#'
#' Analytical primary event censored CDFs for exponential, gamma, normal and
#' lognormal delay distributions with an exponentially tilted primary event
#' window, the [dexpgrowth()] primary distribution with `r` equal to the tilt
#' \eqn{\rho}.
#' They honour `use_numeric`, and use the numerical method of
#' [pcens_cdf.default()] when no closed form applies.
#'
#' @inheritParams pcens_cdf
#'
#' @details
#' With window width \eqn{w}, delay CDF \eqn{F} and the transform
#' \eqn{J(x) = T_f(-\rho; x)} of [tilt_transform], the CDF at \eqn{q} is
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
#' In Stan the numerical method is less accurate in the lower tail of a
#' gamma with shape below 1.
#' The lognormal transform has no closed form.
#' It is evaluated by quadrature for \eqn{\rho > 0} and by a series for
#' \eqn{\rho < 0}, see [tilt_transform_lognormal].
#' The numerical method is used where the transform does not apply or is
#' slower, see `.pcens_tilt_fits()`.
#'
#' The direct form cancels as \eqn{\rho \to 0}.
#' With \eqn{G_k(t) = \int (t - u)^k f(u) du} up to \eqn{t}, two forms
#' replace it.
#' * For \eqn{|\rho| w < 10^{-4}} and \eqn{|\rho| (|q| + w) < 0.1}, the
#'   uniform window limit with its first order correction in \eqn{\rho},
#'   \eqn{\{G_1(q) - G_1(q - w)\} / w +
#'   \rho \{G_2(q) - w G_1(q) - G_2(q - w) - w G_1(q - w)\} / (2 w)}.
#' * For delays on the non-negative reals with \eqn{q < w} and
#'   \eqn{|\rho| q < 10^{-4}},
#'   \eqn{\rho \{G_1(q) + \rho G_2(q) / 2\} / (e^{\rho w} - 1)}.
#'
#' The value of both forms has a truncation error below 1e-9.
#' The Stan gradient in \eqn{\rho} of both has a relative error of up to
#' about 2e-5 at the thresholds.
#'
#' The CDF agrees with numerical integration to a relative difference of
#' about 1e-9 or better, except in the deep lower tail of a normal delay with
#' a small tilt (about 1e-7) and for windows much smaller than the delay,
#' which lose up to about 1e-13 / `pwindow` to the difference of the terms at
#' the two endpoints.
#'
#' @inherit pcens_cdf return
#'
#' @name pcens_cdf_exptilt
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

# The direct form loses about 1e-14 / (|rho| w) and the small tilt forms
# have a truncation error of about (|rho| w)^2 / 12, which balance near 1e-4
.exptilt_small <- 1e-4

# The small window form cancels beyond this |rho| (|q| + w)
.exptilt_small_reach <- 0.1

#' Tilt of the exponentially tilted primary of a pcens object
#'
#' @inheritParams pcens_cdf
#'
#' @return The tilt `r`, a single finite number.
#'
#' @noRd
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
#' @noRd
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
  if (length(q) < min_q && abs(rho) * pwindow <= min_xw) {
    return(pcens_cdf.default(object, q, pwindow, use_numeric))
  }

  result <- rep(NA_real_, length(q))
  result[!is.na(q) & q == Inf] <- 1
  result[!is.na(q) & q == -Inf] <- 0
  finite <- which(is.finite(q))
  if (length(finite) > 0L) {
    # Use the numerical method for the q where the transform does not fit
    fits <- .pcens_tilt_fits(object, -rho, q[finite], pwindow)
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
#' @inheritParams pcens_cdf
#'
#' @param rho The tilt, the `r` of the exponential growth primary.
#'
#' @return Vector of CDFs, clamped to \[0, 1\].
#'
#' @noRd
.exptilt_cdf_finite <- function(object, q, pwindow, rho) {
  if (length(q) == 0L) {
    return(numeric(0))
  }
  lower <- .pcens_tilt_lower(object)
  positive <- is.finite(lower)
  log_cdf <- rep(-Inf, length(q))

  active <- !positive | q > lower
  if (abs(rho) * pwindow < .exptilt_small) {
    small_window <- active &
      abs(rho) * (abs(q) + pwindow) < .exptilt_small_reach
    tiny_delay <- rep(FALSE, length(q))
  } else {
    small_window <- rep(FALSE, length(q))
    tiny_delay <- positive & active & q < pwindow &
      abs(rho) * q < .exptilt_small
  }
  direct <- active & !small_window & !tiny_delay

  if (any(small_window)) {
    log_cdf[small_window] <- .exptilt_lcdf_small_window(
      object, q[small_window], pwindow, rho, lower
    )
  }
  if (any(tiny_delay)) {
    log_cdf[tiny_delay] <- .exptilt_lcdf_tiny_delay(
      object, q[tiny_delay], pwindow, rho
    )
  }
  if (any(direct)) {
    log_cdf[direct] <- .exptilt_lcdf_direct(
      object, q[direct], pwindow, rho, lower
    )
  }
  pmin(1, exp(log_cdf))
}

#' Unique endpoints `q` and `q - pwindow` at which the transforms are needed
#'
#' Endpoints at or below a finite `lower` share one entry.
#'
#' @param q Numeric vector of finite quantiles.
#'
#' @param pwindow Primary event window.
#'
#' @param lower Lower end of the support of the delay, 0 or `-Inf`.
#'
#' @return Sorted numeric vector of unique endpoints.
#'
#' @noRd
.exptilt_endpoints <- function(q, pwindow, lower) {
  endpoints <- sort.int(unique(c(q, q - pwindow)))
  if (is.finite(lower)) {
    endpoints <- c(lower, endpoints[endpoints > lower])
  }
  endpoints
}

#' Position of endpoints in the output of `.exptilt_endpoints()`
#'
#' @param t Numeric vector of endpoints.
#'
#' @param endpoints,lower As for `.exptilt_endpoints()`.
#'
#' @return Integer vector of positions.
#'
#' @noRd
.exptilt_index <- function(t, endpoints, lower) {
  if (is.finite(lower)) {
    t <- pmax(t, lower)
  }
  match(t, endpoints)
}

#' Direct form of the exponentially tilted log CDF
#'
#' The expression of [pcens_cdf_exptilt] on the log scale.
#'
#' @param object A `pcens` object.
#'
#' @param q,pwindow,rho,lower As for `.exptilt_cdf_finite()`.
#'
#' @return Vector of log CDFs.
#'
#' @noRd
.exptilt_lcdf_direct <- function(object, q, pwindow, rho, lower) {
  endpoints <- .exptilt_endpoints(q, pwindow, lower)
  endpoint_terms <- cbind(
    .pcens_tilt_pair(object, endpoints, 0),
    .pcens_tilt_pair(object, endpoints, -rho)
  )
  at_q <- endpoint_terms[
    .exptilt_index(q, endpoints, lower), ,
    drop = FALSE
  ]
  at_y <- endpoint_terms[
    .exptilt_index(q - pwindow, endpoints, lower), ,
    drop = FALSE
  ]
  log_diff_f <- .exptilt_tail_diff(
    at_q[, 1], at_y[, 1], at_q[, 2], at_y[, 2]
  )
  log_diff_j <- .exptilt_tail_diff(
    at_q[, 3], at_y[, 3], at_q[, 4], at_y[, 4]
  )
  if (rho > 0) {
    log_num <- .log_diff_exp(rho * q + log_diff_j, log_diff_f)
    log_den <- .log_diff_exp(rho * pwindow, 0)
  } else {
    log_num <- .log_diff_exp(log_diff_f, rho * q + log_diff_j)
    log_den <- .log1m_exp(rho * pwindow)
  }
  .log_sum_exp(at_y[, 1], log_num - log_den)
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
#' @noRd
.exptilt_tail_diff <- function(lower_q, lower_y, upper_q, upper_y) {
  use_upper <- lower_y - lower_q > upper_q - upper_y
  # Terms that underflow on both sides give NaN, which is a zero difference
  use_upper[is.na(use_upper) | is.infinite(upper_q)] <- FALSE
  out <- .log_diff_exp(lower_q, lower_y)
  if (any(use_upper)) {
    out[use_upper] <- .log_diff_exp(upper_y[use_upper], upper_q[use_upper])
  }
  out
}

#' Small window form of the exponentially tilted log CDF
#'
#' Used for \eqn{|\rho| w < 10^{-4}} and \eqn{|\rho| (|q| + w) < 0.1}, see
#' [pcens_cdf_exptilt].
#'
#' @inheritParams .exptilt_lcdf_direct
#'
#' @noRd
.exptilt_lcdf_small_window <- function(object, q, pwindow, rho, lower) {
  endpoints <- .exptilt_endpoints(q, pwindow, lower)
  moments <- .pcens_tilt_moments(object, endpoints)
  at_q <- moments[.exptilt_index(q, endpoints, lower), , drop = FALSE]
  at_y <- moments[.exptilt_index(q - pwindow, endpoints, lower), ,
    drop = FALSE
  ]
  # Scaled by G_1(q) to avoid underflow for small CDFs
  scale <- at_q[, 1]
  g1_y <- exp(at_y[, 1] - scale)
  g2_q <- exp(at_q[, 2] - scale)
  g2_y <- exp(at_y[, 2] - scale)
  relative <- (1 - g1_y) +
    0.5 * rho * (g2_q - pwindow - g2_y - pwindow * g1_y)
  out <- scale + log(pmax(relative, 0)) - log(pwindow)
  out[!is.finite(scale)] <- -Inf
  out
}

#' Small delay form of the exponentially tilted log CDF
#'
#' Used for delays on the non-negative reals with `q < pwindow` and
#' \eqn{|\rho| q < 10^{-4}}, see [pcens_cdf_exptilt].
#'
#' @inheritParams .exptilt_lcdf_direct
#'
#' @noRd
.exptilt_lcdf_tiny_delay <- function(object, q, pwindow, rho) {
  moments <- .pcens_tilt_moments(object, q)
  if (rho > 0) {
    log_den <- .log_diff_exp(rho * pwindow, 0)
  } else {
    log_den <- .log1m_exp(rho * pwindow)
  }
  moments[, 1] + log(abs(rho)) - log_den +
    log1p(0.5 * rho * exp(moments[, 2] - moments[, 1]))
}
