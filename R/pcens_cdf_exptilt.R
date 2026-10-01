#' Methods for delays with an exponentially tilted primary
#'
#' Analytical primary event censored CDFs for exponential, gamma and normal
#' delay distributions with an exponentially tilted primary event window, the
#' [dexpgrowth()] primary distribution with `r` equal to the tilt \eqn{\rho}.
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
#' The numerical method is less accurate in the lower tail of a gamma with
#' shape below 1, with a relative error of about 1e-5 in R and 1e-4 in Stan.
#'
#' The direct form cancels as \eqn{\rho \to 0}.
#' With \eqn{G_k(t) = \int (t - u)^k f(u) du} up to \eqn{t}, two forms
#' replace it to second order in \eqn{\rho}.
#' * For \eqn{|\rho| w} below \eqn{10^{-2}}, or \eqn{10^{-5}} for delays on
#'   the reals, the uniform window limit with its corrections,
#'   \eqn{\{G_1(q) - G_1(q - w)\} / w +
#'   \rho \{G_2(q) - w G_1(q) - G_2(q - w) - w G_1(q - w)\} / (2 w) +
#'   \rho^2 \{G_3(q) - G_3(q - w) - 3 w [G_2(q) + G_2(q - w)] / 2 +
#'   w^2 [G_1(q) - G_1(q - w)] / 2\} / (6 w)}.
#' * For delays on the non-negative reals with \eqn{q < w} and
#'   \eqn{|\rho| q < 10^{-2}},
#'   \eqn{\rho \{G_1(q) + \rho G_2(q) / 2 + \rho^2 G_3(q) / 6\} /
#'   (e^{\rho w} - 1)}.
#'
#' The truncation error of both forms is about \eqn{(|\rho| w)^4 / 10}, below
#' 3e-8 at their limits.
#' The direct form loses about 1e-14 / \eqn{|\rho| w} for delays on the
#' reals, which is why their limit is lower.
#' For a gamma delay it also loses accuracy in proportion to the shape.
#' The CDF stays within a relative 1e-6 of numerical integration for shapes
#' up to 5000 and \eqn{|\rho|} from 1e-5 and loses accuracy above this.
#'
#' The CDF agrees with numerical integration to a relative difference of
#' about 1e-9 or better in the bulk and the tails, except in the deep lower
#' tail of a normal delay with a small tilt (about 1e-7), for windows much
#' smaller than the delay, which lose up to about 1e-13 / `pwindow` to the
#' difference of the terms at the two endpoints, and for the cases above.
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

# Largest |rho| w for the small tilt forms, see ?pcens_cdf_exptilt
.exptilt_small <- function(lower) {
  if (is.finite(lower)) 1e-2 else 1e-5
}

#' Primary event censored CDF for an exponentially tilted primary
#'
#' Shared implementation of the [pcens_cdf_exptilt] methods.
#'
#' @inheritParams pcens_cdf
#'
#' @return Vector of CDFs.
#'
#' @noRd
.pcens_cdf_exptilt <- function(object, q, pwindow, use_numeric = FALSE) {
  if (isTRUE(use_numeric)) {
    return(pcens_cdf.default(object, q, pwindow, use_numeric))
  }
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
  if (length(pwindow) != 1L || !is.finite(pwindow) || pwindow <= 0 ||
    !.pcens_tilt_available(object, -rho)) {
    return(pcens_cdf.default(object, q, pwindow, use_numeric))
  }

  result <- rep(NA_real_, length(q))
  result[!is.na(q) & q == Inf] <- 1
  result[!is.na(q) & q == -Inf] <- 0
  finite <- which(is.finite(q))
  if (length(finite) > 0L) {
    result[finite] <- .exptilt_cdf_finite(object, q[finite], pwindow, rho)
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
  lower <- .pcens_tilt_lower(object)
  positive <- is.finite(lower)
  log_cdf <- rep(-Inf, length(q))

  active <- !positive | q > lower
  small_limit <- .exptilt_small(lower)
  if (abs(rho) * pwindow < small_limit) {
    small_window <- active
    tiny_delay <- rep(FALSE, length(q))
  } else {
    small_window <- rep(FALSE, length(q))
    tiny_delay <- positive & active & q < pwindow &
      abs(rho) * q < small_limit
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

#' Endpoints and the positions of `q` and `q - pwindow` among them
#'
#' A single `q` skips the de-duplication.
#'
#' @inheritParams .exptilt_endpoints
#'
#' @return List with `endpoints` and the positions `q` and `y`.
#'
#' @noRd
.exptilt_layout <- function(q, pwindow, lower) {
  if (length(q) == 1L) {
    return(list(endpoints = c(q, q - pwindow), q = 1L, y = 2L))
  }
  endpoints <- .exptilt_endpoints(q, pwindow, lower)
  list(
    endpoints = endpoints,
    q = .exptilt_index(q, endpoints, lower),
    y = .exptilt_index(q - pwindow, endpoints, lower)
  )
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
  pos <- .exptilt_layout(q, pwindow, lower)
  endpoints <- pos$endpoints
  endpoint_terms <- cbind(
    .pcens_tilt_transform(object, endpoints, 0),
    .pcens_tilt_transform(object, endpoints, 0, upper = TRUE),
    .pcens_tilt_transform(object, endpoints, -rho),
    .pcens_tilt_transform(object, endpoints, -rho, upper = TRUE)
  )
  at_q <- endpoint_terms[pos$q, , drop = FALSE]
  at_y <- endpoint_terms[pos$y, , drop = FALSE]
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
  # NaN from terms that underflow on both sides is a zero difference
  use_upper[is.na(use_upper)] <- FALSE
  out <- .log_diff_exp(lower_q, lower_y)
  if (any(use_upper)) {
    out[use_upper] <- .log_diff_exp(upper_y[use_upper], upper_q[use_upper])
  }
  out
}

#' Small window form of the exponentially tilted log CDF
#'
#' Used for \eqn{|\rho| w < 10^{-4}}, see [pcens_cdf_exptilt].
#'
#' @inheritParams .exptilt_lcdf_direct
#'
#' @noRd
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
    log1p(rho * exp(moments[, 2] - moments[, 1]) / 2 +
      rho^2 * exp(moments[, 3] - moments[, 1]) / 6)
}
