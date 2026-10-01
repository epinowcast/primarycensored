# Truncated exponential-moment transforms of delay distributions
#
# Internal generics for the delay part of the analytical solutions with an
# exponentially tilted primary event window, built from
# T_f(xi; tau) = int_{-Inf}^{tau} exp(xi u) f(u) du, with lower limit 0 for
# delays on the non-negative reals.
# A delay is added with methods for `.pcens_tilt_lower()`,
# `.pcens_tilt_available()`, `.pcens_tilt_transform()` and
# `.pcens_tilt_moments()`, which dispatch on the delay class of a `pcens`
# object. The Stan equivalents are `check_for_tilt_transform()`,
# `log_tilt_transform_pair()` and `primarycensored_tilt_moments()`.
#   * `.pcens_tilt_transform()`: log T_f at each `t` (up to `t`, or over
#     `(t, Inf)` with `upper = TRUE`) for a tilt `xi`. A window with tilt rho
#     needs xi = -rho.
#   * `.pcens_tilt_available()`: whether the transform is closed form and the
#     tilted delay exists for `xi`, otherwise the numerical method is used.
#   * `.pcens_tilt_lower()`: the lower end of the support, 0 or `-Inf`.
#   * `.pcens_tilt_moments()`: a matrix of the log of the first, second and
#     third moments of the delay about `t`, see ?pcens_cdf_exptilt.
# The methods `.pcens_tilt_ill_conditioned()` and `.pcens_tilt_numeric()` are
# optional and use the numerical method where the closed form loses precision
# to rounding.
#   * `.pcens_tilt_ill_conditioned()`: `TRUE` at each `q` where the closed
#     form, the small window form if `small_window` and otherwise the direct
#     form, is too inaccurate for tilt `rho` and log CDF `log_cdf`. `FALSE` by
#     default.
#   * `.pcens_tilt_numeric()`: the numerical CDF at each `q`, by default
#     `pcens_cdf.default()`.

.pcens_tilt_transform <- function(object, t, xi, upper = FALSE) {
  UseMethod(".pcens_tilt_transform")
}

.pcens_tilt_available <- function(object, xi) {
  UseMethod(".pcens_tilt_available")
}

.pcens_tilt_lower <- function(object) {
  UseMethod(".pcens_tilt_lower")
}

.pcens_tilt_moments <- function(object, t) {
  UseMethod(".pcens_tilt_moments")
}

.pcens_tilt_ill_conditioned <- function(
  object, q, pwindow, rho, log_cdf, small_window
) {
  UseMethod(".pcens_tilt_ill_conditioned")
}

.pcens_tilt_numeric <- function(object, q, pwindow) {
  UseMethod(".pcens_tilt_numeric")
}

#' @exportS3Method
.pcens_tilt_available.default <- function(object, xi) {
  FALSE
}

#' @exportS3Method
.pcens_tilt_ill_conditioned.default <- function(
  object, q, pwindow, rho, log_cdf, small_window
) {
  rep(FALSE, length(q))
}

#' @exportS3Method
.pcens_tilt_numeric.default <- function(object, q, pwindow) {
  pcens_cdf.default(object, q, pwindow)
}

#' Gamma delay parameters of a pcens object
#'
#' The rate is 1 if neither `rate` nor `scale` is given, as in
#' [stats::pgamma()].
#'
#' @param object A `pcens` object.
#'
#' @return A list with `shape` and `rate`.
#'
#' @noRd
.gamma_shape_rate <- function(object) {
  dist_args <- object$args
  shape <- dist_args$shape
  rate <- dist_args$rate
  if (is.null(shape)) {
    stop("shape parameter is required for Gamma distribution", call. = FALSE)
  }
  if (is.null(rate)) {
    rate <- if (is.null(dist_args$scale)) 1 else 1 / dist_args$scale
  }
  list(shape = shape, rate = rate)
}

#' Log-scale arithmetic helpers
#'
#' Vectorised. `.log_diff_exp()` gives `-Inf` rather than `NaN` for a zero or
#' negative difference.
#'
#' @param a,b,x Numeric vectors on the log scale.
#'
#' @return `log(1 - exp(x))`, `log(exp(a) - exp(b))` and
#'   `log(exp(a) + exp(b))` for `.log1m_exp()`, `.log_diff_exp()` and
#'   `.log_sum_exp()`.
#'
#' @noRd
.log1m_exp <- function(x) {
  near_zero <- x > -log(2)
  if (!any(near_zero)) {
    return(log1p(-exp(x)))
  }
  out <- log1p(-exp(x))
  out[near_zero] <- log(-expm1(x[near_zero]))
  out
}

.log_diff_exp <- function(a, b) {
  if (length(a) == length(b) && !anyNA(a) && !anyNA(b) && all(a > b)) {
    return(a + .log1m_exp(b - a))
  }
  n <- max(length(a), length(b))
  a <- rep_len(a, n)
  b <- rep_len(b, n)
  out <- rep(-Inf, n)
  ok <- !is.na(a) & !is.na(b) & a > b
  out[ok] <- a[ok] + .log1m_exp(b[ok] - a[ok])
  out
}

.log_sum_exp <- function(a, b) {
  larger <- pmax(a, b)
  gap <- -abs(a - b)
  # Both -Inf gives NaN
  gap[is.na(gap)] <- 0
  out <- larger + log1p(exp(gap))
  out[larger == -Inf] <- -Inf
  out
}

#' @exportS3Method
.pcens_tilt_lower.pcens_pgamma <- function(object) {
  0
}

#' @exportS3Method
.pcens_tilt_available.pcens_pgamma <- function(object, xi) {
  .gamma_shape_rate(object)$rate - xi > 0
}

#' @exportS3Method
.pcens_tilt_transform.pcens_pgamma <- function(object, t, xi, upper = FALSE) {
  # Gamma CDF with rate - xi times the total (rate / (rate - xi))^shape
  p <- .gamma_shape_rate(object)
  tilted_rate <- p$rate - xi
  log_total <- p$shape * (log(p$rate) - log(tilted_rate))
  out <- rep(if (upper) log_total else -Inf, length(t))
  positive <- t > 0
  out[positive] <- log_total + pgamma(
    t[positive] * tilted_rate,
    shape = p$shape, lower.tail = !upper, log.p = TRUE
  )
  out
}

#' @exportS3Method
.pcens_tilt_lower.pcens_pexp <- function(object) {
  0
}

#' Exponential delay rate of a pcens object, 1 if not given
#'
#' @param object A `pcens` object.
#'
#' @noRd
.exp_rate <- function(object) {
  rate <- object$args$rate
  if (is.null(rate)) 1 else rate
}

#' @exportS3Method
.pcens_tilt_available.pcens_pexp <- function(object, xi) {
  .exp_rate(object) - xi > 0
}

#' @exportS3Method
.pcens_tilt_transform.pcens_pexp <- function(object, t, xi, upper = FALSE) {
  rate <- .exp_rate(object)
  tilted_rate <- rate - xi
  log_total <- log(rate) - log(tilted_rate)
  out <- rep(if (upper) log_total else -Inf, length(t))
  positive <- t > 0
  out[positive] <- log_total + if (upper) {
    -tilted_rate * t[positive]
  } else {
    .log1m_exp(-tilted_rate * t[positive])
  }
  out
}

#' Normal delay mean and sd of a pcens object, defaulting as in
#' [stats::pnorm()]
#'
#' @param object A `pcens` object.
#'
#' @noRd
.norm_mean_sd <- function(object) {
  mu <- object$args$mean
  sigma <- object$args$sd
  list(
    mean = if (is.null(mu)) 0 else mu,
    sd = if (is.null(sigma)) 1 else sigma
  )
}

#' @exportS3Method
.pcens_tilt_lower.pcens_pnorm <- function(object) {
  -Inf
}

#' @exportS3Method
.pcens_tilt_available.pcens_pnorm <- function(object, xi) {
  TRUE
}

#' @exportS3Method
.pcens_tilt_transform.pcens_pnorm <- function(object, t, xi, upper = FALSE) {
  # Completing the square gives a normal with mean mean + xi sd^2
  p <- .norm_mean_sd(object)
  xi * p$mean + 0.5 * xi^2 * p$sd^2 + stats::pnorm(
    t,
    mean = p$mean + xi * p$sd^2, sd = p$sd,
    lower.tail = !upper, log.p = TRUE
  )
}

#' Log moments of a gamma delay about a point
#'
#' The log of \eqn{G_k(t) = \int_0^t (t - u)^k f(u) du} for `k = 1, 2, 3`,
#' from gamma CDFs with the shape raised by `k`.
#'
#' @param t Numeric vector of finite points.
#'
#' @param shape,rate Gamma delay parameters.
#'
#' @return A matrix with columns `G1`, `G2` and `G3`, `-Inf` for `t <= 0`.
#'
#' @noRd
.gamma_moments <- function(t, shape, rate) {
  positive <- t > 0
  tp <- pmax(t, 0)
  log_t <- log(tp)
  log_m0 <- pgamma(tp, shape, rate, log.p = TRUE)
  log_m1 <- log(shape) - log(rate) +
    pgamma(tp, shape + 1, rate, log.p = TRUE)
  log_m2 <- log(shape) + log(shape + 1) - 2 * log(rate) +
    pgamma(tp, shape + 2, rate, log.p = TRUE)
  log_m3 <- log(shape) + log(shape + 1) + log(shape + 2) - 3 * log(rate) +
    pgamma(tp, shape + 3, rate, log.p = TRUE)
  log_g1 <- .log_diff_exp(log_t + log_m0, log_m1)
  log_h <- .log_diff_exp(log_t + log_m1, log_m2)
  log_g2 <- .log_diff_exp(log_t + log_g1, log_h)
  log_a <- .log_diff_exp(log_t + log_m2, log_m3)
  log_b <- .log_diff_exp(log_t + log_h, log_a)
  log_g3 <- .log_diff_exp(log_t + log_g2, log_b)
  cbind(
    G1 = ifelse(positive, log_g1, -Inf),
    G2 = ifelse(positive, log_g2, -Inf),
    G3 = ifelse(positive, log_g3, -Inf)
  )
}

#' @exportS3Method
.pcens_tilt_moments.pcens_pgamma <- function(object, t) {
  p <- .gamma_shape_rate(object)
  .gamma_moments(t, p$shape, p$rate)
}

#' @exportS3Method
.pcens_tilt_moments.pcens_pexp <- function(object, t) {
  # Gamma with shape 1, as the closed forms cancel for small rate * t
  .gamma_moments(t, 1, .exp_rate(object))
}

#' @exportS3Method
.pcens_tilt_moments.pcens_pnorm <- function(object, t) {
  p <- .norm_mean_sd(object)
  z <- (t - p$mean) / p$sd
  log_phi <- stats::dnorm(z, log = TRUE)
  log_Phi <- stats::pnorm(z, log.p = TRUE)
  log_abs_z <- log(abs(z))
  below <- z < 0
  log_g1 <- ifelse(
    below,
    .log_diff_exp(log_phi, log_abs_z + log_Phi),
    .log_sum_exp(log_phi, log_abs_z + log_Phi)
  )
  log_g2 <- ifelse(
    below,
    .log_diff_exp(log(z^2 + 1) + log_Phi, log_abs_z + log_phi),
    .log_sum_exp(log(z^2 + 1) + log_Phi, log_abs_z + log_phi)
  )
  log_g3 <- ifelse(
    below,
    .log_diff_exp(
      log(z^2 + 2) + log_phi, log_abs_z + log(z^2 + 3) + log_Phi
    ),
    .log_sum_exp(
      log(z^2 + 2) + log_phi, log_abs_z + log(z^2 + 3) + log_Phi
    )
  )
  cbind(
    G1 = log(p$sd) + log_g1,
    G2 = 2 * log(p$sd) + log_g2,
    G3 = 3 * log(p$sd) + log_g3
  )
}
