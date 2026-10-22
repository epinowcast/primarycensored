#' Tilt transforms of delay distributions
#'
#' Internal generics for the delay part of the analytical solutions with an
#' exponentially tilted primary event window, see [pcens_cdf_exptilt].
#' They are built from the tilt transform
#' \eqn{T_f(\xi; \tau) = \int e^{\xi u} f(u) du} up to \eqn{\tau}, from 0
#' for delays on the non-negative reals.
#' A delay is added with a method for each generic, dispatched on the delay
#' class of a `pcens` object.
#' `.pcens_tilt_fits()` and `.pcens_tilt_pair()` have defaults.
#' The Stan equivalents are `check_for_tilt_transform()`,
#' `log_tilt_transform_pair()` and `primarycensored_tilt_moments()`.
#'
#' @param object A `pcens` object.
#'
#' @param t Numeric vector of points.
#'
#' @param xi The tilt, a single number or, for `.pcens_tilt_transform()`, a
#'   vector of the length of `t`. A window with tilt \eqn{\rho} needs
#'   `xi = -rho`.
#'
#' @param upper Logical. If `TRUE`, the transform over `(t, Inf)` rather than
#'   up to `t`.
#'
#' @param pwindow The primary event window.
#'
#' @return
#' * `.pcens_tilt_transform()`: the log of the transform at each `t`.
#' * `.pcens_tilt_available()`: whether the transform is closed form and the
#'   tilted delay exists for `xi`. If not, the numerical method is used.
#' * `.pcens_tilt_fits()`: `TRUE` at each `t` where the transform is used
#'   for a window `pwindow`, otherwise the numerical method. Only called where
#'   `.pcens_tilt_available()` is `TRUE`. The default is `TRUE`.
#' * `.pcens_tilt_pair()`: a matrix with columns `lower` and `upper`, the log
#'   transform up to and beyond each `t`. The default calls
#'   `.pcens_tilt_transform()` twice.
#' * `.pcens_tilt_lower()`: the lower end of the support, 0 or `-Inf`.
#' * `.pcens_tilt_moments()`: a matrix with columns `G1`, `G2` and `G3`, the
#'   log of \eqn{G_k(t) = \int (t - u)^k f(u) du} up to `t`, `-Inf` below
#'   the support.
#'
#' @name pcens_tilt
#' @keywords internal
NULL

#' @rdname pcens_tilt
.pcens_tilt_transform <- function(object, t, xi, upper = FALSE) {
  UseMethod(".pcens_tilt_transform")
}

#' @rdname pcens_tilt
.pcens_tilt_available <- function(object, xi) {
  UseMethod(".pcens_tilt_available")
}

#' @rdname pcens_tilt
.pcens_tilt_lower <- function(object) {
  UseMethod(".pcens_tilt_lower")
}

#' @rdname pcens_tilt
.pcens_tilt_moments <- function(object, t) {
  UseMethod(".pcens_tilt_moments")
}

#' @rdname pcens_tilt
.pcens_tilt_fits <- function(object, xi, t, pwindow = 0) {
  UseMethod(".pcens_tilt_fits")
}

#' @rdname pcens_tilt
.pcens_tilt_pair <- function(object, t, xi) {
  UseMethod(".pcens_tilt_pair")
}

#' @exportS3Method
.pcens_tilt_fits.default <- function(object, xi, t, pwindow = 0) {
  rep(TRUE, length(t))
}

#' @exportS3Method
.pcens_tilt_pair.default <- function(object, t, xi) {
  cbind(
    lower = .pcens_tilt_transform(object, t, xi),
    upper = .pcens_tilt_transform(object, t, xi, upper = TRUE)
  )
}

#' @exportS3Method
.pcens_tilt_available.default <- function(object, xi) {
  FALSE
}

#' Gamma delay parameters of a pcens object
#'
#' The rate is 1 if neither `rate` nor `scale` is given, as in
#' [stats::pgamma()].
#'
#' @inheritParams pcens_tilt
#'
#' @return A list with `shape` and `rate`.
#'
#' @keywords internal
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
  tilted_rate <- rep_len(p$rate - xi, length(t))
  log_total <- p$shape * (log(p$rate) - log(tilted_rate))
  out <- if (upper) log_total else rep(-Inf, length(t))
  positive <- t > 0
  out[positive] <- log_total[positive] + pgamma(
    t[positive] * tilted_rate[positive],
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
#' @inheritParams pcens_tilt
#'
#' @return The rate.
#'
#' @keywords internal
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
  tilted_rate <- rep_len(rate - xi, length(t))
  log_total <- log(rate) - log(tilted_rate)
  out <- if (upper) log_total else rep(-Inf, length(t))
  positive <- t > 0
  out[positive] <- log_total[positive] + if (upper) {
    -tilted_rate[positive] * t[positive]
  } else {
    .log1m_exp(-tilted_rate[positive] * t[positive])
  }
  out
}

#' Normal delay mean and sd of a pcens object, defaulting as in
#' [stats::pnorm()]
#'
#' @inheritParams pcens_tilt
#'
#' @return A list with `mean` and `sd`.
#'
#' @keywords internal
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
  # Completing the square gives a normal with mean mean + xi sd^2. The
  # standardised point is formed as below to avoid rounding.
  p <- .norm_mean_sd(object)
  xi * p$mean + 0.5 * xi^2 * p$sd^2 + stats::pnorm(
    (t - p$mean) / p$sd - xi * p$sd,
    lower.tail = !upper, log.p = TRUE
  )
}

#' Log moments about a point from log partial moments
#'
#' The log of \eqn{G_k(t) = \int_0^t (t - u)^k f(u) du} for `k = 1, 2, 3`
#' from the partial moments \eqn{m_k(t) = \int_0^t u^k f(u) du}, with
#' \eqn{G_1 = t m_0 - m_1}, \eqn{G_2 = t G_1 - (t m_1 - m_2)} and
#' \eqn{G_3 = t G_2 - (t (t m_1 - m_2) - (t m_2 - m_3))}, each a difference
#' of positive integrals.
#'
#' @param log_t Log of the points.
#'
#' @param log_m0,log_m1,log_m2,log_m3 Log partial moments at each point.
#'
#' @return As for `.pcens_tilt_moments()`.
#'
#' @keywords internal
.moments_from_partial <- function(log_t, log_m0, log_m1, log_m2, log_m3) {
  log_g1 <- .log_diff_exp(log_t + log_m0, log_m1)
  log_h <- .log_diff_exp(log_t + log_m1, log_m2)
  log_g2 <- .log_diff_exp(log_t + log_g1, log_h)
  log_a <- .log_diff_exp(log_t + log_m2, log_m3)
  log_b <- .log_diff_exp(log_t + log_h, log_a)
  log_g3 <- .log_diff_exp(log_t + log_g2, log_b)
  cbind(G1 = log_g1, G2 = log_g2, G3 = log_g3)
}

#' Log moments of a gamma delay about a point
#'
#' The log of \eqn{G_k(t) = \int_0^t (t - u)^k f(u) du} for `k = 1, 2, 3`,
#' from gamma CDFs with the shape raised by `k`.
#'
#' @inheritParams pcens_tilt
#'
#' @param shape,rate Gamma delay parameters.
#'
#' @return As for `.pcens_tilt_moments()`.
#'
#' @keywords internal
.gamma_moments <- function(t, shape, rate) {
  tp <- pmax(t, 0)
  log_t <- log(tp)
  log_m0 <- pgamma(tp, shape, rate, log.p = TRUE)
  log_m1 <- log(shape) - log(rate) +
    pgamma(tp, shape + 1, rate, log.p = TRUE)
  log_m2 <- log(shape) + log(shape + 1) - 2 * log(rate) +
    pgamma(tp, shape + 2, rate, log.p = TRUE)
  log_m3 <- log(shape) + log(shape + 1) + log(shape + 2) - 3 * log(rate) +
    pgamma(tp, shape + 3, rate, log.p = TRUE)
  .moments_from_partial(log_t, log_m0, log_m1, log_m2, log_m3)
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
