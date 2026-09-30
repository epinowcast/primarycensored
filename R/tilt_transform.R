#' Truncated exponential-moment transforms of delay distributions
#'
#' The primary event censored CDF for several non-uniform primary event
#' windows is a sum of terms built from
#' \deqn{T_f(\xi; \tau) = \int_{-\infty}^{\tau} e^{\xi u} f(u) du,}
#' the truncated exponential-moment transform of the delay density \eqn{f}.
#' The lower limit is 0 for delay distributions on the non-negative reals.
#' The primary event window fixes the tilts \eqn{\xi} and the coefficients.
#' The delay distribution fixes whether \eqn{T_f} is closed form.
#'
#' These internal generics are the extension points for new delay families.
#' Each is dispatched on the delay class of a `pcens` object, for example
#' `pcens_pgamma`, so one set of methods serves every primary event window.
#' A family is added by defining
#' * `.pcens_tilt_lower()`, the lower end of the support.
#' * `.pcens_tilt_available()`, whether the closed form applies for a tilt.
#' * `.pcens_tilt_transform()`, the transform on the log scale.
#' * Optionally `.pcens_tilt_moments()`, which is only needed by the small
#'   tilt forms of [pcens_cdf.pcens_pexp_dexpgrowth()].
#'
#' The Stan equivalents are `check_for_tilt_transform()`,
#' `log_tilt_transform()`, `log_tilt_transform_upper()` and
#' `primarycensored_tilt_moments()`.
#'
#' @param object A `pcens` object as created by [new_pcens()].
#'
#' @param t Numeric vector of finite points at which to evaluate the
#'   transform.
#'
#' @param xi Tilt \eqn{\xi}, a single number. The exponentially tilted window
#'   with tilt \eqn{\rho} needs \eqn{\xi = -\rho}, and \eqn{\xi = 0} gives the
#'   delay CDF.
#'
#' @param upper Logical. If `TRUE` return the transform over \eqn{(t, \infty)}
#'   rather than over the lower end of the support up to `t`. Evaluating the
#'   tail directly keeps precision where the lower transform is close to its
#'   total.
#'
#' @return
#' * `.pcens_tilt_transform()`: the log of the transform at each `t`. It is
#'   `-Inf` below the support for the lower transform.
#' * `.pcens_tilt_available()`: `TRUE` if the transform is closed form and the
#'   tilted delay distribution exists for `xi`, otherwise `FALSE`. Callers use
#'   the numerical method when it is `FALSE`.
#' * `.pcens_tilt_lower()`: the lower end of the support, 0 or `-Inf`.
#' * `.pcens_tilt_moments()`: a matrix with two columns, the log of the first
#'   and second moments of the delay about `t`, see
#'   [pcens_cdf.pcens_pexp_dexpgrowth()].
#'
#' @family tilt
#'
#' @keywords internal
#' @name tilt_transform
NULL

#' @rdname tilt_transform
.pcens_tilt_transform <- function(object, t, xi, upper = FALSE) {
  UseMethod(".pcens_tilt_transform")
}

#' @rdname tilt_transform
.pcens_tilt_available <- function(object, xi) {
  UseMethod(".pcens_tilt_available")
}

#' @rdname tilt_transform
.pcens_tilt_lower <- function(object) {
  UseMethod(".pcens_tilt_lower")
}

#' @rdname tilt_transform
.pcens_tilt_moments <- function(object, t) {
  UseMethod(".pcens_tilt_moments")
}

#' @rdname tilt_transform
.pcens_tilt_available.default <- function(object, xi) {
  FALSE
}

#' @rdname tilt_transform
.pcens_tilt_transform.default <- function(object, t, xi, upper = FALSE) {
  stop(
    "No tilt transform is available for this delay distribution.",
    call. = FALSE
  )
}

#' @rdname tilt_transform
.pcens_tilt_lower.default <- function(object) {
  stop(
    "No tilt transform is available for this delay distribution.",
    call. = FALSE
  )
}

#' @rdname tilt_transform
.pcens_tilt_moments.default <- function(object, t) {
  stop(
    "No tilt transform is available for this delay distribution.",
    call. = FALSE
  )
}

#' Gamma delay parameters of a pcens object
#'
#' Takes the `shape` and the `rate` or `scale` from the delay arguments as
#' the gamma analytical solutions for the uniform primary do.
#'
#' @inheritParams tilt_transform
#'
#' @return A list with `shape` and `rate`.
#'
#' @keywords internal
.gamma_shape_rate <- function(object) {
  shape <- object$args$shape
  scale <- object$args$scale
  rate <- object$args$rate
  if (is.null(shape)) {
    stop("shape parameter is required for Gamma distribution", call. = FALSE)
  }
  if (is.null(rate)) {
    if (is.null(scale)) {
      stop(
        "scale or rate parameter is required for Gamma distribution",
        call. = FALSE
      )
    }
    rate <- 1 / scale
  }
  list(shape = shape, rate = rate)
}

# Log-scale helpers. All are vectorised and return -Inf, rather than NaN,
# when the difference in `.log_diff_exp()` is zero or rounds to a negative.

#' Log-scale arithmetic helpers
#'
#' @param a,b Numeric vectors on the log scale.
#'
#' @param x Numeric vector of values at most 0 on the log scale.
#'
#' @return
#' * `.log1m_exp()`: `log(1 - exp(x))`.
#' * `.log_diff_exp()`: `log(exp(a) - exp(b))`, or `-Inf` where `a <= b`.
#' * `.log_sum_exp()`: `log(exp(a) + exp(b))`.
#'
#' @keywords internal
#' @name log_helpers
.log1m_exp <- function(x) {
  ifelse(x > -log(2), log(-expm1(x)), log1p(-exp(x)))
}

#' @rdname log_helpers
.log_diff_exp <- function(a, b) {
  n <- max(length(a), length(b))
  a <- rep_len(a, n)
  b <- rep_len(b, n)
  out <- rep(-Inf, n)
  ok <- !is.na(a) & !is.na(b) & a > b
  out[ok] <- a[ok] + .log1m_exp(b[ok] - a[ok])
  out
}

#' @rdname log_helpers
.log_sum_exp <- function(a, b) {
  larger <- pmax(a, b)
  gap <- -abs(a - b)
  # Both -Inf gives NaN
  gap[is.na(gap)] <- 0
  ifelse(larger == -Inf, -Inf, larger + log1p(exp(gap)))
}

#' @rdname tilt_transform
.pcens_tilt_lower.pcens_pgamma <- function(object) {
  0
}

#' @rdname tilt_transform
.pcens_tilt_available.pcens_pgamma <- function(object, xi) {
  .gamma_shape_rate(object)$rate - xi > 0
}

#' @rdname tilt_transform
.pcens_tilt_transform.pcens_pgamma <- function(object, t, xi, upper = FALSE) {
  # The tilted density is proportional to a gamma density with rate
  # rate - xi, so the transform is the tilted gamma CDF times the total
  # T_f(xi; Inf) = (rate / (rate - xi))^shape.
  p <- .gamma_shape_rate(object)
  tilted_rate <- p$rate - xi
  log_total <- p$shape * (log(p$rate) - log(tilted_rate))
  out <- rep(if (upper) log_total else -Inf, length(t))
  positive <- t > 0
  out[positive] <- log_total + stats::pgamma(
    t[positive] * tilted_rate,
    shape = p$shape, lower.tail = !upper, log.p = TRUE
  )
  out
}

#' @rdname tilt_transform
.pcens_tilt_lower.pcens_pexp <- function(object) {
  0
}

#' Exponential delay rate of a pcens object
#'
#' @inheritParams tilt_transform
#'
#' @return The rate, 1 if not given as in [stats::pexp()].
#'
#' @keywords internal
.exp_rate <- function(object) {
  rate <- object$args$rate
  if (is.null(rate)) 1 else rate
}

#' @rdname tilt_transform
.pcens_tilt_available.pcens_pexp <- function(object, xi) {
  .exp_rate(object) - xi > 0
}

#' @rdname tilt_transform
.pcens_tilt_transform.pcens_pexp <- function(object, t, xi, upper = FALSE) {
  # T_f(xi; t) = rate / (rate - xi) * (1 - exp(-(rate - xi) t))
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

#' Normal delay parameters of a pcens object
#'
#' @inheritParams tilt_transform
#'
#' @return A list with `mean` and `sd`, defaulting as in [stats::pnorm()].
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

#' @rdname tilt_transform
.pcens_tilt_lower.pcens_pnorm <- function(object) {
  -Inf
}

#' @rdname tilt_transform
.pcens_tilt_available.pcens_pnorm <- function(object, xi) {
  TRUE
}

#' @rdname tilt_transform
.pcens_tilt_transform.pcens_pnorm <- function(object, t, xi, upper = FALSE) {
  # Completing the square gives a normal density with mean
  # mean + xi sd^2, so the transform is
  # exp(xi mean + xi^2 sd^2 / 2) Phi((t - mean - xi sd^2) / sd).
  p <- .norm_mean_sd(object)
  xi * p$mean + 0.5 * xi^2 * p$sd^2 + stats::pnorm(
    t,
    mean = p$mean + xi * p$sd^2, sd = p$sd,
    lower.tail = !upper, log.p = TRUE
  )
}

#' Log moments of a gamma delay about a point
#'
#' For `t > 0` these are the log of
#' \eqn{G_1(t) = \int_0^t (t - u) f(u) du} and
#' \eqn{G_2(t) = \int_0^t (t - u)^2 f(u) du}.
#' They come from the CDFs of gamma distributions with the shape raised by one
#' and two, which give the partial moments of the delay.
#'
#' @param t Numeric vector of finite points.
#'
#' @param shape,rate Gamma delay parameters.
#'
#' @return A matrix with columns `G1` and `G2` on the log scale, `-Inf` for
#'   `t <= 0`.
#'
#' @keywords internal
.gamma_moments <- function(t, shape, rate) {
  positive <- t > 0
  tp <- pmax(t, 0)
  log_t <- log(tp)
  log_m0 <- stats::pgamma(tp, shape, rate, log.p = TRUE)
  # Partial first and second moments of the delay
  log_m1 <- log(shape) - log(rate) +
    stats::pgamma(tp, shape + 1, rate, log.p = TRUE)
  log_m2 <- log(shape) + log(shape + 1) - 2 * log(rate) +
    stats::pgamma(tp, shape + 2, rate, log.p = TRUE)
  # All differences are of positive integrals, for example
  # G_1 = t F - m_1 and G_2 = t G_1 - (t m_1 - m_2)
  log_g1 <- .log_diff_exp(log_t + log_m0, log_m1)
  log_h <- .log_diff_exp(log_t + log_m1, log_m2)
  log_g2 <- .log_diff_exp(log_t + log_g1, log_h)
  cbind(
    G1 = ifelse(positive, log_g1, -Inf),
    G2 = ifelse(positive, log_g2, -Inf)
  )
}

#' @rdname tilt_transform
.pcens_tilt_moments.pcens_pgamma <- function(object, t) {
  p <- .gamma_shape_rate(object)
  .gamma_moments(t, p$shape, p$rate)
}

#' @rdname tilt_transform
.pcens_tilt_moments.pcens_pexp <- function(object, t) {
  # The exponential is the gamma distribution with shape 1. Its own closed
  # forms cancel when the rate times t is small.
  .gamma_moments(t, 1, .exp_rate(object))
}

#' @rdname tilt_transform
.pcens_tilt_moments.pcens_pnorm <- function(object, t) {
  # With z = (t - mean) / sd,
  # G_1 = sd (phi(z) + z Phi(z)) and G_2 = sd^2 ((z^2 + 1) Phi(z) + z phi(z))
  p <- .norm_mean_sd(object)
  z <- (t - p$mean) / p$sd
  log_phi <- stats::dnorm(z, log = TRUE)
  log_Phi <- stats::pnorm(z, log.p = TRUE)
  log_abs_z <- log(abs(z))
  below <- z < 0
  # For z < 0 both terms of each sum have opposite signs and the first is
  # the larger
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
  cbind(G1 = log(p$sd) + log_g1, G2 = 2 * log(p$sd) + log_g2)
}
