#' Methods for log-logistic delays
#'
#' Analytical primary event censored CDFs for a log-logistic delay with an
#' exponentially tilted ([dexpgrowth()]) or a uniform primary event window.
#' They honour `use_numeric`. The numerical method is used where the closed
#' form does not apply or is ill conditioned.
#'
#' @inheritParams pcens_cdf
#'
#' @details
#' The delay has \eqn{F(x) = 1 / (1 + (x / \lambda)^{-k})}, given as in
#' `flexsurv::pllogis()` and `actuar::pllogis()` with `shape` and `scale`,
#' or `rate` for the latter. The Stan functions use `[scale, shape]`, the
#' order of Stan's `loglogistic_cdf()`.
#'
#' **Uniform window.** With \eqn{G_1(t) = \int_0^t F(u) du} the CDF is
#' \eqn{\{G_1(q) - G_1(q - w)\} / w}, where
#' \eqn{G_1(t) = t F(t) (1 - r_{1/k}(A))} and \eqn{A = (t / \lambda)^k}.
#' Each unique endpoint is evaluated once.
#'
#' **Tilted window.** The transform of [tilt_transform] is
#' \eqn{T_f(\xi; t) = F(t) \sum_n (\xi t)^n / n!\, r_{n/k}(A)}, with
#' \deqn{r_a(A) = \frac{1 + A}{A^{a + 1}} \int_0^A \frac{y^a}{(1 + y)^2} dy.}
#' The ratio is at most \eqn{1 / (a + 1)}.
#' Base R has no Gauss hypergeometric function, and the incomplete beta form
#' holds only for \eqn{a < 1}.
#' `.loglogistic_ratio()` therefore uses three series that do not cancel,
#' for \eqn{A \le 1}, \eqn{1 < A \le 3} and \eqn{A > 3}.
#' The upper transform is `NaN` for \eqn{\xi \ne 0}, so the combination of
#' terms takes differences of the lower transform.
#'
#' **Fallback to the numerical method.** The series is used where
#' \eqn{|\xi| t \le 10} and the shape is at least 0.2. The uniform
#' solution and the small tilt form are used for a shape of at least 0.01.
#' Elsewhere the numerical method is used.
#' Within these limits a bound on the rounding error of the closed form is
#' compared with the smaller of the CDF and the survival. The numerical
#' method is used where it exceeds 1e-8 of that tail, or 1e-6 of the density
#' of the CDF for the solutions that difference \eqn{G_1}.
#' The bound is not taken below 1e-15 in the upper tail.
#' The same rule is used in Stan.
#' The numerical method integrates the delay CDF, or the survival where the
#' delay CDF at the mean primary event time is above a half, with a relative
#' tolerance of 1e-13.
#' The switch between methods makes the CDF discontinuous by about 1e-8 of
#' the smaller tail.
#'
#' @inherit pcens_cdf return
#'
#' @name pcens_cdf_loglogistic
#'
#' @concept pcens
#'
#' @examples
#' # Log-logistic delay with a growing primary event process
#' pllogis <- add_name_attribute(
#'   function(q, shape, scale) {
#'     plogis(shape * (log(pmax(q, 0)) - log(scale)))
#'   },
#'   "pllogis"
#' )
#' obj <- new_pcens(
#'   pdist = pllogis, dprimary = dexpgrowth,
#'   primary_args = list(r = 0.3), shape = 2, scale = 5
#' )
#' pcens_cdf(obj, q = c(1, 4, 8), pwindow = 2)
#'
#' # The same delay with a uniform primary
#' obj <- new_pcens(
#'   pdist = pllogis, dprimary = dunif,
#'   primary_args = list(), shape = 2, scale = 5
#' )
#' pcens_cdf(obj, q = c(1, 4, 8), pwindow = 2)
NULL

# Smallest shape for which any solution is used, as the ratios involve
# powers of 1 / shape that overflow when tiny
.loglogistic_floor_shape <- 0.01

# The numerical method is used where an error bound exceeds this fraction of
# the smaller of the CDF and the survival, or of the density of the CDF
.loglogistic_precision <- 1e-8
.loglogistic_pmf_precision <- 1e-6

# Rounding error of a difference of G_1 relative to 4e-16 t F(t)
.loglogistic_error_constant <- 4e-16

# A above which the tail series of the partial moments is used
.loglogistic_tail_start <- 3

#' Shape and scale of a log-logistic pcens object
#'
#' The scale is 1 if neither `scale` nor `rate` is given.
#'
#' @noRd
.loglogistic_shape_scale <- function(object) {
  shape <- object$args$shape
  scale <- object$args$scale
  rate <- object$args$rate
  if (is.null(shape)) {
    stop(
      "shape parameter is required for the log-logistic distribution",
      call. = FALSE
    )
  }
  if (is.null(scale)) {
    scale <- if (is.null(rate)) 1 else 1 / rate
  }
  list(shape = shape, scale = scale)
}

#' Log of one plus an exponential without overflow
#'
#' @noRd
.log1p_exp <- function(x) {
  pmax(x, 0) + log1p(exp(-abs(x)))
}

#' Partial moment ratio for A up to 1, see `.loglogistic_ratio()`
#'
#' @noRd
.loglogistic_ratio_small <- function(a, log_A) {
  A <- exp(log_A)
  nr <- length(A)
  nc <- length(a)
  a_mat <- matrix(a, nr, nc, byrow = TRUE)
  b_mat <- matrix(A / (1 + A), nr, nc)
  weight <- exp(-a_mat * matrix(log1p(A), nr, nc))
  total <- weight / (a_mat + 1)
  m <- 0
  repeat {
    weight <- weight * b_mat * (a_mat + m) / (m + 1)
    m <- m + 1
    term <- weight / (a_mat + m + 1)
    total <- total + term
    if (all(term <= 1e-16 * total) || m >= 10000) {
      break
    }
  }
  total
}

#' Partial moment M_a(1) from the digamma function
#'
#' @noRd
.loglogistic_partial_at_one <- function(a) {
  out <- rep(0.5, length(a))
  positive <- a > 0
  ap <- a[positive]
  out[positive] <- 0.5 *
    (ap * (digamma((ap + 1) / 2) - digamma(ap / 2)) - 1)
  out
}

#' Power series of M_a(A) - M_a(1) in tau = (A - 1) / (A + 1)
#'
#' The recurrence runs on the terms times tau^j so large `a` does not
#' overflow.
#'
#' @noRd
.loglogistic_tau_integral <- function(a, tau) {
  nr <- length(tau)
  nc <- length(a)
  tau_mat <- matrix(tau, nr, nc)
  shift <- matrix(2 * a, nr, nc, byrow = TRUE) * tau_mat
  tau_sq <- tau_mat^2
  previous <- matrix(0, nr, nc)
  current <- matrix(1, nr, nc)
  total <- current
  j <- 0
  repeat {
    following <- (shift * current + (j - 1) * tau_sq * previous) / (j + 1)
    previous <- current
    current <- following
    j <- j + 1
    term <- current / (j + 1)
    total <- total + term
    if (all(term <= 1e-16 * total) || j >= 10000) {
      break
    }
  }
  0.5 * tau_mat * total
}

#' Stable (1 - exp(-z)) / z
#'
#' @noRd
.exprel_neg <- function(z) {
  out <- z
  out[] <- 1
  nonzero <- z != 0
  out[nonzero] <- -expm1(-z[nonzero]) / z[nonzero]
  out
}

#' Partial moment ratio for A above the tail start
#'
#' Adds the integral from 3 to `A` to M_a(3), see `.loglogistic_ratio()`.
#'
#' @noRd
.loglogistic_ratio_large <- function(a, log_A) {
  nr <- length(log_A)
  nc <- length(a)
  log_start <- log(.loglogistic_tail_start)
  a_mat <- matrix(a, nr, nc, byrow = TRUE)
  ell <- matrix(log_A - log_start, nr, nc)
  at_start <- (.loglogistic_partial_at_one(a) +
    .loglogistic_tau_integral(a, 0.5)[1L, ]) * exp(-a * log_start)
  shrink <- exp(-a_mat * ell)
  total <- shrink * matrix(at_start, nr, nc, byrow = TRUE)
  decay <- exp(-log_A)
  decay_m <- decay
  start_m <- 1 / .loglogistic_tail_start
  m <- 0L
  repeat {
    x <- a_mat - m - 1
    term <- (matrix(decay_m, nr, nc) - shrink * start_m) / x
    near_zero <- which(abs(x * ell) < 1)
    if (length(near_zero) > 0L) {
      term[near_zero] <- matrix(decay_m, nr, nc)[near_zero] *
        ell[near_zero] * .exprel_neg((x * ell)[near_zero])
    }
    term <- (if (m %% 2L == 0L) 1 else -1) * (m + 1) * term
    total <- total + term
    m <- m + 1L
    decay_m <- decay_m * decay
    start_m <- start_m / .loglogistic_tail_start
    if (all(abs(term) <= 1e-17 * abs(total)) || m >= 200L) {
      break
    }
  }
  (1 + decay) * total
}

#' Partial moment ratio r_a(A) of the log-logistic
#'
#' Rows are `log_A` and columns are `a`, see [pcens_cdf_loglogistic].
#'
#' @noRd
.loglogistic_ratio <- function(a, log_A) {
  out <- matrix(NA_real_, length(log_A), length(a))
  small <- log_A <= 0
  large <- log_A > log(.loglogistic_tail_start)
  mid <- !small & !large
  if (any(small)) {
    out[small, ] <- .loglogistic_ratio_small(a, log_A[small])
  }
  if (any(mid)) {
    A <- exp(log_A[mid])
    moment <- matrix(
      .loglogistic_partial_at_one(a), length(A), length(a),
      byrow = TRUE
    ) + .loglogistic_tau_integral(a, (A - 1) / (A + 1))
    out[mid, ] <- moment * (1 + A) / A *
      exp(-outer(log_A[mid], a))
  }
  if (any(large)) {
    out[large, ] <- .loglogistic_ratio_large(a, log_A[large])
  }
  out
}

#' Mean of the primary event time on the window, as in Stan
#'
#' The uniform mean is used for |rho| w of at most 1e-3.
#'
#' @noRd
.loglogistic_primary_mean <- function(rho, pwindow) {
  if (is.null(rho) || abs(rho) * pwindow <= 1e-3) {
    return(pwindow / 2)
  }
  a <- abs(rho) * pwindow
  mean_up <- pwindow * (1 / -expm1(-a) - 1 / a)
  if (rho > 0) mean_up else pwindow - mean_up
}

#' Numerical primary event censored CDF of a log-logistic delay
#'
#' Integrates the delay CDF, or the survival where the delay CDF at the mean
#' primary event time is above a half.
#' Breaks are placed where the delay CDF moves so `integrate()` does not miss
#' a sharp step.
#'
#' @noRd
.loglogistic_numeric_cdf <- function(object, q, pwindow) {
  p <- .loglogistic_shape_scale(object)
  log_scale <- log(p$scale)
  primary <- function(x) {
    do.call(
      object$dprimary,
      c(list(x = x, min = 0, max = pwindow), object$primary_args)
    )
  }
  delay <- function(u, upper) {
    stats::plogis(
      p$shape * (log(pmax(u, 0)) - log_scale),
      lower.tail = !upper
    )
  }
  offsets <- p$scale * exp(seq(-30, 30) / p$shape)
  rho <- object$primary_args$r
  mean_primary <- .loglogistic_primary_mean(rho, pwindow)
  result <- vapply(q, function(d) {
    if (is.na(d) || d <= 0) {
      return(0)
    }
    breaks <- c(0, pwindow, if (d < pwindow) d, d - offsets)
    breaks <- sort(unique(breaks[is.finite(breaks)]))
    breaks <- breaks[breaks >= 0 & breaks <= pwindow]
    integral <- function(upper) {
      sum(vapply(
        seq_len(length(breaks) - 1L),
        function(i) {
          stats::integrate(
            function(x) delay(d - x, upper) * primary(x),
            lower = breaks[i], upper = breaks[i + 1L],
            rel.tol = 1e-13, abs.tol = 0, subdivisions = 1000L,
            stop.on.error = FALSE
          )$value
        },
        numeric(1)
      ))
    }
    # Survival where the delay CDF at the mean primary time is above a half
    if (d - mean_primary > p$scale) 1 - integral(TRUE) else integral(FALSE)
  }, numeric(1))
  pmin(1, pmax(0, result))
}

#' @rdname pcens_cdf_loglogistic
#' @export
pcens_cdf.pcens_pllogis_dunif <- function(
  object,
  q,
  pwindow,
  use_numeric = FALSE
) {
  p <- .loglogistic_shape_scale(object)
  if (isTRUE(use_numeric) || length(pwindow) != 1L || !is.finite(pwindow) ||
    pwindow <= 0 || p$shape < .loglogistic_floor_shape) {
    return(pcens_cdf.default(object, q, pwindow, use_numeric))
  }
  result <- rep(NA_real_, length(q))
  result[!is.na(q) & q == Inf] <- 1
  result[!is.na(q) & q == -Inf] <- 0
  finite <- which(is.finite(q))
  if (length(finite) > 0L) {
    uniform <- .loglogistic_uniform_terms(object, q[finite], pwindow)
    result[finite] <- uniform$cdf
    ill <- .loglogistic_uniform_ill(uniform)
    if (any(ill)) {
      result[finite[ill]] <- .loglogistic_numeric_cdf(
        object, q[finite[ill]], pwindow
      )
    }
  }
  result
}

#' Uniform CDF of a log-logistic delay and its error bound
#'
#' @noRd
.loglogistic_uniform_terms <- function(object, q, pwindow) {
  p <- .loglogistic_shape_scale(object)
  lower <- 0
  endpoints <- .exptilt_endpoints(q, pwindow, lower)
  log_g1 <- .loglogistic_log_moments(endpoints, p$shape, p$scale, 1L)[, 1L]
  at_q <- log_g1[.exptilt_index(q, endpoints, lower)]
  at_y <- log_g1[.exptilt_index(q - pwindow, endpoints, lower)]
  bounds <- .loglogistic_window_error(object, q, pwindow)
  list(
    cdf = pmin(1, exp(.log_diff_exp(at_q, at_y) - log(pwindow))),
    log_error = bounds$log_error,
    log_density = bounds$log_density
  )
}

#' Error bound of a difference of G_1 over a window
#'
#' The bound is 4e-16 (1 + |log tF|) tF / w at each endpoint, with the
#' density of the CDF for the PMF check.
#'
#' @noRd
.loglogistic_window_error <- function(object, q, pwindow) {
  p <- .loglogistic_shape_scale(object)
  log_tf <- function(t) {
    out <- rep(-Inf, length(t))
    positive <- t > 0
    out[positive] <- log(t[positive]) -
      .log1p_exp(-p$shape * (log(t[positive]) - log(p$scale)))
    out
  }
  log_cdf <- function(t) {
    out <- rep(-Inf, length(t))
    positive <- t > 0
    out[positive] <- -.log1p_exp(-p$shape * (log(t[positive]) - log(p$scale)))
    out
  }
  list(
    log_error = log(.loglogistic_error_constant) - log(pwindow) +
      .log_sum_exp(
        .log_error_scale(log_tf(q)), .log_error_scale(log_tf(q - pwindow))
      ),
    log_density = .log_diff_exp(log_cdf(q), log_cdf(q - pwindow)) -
      log(pwindow)
  )
}

#' Rounding error scale of a value held on the log scale
#'
#' @noRd
.log_error_scale <- function(x) {
  ifelse(is.finite(x), x + log1p(abs(x)), -Inf)
}

#' Check if the uniform solution is ill conditioned
#'
#' @noRd
.loglogistic_uniform_ill <- function(uniform) {
  .log_error_exceeds(
    uniform$log_error, log(uniform$cdf), uniform$log_density
  )
}

#' Test an error bound against the smaller tail and the PMF
#'
#' The bound is too large if it exceeds `.loglogistic_precision` times the
#' smaller of the CDF and the survival, not taken below 1e-15 in the upper
#' tail. With a density it must also be below `.loglogistic_pmf_precision`
#' times the density, not taken below 1e-13 times the smaller tail.
#'
#' @noRd
.log_error_exceeds <- function(log_error, log_cdf, log_density = NULL) {
  log_cdf <- pmin(log_cdf, 0)
  log_tail <- pmin(log_cdf, .log1m_exp(log_cdf))
  log_limit <- log(.loglogistic_precision) + log_tail
  upper <- log_cdf > -log(2)
  log_limit[upper] <- pmax(log_limit[upper], log(1e-15))
  if (!is.null(log_density)) {
    log_limit <- pmin(
      log_limit,
      pmax(
        log(.loglogistic_pmf_precision) + log_density,
        log(1e-13) + log_tail
      )
    )
  }
  ill <- log_error > log_limit
  ill[is.na(ill)] <- TRUE
  ill
}
