#' Methods for log-logistic delays
#'
#' Analytical primary event censored CDFs for a log-logistic delay with a
#' uniform primary event window.
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
#' With \eqn{G_1(t) = \int_0^t F(u) du} the CDF is
#' \eqn{\{G_1(q) - G_1(q - w)\} / w}, evaluated with series for a shape of
#' at least 0.01.
#' The numerical method is used where this difference loses precision, so
#' the CDF can step by about 1e-8 of the smaller tail at the switch.
#' See `vignette("analytic-solutions")` for the derivation and the switching
#' rule.
#'
#' @inherit pcens_cdf return
#'
#' @name pcens_cdf_loglogistic
#'
#' @concept pcens
#'
#' @examplesIf requireNamespace("flexsurv", quietly = TRUE)
#' # Log-logistic delay with a uniform primary
#' obj <- new_pcens(
#'   pdist = flexsurv::pllogis, dprimary = dunif,
#'   primary_args = list(), shape = 2, scale = 5
#' )
#' pcens_cdf(obj, q = c(1, 4, 8), pwindow = 2)
NULL

# Smallest shape for which the solution is used, as the ratios involve
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

# Relative size of the last series term, and the most terms summed, as in
# the Stan functions. The tail series converges geometrically from A = 3, so
# needs fewer terms.
.loglogistic_series_tol <- 1e-16
.loglogistic_max_terms <- 5000L
.loglogistic_tail_max_terms <- 200L

# Floors on the error limit, in the upper tail and relative to the smaller
# tail, near the resolution of a CDF close to 1
.loglogistic_upper_floor <- 1e-15
.loglogistic_tail_floor <- 1e-13

# Relative tolerance of the numerical method
.loglogistic_numeric_tol <- 1e-13

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
    if (all(term <= .loglogistic_series_tol * total) ||
      m >= .loglogistic_max_terms) {
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
    if (all(term <= .loglogistic_series_tol * total) ||
      j >= .loglogistic_max_terms) {
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
    if (all(abs(term) <= .loglogistic_series_tol * abs(total)) ||
      m >= .loglogistic_tail_max_terms) {
      break
    }
  }
  (1 + decay) * total
}

#' Partial moment ratio r_a(A) of the log-logistic
#'
#' Rows are `log_A` and columns are `a`. The three series are described in
#' `vignette("analytic-solutions")`.
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

#' Log of G_1(t) = int_0^t F(u) du for the delay
#'
#' @noRd
.loglogistic_log_g1 <- function(t, shape, scale) {
  out <- rep(-Inf, length(t))
  positive <- t > 0
  if (!any(positive)) {
    return(out)
  }
  log_t <- log(t[positive])
  log_ratio <- shape * (log_t - log(scale))
  r <- .loglogistic_ratio(1 / shape, log_ratio)[, 1L]
  out[positive] <- log_t + stats::plogis(log_ratio, log.p = TRUE) +
    log1p(-r)
  out
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
  result <- vapply(q, function(d) {
    if (is.na(d) || d <= 0) {
      return(0)
    }
    breaks <- c(0, pwindow, if (d < pwindow) d, d - offsets)
    breaks <- sort(unique(breaks[is.finite(breaks)]))
    breaks <- breaks[breaks >= 0 & breaks <= pwindow]
    # integrate() often reports roundoff at this tolerance while its value is
    # still accurate, so its message is not used. Its summed error estimate
    # is checked against the limit used for the analytical solution instead.
    integral <- function(upper) {
      pieces <- vapply(
        seq_len(length(breaks) - 1L),
        function(i) {
          fit <- stats::integrate(
            function(x) delay(d - x, upper) * primary(x),
            lower = breaks[i], upper = breaks[i + 1L],
            rel.tol = .loglogistic_numeric_tol, abs.tol = 0,
            subdivisions = 1000L, stop.on.error = FALSE
          )
          c(fit$value, fit$abs.error)
        },
        numeric(2)
      )
      value <- sum(pieces[1L, ])
      # A survival is resolved only to 1e-15 once taken from 1, and a CDF
      # below the smallest normal double has no relative precision
      min_limit <- if (upper) {
        .loglogistic_upper_floor
      } else {
        .Machine$double.xmin
      }
      limit <- max(.loglogistic_precision * value, min_limit)
      if (!isTRUE(sum(pieces[2L, ]) <= limit)) {
        warning(
          "Numerical integration of the log-logistic CDF may be inaccurate ",
          "at q = ", d, call. = FALSE
        )
      }
      value
    }
    # Survival where the delay CDF at the mean primary time is above a half
    if (d - pwindow / 2 > p$scale) 1 - integral(TRUE) else integral(FALSE)
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
    q_finite <- q[finite]
    cdf <- .loglogistic_uniform_cdf(q_finite, pwindow, p)
    ill <- .loglogistic_uniform_ill(q_finite, pwindow, p, cdf)
    if (any(ill)) {
      cdf[ill] <- .loglogistic_numeric_cdf(object, q_finite[ill], pwindow)
    }
    result[finite] <- cdf
  }
  result
}

#' Log CDF of the log-logistic delay, `-Inf` for t <= 0
#'
#' @noRd
.loglogistic_lcdf <- function(t, shape, scale) {
  out <- rep(-Inf, length(t))
  positive <- t > 0
  out[positive] <- stats::plogis(
    shape * (log(t[positive]) - log(scale)),
    log.p = TRUE
  )
  out
}

#' Uniform primary CDF of a log-logistic delay at finite `q`
#'
#' Each unique endpoint is evaluated once.
#'
#' @param p Shape and scale, from `.loglogistic_shape_scale()`.
#'
#' @noRd
.loglogistic_uniform_cdf <- function(q, pwindow, p) {
  pos <- .exptilt_layout(q, pwindow, 0)
  log_g1 <- .loglogistic_log_g1(pos$endpoints, p$shape, p$scale)
  pmin(1, exp(.log_diff_exp(log_g1[pos$q], log_g1[pos$y]) - log(pwindow)))
}

#' Check if the uniform solution is ill conditioned
#'
#' The error of the difference of G_1 is bounded by
#' 4e-16 (1 + |log tF|) tF / w at each endpoint.
#' The bound is too large if it exceeds `.loglogistic_precision` times the
#' smaller of the CDF and the survival, not taken below
#' `.loglogistic_upper_floor` in the upper
#' tail, or `.loglogistic_pmf_precision` times the density of the CDF, not
#' taken below `.loglogistic_tail_floor` times the smaller tail.
#'
#' @param cdf The uniform CDF at `q`.
#'
#' @noRd
.loglogistic_uniform_ill <- function(q, pwindow, p, cdf) {
  log_cdf_q <- .loglogistic_lcdf(q, p$shape, p$scale)
  log_cdf_y <- .loglogistic_lcdf(q - pwindow, p$shape, p$scale)
  # Rounding error scale of log(t F(t)), held on the log scale
  log_error_scale <- function(t, log_cdf) {
    x <- log(pmax(t, 0)) + log_cdf
    ifelse(is.finite(x), x + log1p(abs(x)), -Inf)
  }
  log_error <- log(.loglogistic_error_constant) - log(pwindow) +
    .log_sum_exp(
      log_error_scale(q, log_cdf_q), log_error_scale(q - pwindow, log_cdf_y)
    )
  log_density <- .log_diff_exp(log_cdf_q, log_cdf_y) - log(pwindow)
  log_cdf <- pmin(log(cdf), 0)
  log_tail <- pmin(log_cdf, .log1m_exp(log_cdf))
  log_limit <- log(.loglogistic_precision) + log_tail
  upper <- log_cdf > -log(2)
  log_limit[upper] <- pmax(log_limit[upper], log(.loglogistic_upper_floor))
  log_limit <- pmin(
    log_limit,
    pmax(
      log(.loglogistic_pmf_precision) + log_density,
      log(.loglogistic_tail_floor) + log_tail
    )
  )
  ill <- log_error > log_limit
  ill[is.na(ill)] <- TRUE
  ill
}
