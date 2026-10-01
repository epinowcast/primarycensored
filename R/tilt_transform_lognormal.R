#' Tilt transform of a lognormal delay
#'
#' The truncated exponential-moment transform of the lognormal has no closed
#' form. With \eqn{z_t = (\log t - \mu) / \sigma},
#' \deqn{T_f(\xi; t) = \int_{-\infty}^{z_t} e^{\xi e^{\mu + \sigma z}}
#'   \phi(z) dz.}
#' * For \eqn{\xi = 0} it is the lognormal CDF.
#' * For \eqn{\xi < 0} the integrand is log-concave and is integrated by
#'   Gauss-Legendre quadrature on panels around its mode. The upper
#'   transform is integrated in the same way, so it does not cancel.
#' * For \eqn{\xi > 0} the transform is the series
#'   \eqn{\sum_k \xi^k m_k(t) / k!} of positive terms in the partial moments
#'   \eqn{m_k(t) = e^{k \mu + k^2 \sigma^2 / 2} \Phi(z_t - k \sigma)}.
#'   The total diverges, so the upper transform is `Inf`.
#'
#' The series needs about \eqn{\xi t} terms per point. `.pcens_tilt_fits()`
#' is `FALSE` beyond \eqn{\xi t} of 200, unless \eqn{\xi w} is above 2, and
#' beyond 20000 terms.
#'
#' @inheritParams tilt_transform
#'
#' @param meanlog,sdlog Mean and standard deviation of the log of the delay.
#'
#' @return
#' * `.lnorm_meanlog_sdlog()`: a list with `meanlog` and `sdlog`, defaulting
#'   as in [stats::plnorm()].
#' * `.lnorm_tilt_pair()`: a matrix with columns `lower` and `upper`, the log
#'   of the transform over the lower and the upper part of the support.
#'
#' @keywords internal
#' @name tilt_transform_lognormal
.lnorm_meanlog_sdlog <- function(object) {
  mu <- object$args$meanlog
  sigma <- object$args$sdlog
  list(
    meanlog = if (is.null(mu)) 0 else mu,
    sdlog = if (is.null(sigma)) 1 else sigma
  )
}

# Gauss-Legendre rule shared by all panels, with as many nodes as Stan
.lnorm_rule_size <- 32L
.lnorm_cache <- new.env(parent = emptyenv())

# Gauss-Legendre nodes `x`, weights `w` and `log_w` on [-1, 1], cached
.lnorm_rule <- function(n = .lnorm_rule_size) {
  key <- as.character(n)
  if (!is.null(.lnorm_cache[[key]])) {
    return(.lnorm_cache[[key]])
  }
  # Eigenvalues of the Jacobi matrix, polished by Newton
  i <- seq_len(n - 1L)
  jacobi <- matrix(0, n, n)
  off <- i / sqrt(4 * i^2 - 1)
  jacobi[cbind(i, i + 1L)] <- off
  jacobi[cbind(i + 1L, i)] <- off
  x <- sort(eigen(jacobi, symmetric = TRUE, only.values = TRUE)$values)
  legendre <- function(x) {
    p_prev <- rep(1, n)
    p <- x
    for (j in 2:n) {
      p_next <- ((2 * j - 1) * x * p - (j - 1) * p_prev) / j
      p_prev <- p
      p <- p_next
    }
    list(p = p, dp = n * (x * p - p_prev) / (x^2 - 1))
  }
  for (iter in seq_len(3)) {
    lg <- legendre(x)
    x <- x - lg$p / lg$dp
  }
  lg <- legendre(x)
  w <- 2 / ((1 - x^2) * lg$dp^2)
  rule <- list(x = x, w = w, log_w = log(w))
  .lnorm_cache[[key]] <- rule
  rule
}

# Principal branch of the Lambert W function for x >= 0, by Halley's method
.lambert_w0 <- function(x) {
  if (x == 0) {
    return(0)
  }
  w <- log1p(x)
  for (iter in seq_len(20)) {
    ew <- exp(w)
    f <- w * ew - x
    delta <- f / (ew * (w + 1) - (w + 2) * f / (2 * w + 2))
    w <- w - delta
    if (abs(delta) <= 1e-16 * abs(w)) {
      break
    }
  }
  w
}

# Log of the integral of exp(-rho exp(meanlog + sdlog z)) dnorm(z) over each
# range (a, b), by Gauss-Legendre quadrature on `.lnorm_n_panels()` panels
.lnorm_panel <- function(a, b, meanlog, sdlog, rho, rule,
                         n_panels = .lnorm_n_panels(sdlog)) {
  if (n_panels == 1L) {
    return(.lnorm_one_panel(a, b, meanlog, sdlog, rho, rule))
  }
  width <- (b - a) / n_panels
  Reduce(.log_sum_exp, lapply(seq_len(n_panels), function(j) {
    .lnorm_one_panel(
      a + (j - 1L) * width, a + j * width, meanlog, sdlog, rho, rule
    )
  }))
}

.lnorm_one_panel <- function(a, b, meanlog, sdlog, rho, rule) {
  half <- (b - a) / 2
  z <- outer(half, rule$x) + (a + b) / 2
  log_f <- sweep(
    -rho * exp(meanlog + sdlog * z) - 0.5 * z^2, 2L, rule$log_w, "+"
  )
  peak <- log_f[cbind(seq_along(a), max.col(log_f, ties.method = "first"))]
  peak + log(rowSums(exp(log_f - peak))) + log(half) - 0.5 * log(2 * pi)
}

# One 32 point panel loses accuracy for sdlog above 1.8, so each range is
# split into ceiling(sdlog / 1.8) panels
.lnorm_n_panels <- function(sdlog) {
  max(1L, as.integer(ceiling(sdlog / 1.8)))
}

# Mode z0 = -W_0(rho sdlog^2 exp(meanlog)) / sdlog of the concave log
# integrand, limits lo and hi beyond which it is below exp(-40) of its peak,
# and the log total, for rho = -xi > 0
.lnorm_bump <- function(meanlog, sdlog, rho, rule) {
  w0 <- .lambert_w0(exp(log(rho) + 2 * log(sdlog) + meanlog))
  z0 <- -w0 / sdlog
  u0 <- exp(meanlog - w0)
  hi <- min(
    z0 + sqrt(80 / (1 + w0)),
    max((log(u0 + 40 / rho) - meanlog) / sdlog, -z0)
  )
  lo <- max(z0 - sqrt(80), -sqrt(z0^2 + 2 * w0 / sdlog^2 + 80))
  panels <- .lnorm_panel(
    c(lo, z0), c(z0, hi), meanlog, sdlog, rho, rule
  )
  list(z0 = z0, lo = lo, hi = hi, total = .log_sum_exp(panels[1], panels[2]))
}

# Lower and upper transform for xi < 0. The side of the mode that t is on is
# one panel ending at t. The other is the difference from the total, which
# does not cancel as it holds the mass on the far side of the mode.
.lnorm_tilt_quadrature <- function(t, meanlog, sdlog, xi,
                                   rule = .lnorm_rule()) {
  rho <- -xi
  bump <- .lnorm_bump(meanlog, sdlog, rho, rule)
  n <- length(t)
  lower <- rep(-Inf, n)
  upper <- rep(bump$total, n)
  positive <- which(t > 0)
  if (length(positive) == 0L) {
    return(cbind(lower = lower, upper = upper))
  }
  tp <- t[positive]
  z <- (log(tp) - meanlog) / sdlog
  # Slope of the log integrand at the point, positive left of the mode
  slope <- -rho * sdlog * tp - z
  left <- z < bump$z0
  curvature <- 1 + rho * sdlog^2 * tp
  width <- ifelse(
    left,
    80 / (slope + sqrt(slope^2 + 80)),
    80 / (abs(slope) + sqrt(slope^2 + 80 * curvature))
  )
  # Right of the mode, stop where the tilt alone has decayed by 40
  width <- ifelse(
    left, width,
    pmin(width, pmax(abs(z), (log(tp + 40 / rho) - meanlog) / sdlog) - z)
  )
  panel <- .lnorm_panel(
    ifelse(left, z - width, z), ifelse(left, z, z + width),
    meanlog, sdlog, rho, rule
  )
  rest <- .log_diff_exp(bump$total, panel)
  lower[positive] <- ifelse(left, panel, rest)
  upper[positive] <- ifelse(left, rest, panel)
  cbind(lower = lower, upper = upper)
}

# Log lower transform for xi > 0 as the sum of positive terms
# xi^k m_k(t) / k!. The number of terms doubles up to `.lnorm_max_terms`
# until the last is below exp(-40) of the largest. Points are taken in blocks
# of at most `cells` terms to bound the memory.
.lnorm_tilt_series <- function(t, meanlog, sdlog, xi, cells = 1e6) {
  out <- rep(-Inf, length(t))
  positive <- which(t > 0)
  if (length(positive) == 0L) {
    return(out)
  }
  z <- (log(t[positive]) - meanlog) / sdlog
  n_terms <- .lnorm_series_terms(xi, max(t[positive]))
  out[positive] <- .lnorm_series_blocks(z, meanlog, sdlog, xi, n_terms, cells)
  out
}

.lnorm_series_blocks <- function(z, meanlog, sdlog, xi, n_terms, cells) {
  if (n_terms > .lnorm_max_terms) {
    stop(
      "The lognormal tilt transform needs more than ", .lnorm_max_terms,
      " terms. Use use_numeric = TRUE.",
      call. = FALSE
    )
  }
  size <- max(1L, cells %/% (n_terms + 1L))
  out <- numeric(length(z))
  for (start in seq(1L, length(z), by = size)) {
    i <- start:min(start + size - 1L, length(z))
    out[i] <- .lnorm_series_block(z[i], meanlog, sdlog, xi, n_terms, cells)
  }
  out
}

.lnorm_series_block <- function(z, meanlog, sdlog, xi, n_terms, cells) {
  k <- seq_len(n_terms + 1L) - 1L
  log_terms <- .lnorm_log_pnorm(z, k * sdlog) +
    rep(
      k * (log(xi) + meanlog) + 0.5 * k^2 * sdlog^2 - lgamma(k + 1),
      each = length(z)
    )
  peak <- log_terms[cbind(
    seq_along(z), max.col(log_terms, ties.method = "first")
  )]
  if (all(log_terms[, n_terms + 1L] - peak < -40)) {
    return(peak + log(rowSums(exp(log_terms - peak))))
  }
  more <- if (n_terms >= .lnorm_max_terms) {
    .lnorm_max_terms + 1L
  } else {
    min(2 * n_terms, .lnorm_max_terms)
  }
  .lnorm_series_blocks(z, meanlog, sdlog, xi, more, cells)
}

.lnorm_max_terms <- 20000L

# Largest xi t for the series, a speed cut-off against the numerical method
# set separately from the Stan cut-off, which is set for the ODE
.lnorm_series_max_xt <- 200

# The numerical method loses accuracy in the lower tail where xi w is above 2,
# so the series is kept there up to `.lnorm_max_terms`
.lnorm_series_min_xw <- 2

# Number of terms the series starts with, see `.lnorm_tilt_series()`
.lnorm_series_terms <- function(xi, t) {
  lambda <- xi * pmax(t, 0)
  ceiling(lambda + 9 * sqrt(lambda) + 30)
}

# The transform has a fixed cost that the numerical method beats for fewer
# than about 10 quantiles, where |r| w is at most 1 and it is accurate
.lnorm_exptilt_min_q <- 10L
.lnorm_exptilt_min_xw <- 1

# Matrix of log(Phi(z_i - shift_j))
.lnorm_log_pnorm <- function(z, shift) {
  stats::pnorm(outer(z, shift, "-"), log.p = TRUE)
}

#' @rdname tilt_transform_lognormal
.lnorm_tilt_pair <- function(t, meanlog, sdlog, xi) {
  if (xi == 0) {
    z <- (log(pmax(t, 0)) - meanlog) / sdlog
    return(cbind(
      lower = stats::pnorm(z, log.p = TRUE),
      upper = stats::pnorm(z, lower.tail = FALSE, log.p = TRUE)
    ))
  }
  if (xi > 0) {
    return(cbind(
      lower = .lnorm_tilt_series(t, meanlog, sdlog, xi),
      upper = rep(Inf, length(t))
    ))
  }
  .lnorm_tilt_quadrature(t, meanlog, sdlog, xi)
}

#' @rdname tilt_transform
#' @exportS3Method
.pcens_tilt_lower.pcens_plnorm <- function(object) {
  0
}

#' @rdname tilt_transform
#' @exportS3Method
.pcens_tilt_available.pcens_plnorm <- function(object, xi) {
  p <- .lnorm_meanlog_sdlog(object)
  if (!is.finite(p$meanlog) || !is.finite(p$sdlog) || p$sdlog <= 0) {
    return(FALSE)
  }
  xi >= 0 || log(-xi) + 2 * log(p$sdlog) + p$meanlog < 690
}

#' @rdname tilt_transform
#' @exportS3Method
.pcens_tilt_fits.pcens_plnorm <- function(object, xi, t, pwindow = 0) {
  if (xi <= 0) {
    return(rep(TRUE, length(t)))
  }
  (xi * pmax(t, 0) <= .lnorm_series_max_xt |
    xi * pwindow > .lnorm_series_min_xw) &
    .lnorm_series_terms(xi, t) <= .lnorm_max_terms
}

#' @rdname tilt_transform
#' @exportS3Method
.pcens_tilt_transform.pcens_plnorm <- function(object, t, xi, upper = FALSE) {
  p <- .lnorm_meanlog_sdlog(object)
  .lnorm_tilt_pair(t, p$meanlog, p$sdlog, xi)[, if (upper) 2L else 1L]
}

#' @rdname tilt_transform
#' @exportS3Method
.pcens_tilt_pair.pcens_plnorm <- function(object, t, xi) {
  p <- .lnorm_meanlog_sdlog(object)
  .lnorm_tilt_pair(t, p$meanlog, p$sdlog, xi)
}

#' @rdname tilt_transform
#' @exportS3Method
.pcens_tilt_moments.pcens_plnorm <- function(object, t) {
  # G_1 = t m_0 - m_1 and G_2 = t G_1 - (t m_1 - m_2) from the partial
  # moments m_k, as differences of positive integrals
  p <- .lnorm_meanlog_sdlog(object)
  positive <- t > 0
  tp <- pmax(t, 0)
  log_t <- log(tp)
  z <- (log_t - p$meanlog) / p$sdlog
  log_m0 <- stats::pnorm(z, log.p = TRUE)
  log_m1 <- p$meanlog + 0.5 * p$sdlog^2 +
    stats::pnorm(z - p$sdlog, log.p = TRUE)
  log_m2 <- 2 * p$meanlog + 2 * p$sdlog^2 +
    stats::pnorm(z - 2 * p$sdlog, log.p = TRUE)
  log_g1 <- .log_diff_exp(log_t + log_m0, log_m1)
  log_h <- .log_diff_exp(log_t + log_m1, log_m2)
  log_g2 <- .log_diff_exp(log_t + log_g1, log_h)
  cbind(
    G1 = ifelse(positive, log_g1, -Inf),
    G2 = ifelse(positive, log_g2, -Inf)
  )
}
