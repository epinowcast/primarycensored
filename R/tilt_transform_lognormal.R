#' Tilt transform of a lognormal delay
#'
#' The truncated exponential-moment transform of the lognormal has no closed
#' form, so it is evaluated by two methods that depend on the sign of the tilt
#' \eqn{\xi}. With \eqn{u = e^{\mu + \sigma z}} and \eqn{z} standard normal,
#' \deqn{T_f(\xi; t) = \int_{-\infty}^{z_t} e^{\xi e^{\mu + \sigma z}}
#'   \phi(z) dz, \qquad z_t = (\log t - \mu) / \sigma.}
#' * For \eqn{\xi = 0} it is the lognormal CDF.
#' * For \eqn{\xi < 0}, the tilted primary with \eqn{\rho = -\xi > 0}, the
#'   integrand is log-concave and is integrated by Gauss-Legendre quadrature
#'   on panels chosen from its mode, see `.lnorm_tilt_quadrature()`.
#'   The tail transform \eqn{T_f(\xi; \infty) - T_f(\xi; t)} is integrated in
#'   the same way, so it does not cancel.
#' * For \eqn{\xi > 0}, the tilted primary with \eqn{\rho < 0}, the transform
#'   is the series \eqn{\sum_k \xi^k m_k(t) / k!} of positive terms with the
#'   partial moments \eqn{m_k(t) = e^{k \mu + k^2 \sigma^2 / 2}
#'   \Phi(z_t - k \sigma)}, see `.lnorm_tilt_series()`.
#'   The total diverges, so the tail transform is `Inf` on the log scale.
#'   The series costs about \eqn{\xi t} terms per point. `.pcens_tilt_fits()`
#'   is `FALSE` beyond \eqn{\xi t} of 200, where the numerical method is
#'   faster, unless the window is wide in tilt terms, see
#'   `.lnorm_series_max_xt`. It is `FALSE` beyond 20000 terms, about
#'   \eqn{\xi t} of 18700, whatever the window.
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

# Gauss-Legendre rule shared by all panels. The Stan functions use the same
# number of nodes, see `primarycensored_lognormal_tilt_rule()`. The nodes are
# found by Newton's method on the Legendre polynomial and cached.
.lnorm_rule_size <- 32L
.lnorm_cache <- new.env(parent = emptyenv())

#' Gauss-Legendre rule on \[-1, 1\]
#'
#' @param n Number of nodes.
#'
#' @return A list with the nodes `x`, the weights `w` and their log `log_w`.
#'
#' @keywords internal
.lnorm_rule <- function(n = .lnorm_rule_size) {
  key <- as.character(n)
  if (!is.null(.lnorm_cache[[key]])) {
    return(.lnorm_cache[[key]])
  }
  # The nodes are the eigenvalues of the Jacobi matrix, polished by Newton
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

#' Principal branch of the Lambert W function
#'
#' Solves \eqn{w e^w = x} for \eqn{x \ge 0} by Halley's method.
#'
#' @param x Non-negative number.
#'
#' @return \eqn{W_0(x)}.
#'
#' @keywords internal
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

#' Log of a Gaussian weighted integral over panels
#'
#' The log of
#' \eqn{\int_a^b e^{-\rho e^{\mu + \sigma z}} \phi(z) dz}
#' for each pair of limits by Gauss-Legendre quadrature, where \eqn{\phi} is
#' the standard normal density. The integrand falls from its plateau to zero
#' over a width of about \eqn{1 / \sigma}, so the range of each pair is split
#' into `.lnorm_n_panels()` equal panels, one for \eqn{\sigma} up to 1.8.
#'
#' @param a,b Numeric vectors of lower and upper limits of equal length.
#'
#' @inheritParams tilt_transform_lognormal
#'
#' @param rho Tilt, \eqn{\rho = -\xi}, at least 0.
#'
#' @param rule Output of `.lnorm_rule()`.
#'
#' @param n_panels Number of equal panels each range is split into.
#'
#' @return Numeric vector of the log integrals.
#'
#' @keywords internal
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

#' One Gauss-Legendre panel of a Gaussian weighted integral
#'
#' @inheritParams .lnorm_panel
#'
#' @return Numeric vector of the log integrals.
#'
#' @keywords internal
.lnorm_one_panel <- function(a, b, meanlog, sdlog, rho, rule) {
  half <- (b - a) / 2
  z <- outer(half, rule$x) + (a + b) / 2
  log_f <- sweep(
    -rho * exp(meanlog + sdlog * z) - 0.5 * z^2, 2L, rule$log_w, "+"
  )
  peak <- log_f[cbind(seq_along(a), max.col(log_f, ties.method = "first"))]
  peak + log(rowSums(exp(log_f - peak))) + log(half) - 0.5 * log(2 * pi)
}

#' Number of panels for a lognormal tilt integral
#'
#' The 32 point rule on one panel per side is accurate to about 1e-11 in the
#' log transform for `sdlog` up to 1.8, and loses accuracy beyond, to 1e-7 in
#' the CDF at 4 and 3e-5 at 15. Splitting each range into
#' \eqn{\lceil \sigma / 1.8 \rceil} panels keeps the width of a panel
#' in units of \eqn{1 / \sigma} within the range that was tested.
#'
#' @inheritParams tilt_transform_lognormal
#'
#' @return Integer number of panels, at least 1.
#'
#' @keywords internal
.lnorm_n_panels <- function(sdlog) {
  max(1L, as.integer(ceiling(sdlog / 1.8)))
}

#' Mode and total of the tilted lognormal integrand
#'
#' With \eqn{\rho = -\xi > 0} the log of the integrand is
#' \eqn{\ell(z) = -\rho e^{\mu + \sigma z} - z^2 / 2}. It is concave with
#' \eqn{\ell'' \le -1}, so it has one mode \eqn{z_0 = -W_0(\rho \sigma^2
#' e^\mu) / \sigma} and \eqn{\ell'' (z_0) = -(1 + W_0)}. Beyond the limits
#' below the integrand is less than \eqn{e^{-40}} of its peak:
#' * Left, from the two bounds \eqn{\ell(z) - \ell(z_0) \le -(z - z_0)^2 / 2}
#'   and \eqn{\le -(z^2 - z_0^2) / 2 + \rho e^{\mu + \sigma z_0}}.
#' * Right, from the curvature at the mode and from
#'   \eqn{\rho \{u(z) - u(z_0)\} > 40} for \eqn{z \ge |z_0|}.
#'
#' The total over the whole line is the sum over the two panels that meet at
#' the mode. The mode and the total do not depend on the point, so they are
#' computed once for all points.
#'
#' @inheritParams .lnorm_panel
#'
#' @return A list with the mode `z0`, the limits `lo` and `hi`, and the log
#'   `total` integral over `[lo, hi]`.
#'
#' @keywords internal
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

#' Tilted lognormal transform by quadrature
#'
#' The lower and the upper transform for a tilt \eqn{\xi < 0}. Left of the
#' mode the lower transform is one panel that ends at the point, bounded on
#' the left from the slope of the integrand there. Right of the mode the
#' upper transform is one panel that starts at the point, bounded on the
#' right from the slope and the curvature there. The other transform is the
#' difference from the total of the bump, see `.lnorm_bump()`. It is at least
#' the mass on the other side of the mode, which is not small, so the
#' difference does not cancel.
#'
#' The panels are accurate to an absolute difference of the log transform of
#' about 1e-11 for sdlog up to 1.8 and of about 1e-13 for sdlog of 1 or
#' below. Each range is split into more panels for a larger sdlog, see
#' `.lnorm_n_panels()`, which keeps the same accuracy to an sdlog of 15.
#'
#' @inheritParams .lnorm_panel
#'
#' @param t Numeric vector of finite points.
#'
#' @param xi Tilt, negative.
#'
#' @return A matrix with columns `lower` and `upper`, the log of the
#'   transform over the lower and the upper part of the support.
#'
#' @keywords internal
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
  # Right of the mode the integrand is below e^-40 of its value at the point
  # beyond the point where the tilt alone has decayed by 40
  width <- ifelse(
    left, width,
    pmin(width, pmax(abs(z), (log(tp + 40 / rho) - meanlog) / sdlog) - z)
  )
  panel <- .lnorm_panel(
    ifelse(left, z - width, z), ifelse(left, z, z + width),
    meanlog, sdlog, rho, rule
  )
  # The panel is the lower transform left of the mode and the upper
  # transform right of it. The other one is the rest of the total.
  rest <- .log_diff_exp(bump$total, panel)
  lower[positive] <- ifelse(left, panel, rest)
  upper[positive] <- ifelse(left, rest, panel)
  cbind(lower = lower, upper = upper)
}

#' Tilted lognormal transform by series
#'
#' The lower transform for a tilt \eqn{\xi > 0} as the sum of the positive
#' terms \eqn{\xi^k m_k(t) / k!} over the partial moments
#' \eqn{m_k(t) = \int_0^t u^k f(u) du =
#' e^{k \mu + k^2 \sigma^2 / 2} \Phi(z_t - k \sigma)}.
#' Past \eqn{k = \xi t} the terms fall by at least a factor
#' \eqn{\xi t / (k + 1)} each. The number of terms starts at about
#' \eqn{\xi t + 9 \sqrt{\xi t} + 30} and doubles, up to `.lnorm_max_terms`,
#' until the last term is below \eqn{e^{-40}} of the largest for every point.
#' It stops with an error where that needs more than `.lnorm_max_terms`, so
#' callers check `.pcens_tilt_fits()` first, which also applies the cut-off
#' of `.lnorm_series_max_xt`.
#'
#' @inheritParams .lnorm_tilt_quadrature
#'
#' @param xi Tilt, positive.
#'
#' @return Numeric vector of the log transform, `-Inf` for `t <= 0`.
#'
#' @keywords internal
.lnorm_tilt_series <- function(t, meanlog, sdlog, xi) {
  out <- rep(-Inf, length(t))
  positive <- which(t > 0)
  if (length(positive) == 0L) {
    return(out)
  }
  z <- (log(t[positive]) - meanlog) / sdlog
  # The terms fall below e^-40 of the largest within about 9 sqrt(xi t)
  # of xi t, which the check below confirms
  n_terms <- .lnorm_series_terms(xi, max(t[positive]))
  too_long <- function() {
    stop(
      "The lognormal tilt transform needs more than ", .lnorm_max_terms,
      " terms. Use use_numeric = TRUE.",
      call. = FALSE
    )
  }
  if (n_terms > .lnorm_max_terms) {
    too_long()
  }
  repeat {
    k <- seq_len(n_terms + 1L) - 1L
    log_terms <- .lnorm_log_pnorm(z, k * sdlog) +
      rep(
        k * (log(xi) + meanlog) + 0.5 * k^2 * sdlog^2 - lgamma(k + 1),
        each = length(z)
      )
    peak <- log_terms[cbind(
      seq_along(z), max.col(log_terms, ties.method = "first")
    )]
    # Terms of the last block relative to the largest
    last <- log_terms[, n_terms + 1L] - peak
    if (all(last < -40)) {
      break
    }
    if (n_terms >= .lnorm_max_terms) {
      too_long()
    }
    n_terms <- min(2 * n_terms, .lnorm_max_terms)
  }
  out[positive] <- peak + log(rowSums(exp(log_terms - peak)))
  out
}

.lnorm_max_terms <- 20000L

# The series costs about xi t terms per quantile, about 0.05 microseconds
# each, and the numerical method of `pcens_cdf.default()` a fixed cost of
# about 20 microseconds per quantile. They cross at xi t of about 200 for
# 12 to 200 quantiles (ratio 1.0 at 200, 0.7 at 100, 1.4 at 300 and 2 at
# 500 for 200 quantiles), so the series is kept up to 200.
.lnorm_series_max_xt <- 200

# The numerical method integrates the delay CDF against the window density.
# Where xi w is above 2 the density falls by more than e^2 across the window
# and the integrator loses accuracy in the lower tail, where the CDF is
# small and varies fast. The relative error was 6e-4 to 5e-3 at xi w of 600
# for a CDF of 1e-7 to 1e-4, while the series stays accurate to 1e-13. The
# series is kept there up to `.lnorm_max_terms`.
.lnorm_series_min_xw <- 2

#' Number of terms the lognormal series starts with
#'
#' The terms \eqn{\xi^k m_k(t) / k!} are below \eqn{e^{-40}} of the largest
#' for \eqn{k} beyond \eqn{\xi t + 9 \sqrt{\xi t} + 30}, see
#' `.lnorm_tilt_series()`. The series is used only where this is at most
#' `.lnorm_max_terms`, which is for \eqn{\xi t} up to about 18700, and where
#' \eqn{\xi t} is at most `.lnorm_series_max_xt` or \eqn{\xi w} is above
#' `.lnorm_series_min_xw`.
#'
#' @param xi Tilt, positive.
#'
#' @param t Numeric vector of points.
#'
#' @return Numeric vector of the number of terms.
#'
#' @keywords internal
.lnorm_series_terms <- function(xi, t) {
  lambda <- xi * pmax(t, 0)
  ceiling(lambda + 9 * sqrt(lambda) + 30)
}

# The quadrature and the series have a fixed cost of about 0.2 ms that the
# numerical method of `pcens_cdf.default()` beats for fewer than about 10
# quantiles, and they are faster beyond that, see the benchmarks in NEWS.md.
# The numerical method is used for fewer quantiles only where it is accurate,
# which is for |r| w up to 1. Its relative error is about 1e-6 there, 1e-4 at
# 50, and it fails at 1000, where the transform is accurate to 1e-9.
.lnorm_exptilt_min_q <- 10L
.lnorm_exptilt_min_xw <- 1

#' Log standard normal CDF at every point less every shift
#'
#' @param z,shift Numeric vectors.
#'
#' @return A matrix with `length(z)` rows and `length(shift)` columns of
#'   `log(Phi(z_i - shift_j))`.
#'
#' @keywords internal
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
  # The mode of the integrand needs rho sdlog^2 exp(meanlog) to be finite
  xi >= 0 || log(-xi) + 2 * log(p$sdlog) + p$meanlog < 690
}

#' @rdname tilt_transform
#' @exportS3Method
.pcens_tilt_fits.pcens_plnorm <- function(object, xi, t, pwindow = 0) {
  # The series for a positive tilt is used where it is faster than the
  # numerical method or the numerical method is not accurate, up to a limit
  # on the number of terms
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
  # The partial moments m_k(t) = e^{k mu + k^2 sigma^2 / 2} Phi(z - k sigma)
  # give G_1 = t m_0 - m_1 and G_2 = t G_1 - (t m_1 - m_2), every difference
  # being of positive integrals as for the gamma
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
