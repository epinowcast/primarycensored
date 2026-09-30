# References for the tests in test-stan-gamma-tail.R and
# test-stan-gamma-tail-gradient.R

# log F_T(t) for a Gamma delay via R
ref_lgamma_delay <- function(t, shape, rate) {
  pgamma(t * rate, shape = shape, log.p = TRUE)
}

# log of the uniform primary event censored CDF, by integrating F_T over
# [max(d - pwindow, 0), d] scaled by its maximum so that it does not
# underflow. The range is cut where the integrand is below exp(-40).
ref_lcdf_unif_gamma <- function(d, pwindow, shape, rate) {
  q_lo <- max(d - pwindow, 0)
  log_max <- ref_lgamma_delay(d, shape, rate)
  target <- log_max - 40
  lower <- q_lo
  if (ref_lgamma_delay(q_lo, shape, rate) < target) {
    lower <- uniroot(
      function(u) ref_lgamma_delay(u, shape, rate) - target,
      lower = q_lo, upper = d, tol = 1e-14
    )$root
  }
  scaled <- integrate(
    function(u) exp(ref_lgamma_delay(u, shape, rate) - log_max),
    lower = lower, upper = d, rel.tol = 1e-13, subdivisions = 1000L
  )$value
  log_max + log(scaled) - log(pwindow)
}

# log of the truncated uniform primary event censored CDF,
# (F(d) - F(L)) / (F(D) - F(L)), where F is the untruncated CDF
ref_lcdf_unif_gamma_trunc <- function(d, pwindow, shape, rate, L = 0,
                                      D = Inf) {
  log_diff <- function(a, b) a + log1p(-exp(b - a))
  lcdf <- function(x) ref_lcdf_unif_gamma(x, pwindow, shape, rate)
  log_lower <- if (L > 0) lcdf(L) else -Inf
  numerator <- if (L > 0) log_diff(lcdf(d), log_lower) else lcdf(d)
  if (is.infinite(D)) {
    return(numerator)
  }
  denominator <- if (L > 0) log_diff(lcdf(D), log_lower) else lcdf(D)
  numerator - denominator
}

# Gradient of `ref_lcdf_unif_gamma_trunc()` with respect to the log of
# (shape, rate), plus 1 for the log Jacobian of a lower bound of 0, by a
# fifth order central difference
ref_gamma_delay_gradient <- function(d, params, pwindow, L = 0, D = Inf) {
  f <- function(shape, rate) {
    ref_lcdf_unif_gamma_trunc(d, pwindow, shape, rate, L, D)
  }
  vapply(seq_along(params), function(i) {
    h <- 1e-3
    g <- function(e) {
      p <- params
      p[i] <- p[i] * exp(e)
      f(p[1], p[2])
    }
    (-g(2 * h) + 8 * g(h) - 8 * g(-h) + g(-2 * h)) / (12 * h) + 1
  }, numeric(1))
}
