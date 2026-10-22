# References for the tests in test-stan-gamma-tail.R,
# test-stan-gamma-tail-gradient.R and test-stan-gengamma-tail.R

# log F_T(t) for the generalised gamma (Stacy parameterisation) via R
ref_lgengamma <- function(t, shape, scale, k) {
  pgamma((t / scale)^shape, shape = k, log.p = TRUE)
}

# Log uniform primary censored CDF, (1 / w) int_{max(d - w, 0)}^{d} F_T(u) du.
# The integrand is scaled by its value at d and cut once it is below exp(-40).
ref_lcdf_unif <- function(d, pwindow, shape, scale, k) {
  q_lo <- max(d - pwindow, 0)
  log_max <- ref_lgengamma(d, shape, scale, k)
  target <- log_max - 40
  lower <- q_lo
  if (ref_lgengamma(q_lo, shape, scale, k) < target) {
    lower <- uniroot(
      function(u) ref_lgengamma(u, shape, scale, k) - target,
      lower = q_lo, upper = d, tol = 1e-14
    )$root
  }
  scaled <- integrate(
    function(u) exp(ref_lgengamma(u, shape, scale, k) - log_max),
    lower = lower, upper = d, rel.tol = 1e-13, subdivisions = 1000L
  )$value
  log_max + log(scaled) - log(pwindow)
}

# The same for a Gamma delay, the generalised gamma with shape 1
ref_lcdf_unif_gamma <- function(d, pwindow, shape, rate) {
  ref_lcdf_unif(d, pwindow, 1, 1 / rate, shape)
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
