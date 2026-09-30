# Helpers for the exponentially tilted primary event window tests.
#
# The reference integrates the delay CDF against the tilted window density
# with tight tolerances, to about 1e-12 relative.

# Density of the tilted window on [0, w] written to avoid overflow and
# cancellation: rho exp(rho z) / (exp(rho w) - 1).
exptilt_window_density <- function(z, w, rho) {
  if (abs(rho * w) < 1e-12) {
    return((1 + rho * (z - w / 2)) / w)
  }
  if (rho > 0) {
    rho * exp(rho * (z - w)) / -expm1(-rho * w)
  } else {
    rho * exp(rho * z) / expm1(rho * w)
  }
}

# Reference primary event censored CDF for any delay CDF `cdf(x)`. The
# integral is split at the kink where the delay CDF leaves zero.
exptilt_reference <- function(q, pwindow, rho, cdf) {
  vapply(q, function(qq) {
    integrand <- function(z) {
      cdf(qq - z) * exptilt_window_density(z, pwindow, rho)
    }
    breaks <- c(0, if (qq > 0 && qq < pwindow) qq, pwindow)
    sum(vapply(seq_len(length(breaks) - 1L), function(i) {
      stats::integrate(
        integrand, breaks[i], breaks[i + 1L],
        rel.tol = 1e-13, abs.tol = 0, subdivisions = 2000L
      )$value
    }, numeric(1)))
  }, numeric(1))
}

# Delay families with a closed form transform, and a parameter set per
# family covering small and large shapes and negative means.
exptilt_families <- function() {
  list(
    list(
      label = "exponential rate 2", pdist = pexp, rdist = rexp,
      args = list(rate = 2), rate = 2, positive = TRUE
    ),
    list(
      label = "exponential rate 0.3", pdist = pexp, rdist = rexp,
      args = list(rate = 0.3), rate = 0.3, positive = TRUE
    ),
    list(
      label = "gamma shape 0.6", pdist = pgamma, rdist = rgamma,
      args = list(shape = 0.6, rate = 1.3), rate = 1.3, positive = TRUE
    ),
    list(
      label = "gamma shape 2.5", pdist = pgamma, rdist = rgamma,
      args = list(shape = 2.5, scale = 2.5), rate = 0.4, positive = TRUE
    ),
    list(
      label = "gamma shape 20", pdist = pgamma, rdist = rgamma,
      args = list(shape = 20, rate = 4), rate = 4, positive = TRUE
    ),
    list(
      label = "normal mean 3", pdist = pnorm, rdist = rnorm,
      args = list(mean = 3, sd = 2), rate = Inf, positive = FALSE
    ),
    list(
      label = "normal mean -1", pdist = pnorm, rdist = rnorm,
      args = list(mean = -1, sd = 3), rate = Inf, positive = FALSE
    )
  )
}

# Tilts for which the exponential and gamma forms are admissible
# (`rate + rho > 0`) for a family.
exptilt_admissible <- function(family, rho) {
  family$rate + rho > 0
}

# Delay CDF of a family as a function of the quantile alone.
exptilt_cdf <- function(family) {
  function(x) do.call(family$pdist, c(list(x), family$args))
}

exptilt_object <- function(family, rho) {
  do.call(
    new_pcens,
    c(
      list(
        pdist = family$pdist, dprimary = dexpgrowth,
        primary_args = list(r = rho)
      ),
      family$args
    )
  )
}

# Label for a failed expectation.
exptilt_label <- function(family, pwindow, rho) {
  sprintf("%s, pwindow = %g, r = %g", family$label, pwindow, rho)
}

# Largest relative difference, ignoring values that underflow.
max_rel_diff <- function(actual, expected) {
  keep <- expected > 1e-300
  max(abs(actual[keep] / expected[keep] - 1))
}
