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

# Log reference for a gamma delay in the lower tail, where the delay
# density varies over many orders of magnitude within the window. The
# integral over the window is scaled by the density at q and split at
# geometric steps from q.
exptilt_gamma_log_reference <- function(q, pwindow, rho, shape, rate) {
  vapply(q, function(qq) {
    scale <- stats::dgamma(qq, shape, rate, log = TRUE)
    top <- min(pwindow, qq)
    breaks <- sort(unique(c(0, top * 10^(-8:0))))
    integral <- sum(vapply(seq_len(length(breaks) - 1L), function(i) {
      stats::integrate(
        function(x) {
          exp(stats::dgamma(qq - x, shape, rate, log = TRUE) - scale) *
            expm1(rho * x) / expm1(rho * pwindow)
        },
        breaks[i], breaks[i + 1L],
        rel.tol = 1e-13, abs.tol = 0, subdivisions = 2000L
      )$value
    }, numeric(1)))
    below <- stats::pgamma(qq - pwindow, shape, rate, log.p = TRUE) - scale
    scale + log(integral + exp(below))
  }, numeric(1))
}

# Log reference for a normal delay, scaled to keep a CDF far below the
# smallest double.
exptilt_normal_log_reference <- function(d, pwindow, rho, mu, sigma) {
  vapply(d, function(dd) {
    log_integrand <- function(z) {
      stats::pnorm((dd - z - mu) / sigma, log.p = TRUE) +
        log(exptilt_window_density(z, pwindow, rho))
    }
    shift <- max(log_integrand(0), log_integrand(pwindow))
    integral <- stats::integrate(
      function(z) exp(vapply(z, log_integrand, numeric(1)) - shift),
      0, pwindow,
      rel.tol = 1e-13, abs.tol = 0, subdivisions = 2000L
    )$value
    shift + log(integral)
  }, numeric(1))
}

# Normal delay with mean -4 and sd 0.3 from 7 to 37 sds below its mean, with
# tilts from |rho| w of 1e-5 to 1e-2 on both sides of the small window limit.
exptilt_normal_tail_grid <- function() {
  tail_points <- expand.grid(
    z = c(-37, -27, -13, -7), pwindow = c(1, 7),
    scaled = c(1e-5, 1e-4, 1e-3, 1e-2), sign = c(-1, 1)
  )
  data.frame(
    d = -4 + 0.3 * tail_points$z, pwindow = tail_points$pwindow,
    rho = tail_points$sign * tail_points$scaled / tail_points$pwindow
  )
}

# Log reference for each row of the grid
exptilt_normal_tail_reference <- function(grid) {
  vapply(seq_len(nrow(grid)), function(i) {
    exptilt_normal_log_reference(
      grid$d[i], grid$pwindow[i], grid$rho[i], -4, 0.3
    )
  }, numeric(1))
}

# Log reference for one delay, or with `vectorised` the sum of the log PMF
# over the delays 0 to d, for an exponential, gamma or normal delay.
exptilt_log_reference <- function(dist_id, params, d, pwindow, rho,
                                  vectorised = FALSE) {
  if (vectorised) {
    cdf <- switch(as.character(dist_id),
      "2" = function(x) stats::pgamma(x, params[1], params[2]),
      "4" = function(x) stats::pexp(x, params[1]),
      "18" = function(x) stats::pnorm(x, params[1], params[2])
    )
    return(sum(log(diff(
      exptilt_reference(0:(d + 1), pwindow, rho, cdf)
    ))))
  }
  switch(as.character(dist_id),
    "2" = exptilt_gamma_log_reference(
      d, pwindow, rho, params[1], params[2]
    ),
    "4" = exptilt_gamma_log_reference(d, pwindow, rho, 1, params[1]),
    "18" = exptilt_normal_log_reference(
      d, pwindow, rho, params[1], params[2]
    )
  )
}

# Five point differences of the log reference in the parameters of the
# gradient model: the first delay parameter, the log of the second and the
# tilt.
exptilt_log_gradient <- function(dist_id, params, d, pwindow, rho,
                                 vectorised = FALSE) {
  theta <- c(params, rho)
  steps <- c(1e-5 * pmax(abs(params), 1), 1e-3 / pwindow)
  grad <- vapply(1:3, function(i) {
    at <- function(k) {
      shifted <- theta
      shifted[i] <- theta[i] + k * steps[i]
      exptilt_log_reference(
        dist_id, shifted[1:2], d, pwindow, shifted[3], vectorised
      )
    }
    (-at(2) + 8 * at(1) - 8 * at(-1) + at(-2)) / (12 * steps[i])
  }, numeric(1))
  grad[2] <- grad[2] * theta[2] + 1
  grad
}
