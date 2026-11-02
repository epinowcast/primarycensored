# Helpers for the lognormal tilt transform tests. The reference integrates
# the tilted density in z = (log u - meanlog) / sdlog over many pieces around
# the part within e^-60 of its peak.

# Lognormal delay parameter sets with a small, a typical and a large sdlog
# and a negative and a positive meanlog.
exptilt_lnorm_cases <- function() {
  list(
    list(meanlog = 1.6, sdlog = 0.5),
    list(meanlog = 0, sdlog = 1),
    list(meanlog = 2, sdlog = 0.25),
    list(meanlog = -1, sdlog = 1.5),
    list(meanlog = 1, sdlog = 1.8)
  )
}

# Delay CDF of a lognormal case as a function of the quantile alone.
exptilt_lnorm_cdf <- function(case) {
  function(x) stats::plnorm(x, case$meanlog, case$sdlog)
}

# Family entry in the format of `exptilt_families()`, for the helpers that
# build objects and labels.
exptilt_lnorm_family <- function(case) {
  list(
    label = sprintf("lognormal %g, %g", case$meanlog, case$sdlog),
    pdist = stats::plnorm,
    args = list(meanlog = case$meanlog, sdlog = case$sdlog),
    rate = Inf, positive = TRUE
  )
}

# Log of the transform of the lognormal density between `lower` and `upper`
# on the z scale, z = (log u - meanlog) / sdlog, for the tilt xi.
lnorm_tilt_integral <- function(lower, upper, meanlog, sdlog, xi) {
  log_integrand <- function(z) {
    xi * exp(meanlog + sdlog * z) + stats::dnorm(z, log = TRUE)
  }
  z_grid <- seq(lower, upper, length.out = 20001)
  values <- log_integrand(z_grid)
  values[is.na(values)] <- -Inf
  peak <- max(values)
  if (!is.finite(peak)) {
    return(-Inf)
  }
  keep <- which(values > peak - 60)
  lower <- z_grid[max(1L, min(keep) - 1L)]
  upper <- z_grid[min(length(z_grid), max(keep) + 1L)]
  if (upper <= lower) {
    return(peak)
  }
  pieces <- seq(lower, upper, length.out = 81)
  total <- 0
  for (i in seq_len(80)) {
    total <- total + stats::integrate(
      function(z) exp(log_integrand(z) - peak),
      pieces[i], pieces[i + 1L],
      rel.tol = 1e-13, abs.tol = 0, subdivisions = 200L,
      stop.on.error = FALSE
    )$value
  }
  peak + log(total)
}

# Reference log transform at each `t`. The upper transform needs xi <= 0.
lnorm_tilt_reference <- function(t, meanlog, sdlog, xi, upper = FALSE) {
  vapply(t, function(tt) {
    if (tt <= 0) {
      if (upper) {
        return(lnorm_tilt_integral(-45, 45, meanlog, sdlog, xi))
      }
      return(-Inf)
    }
    zt <- (log(tt) - meanlog) / sdlog
    if (upper) {
      lnorm_tilt_integral(zt, max(zt, 10) + 45, meanlog, sdlog, xi)
    } else {
      lnorm_tilt_integral(min(zt - 15, -45), zt, meanlog, sdlog, xi)
    }
  }, numeric(1))
}
