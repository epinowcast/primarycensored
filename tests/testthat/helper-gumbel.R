# Helpers for the truncated Gumbel primary event window tests

# Reference primary event censored CDF for any delay CDF `cdf(x)` with tight
# tolerances. The integral is taken in u = s(z) - s(w), with
# s(z) = exp(-(z - mu) / beta), where the window density is
# exp(-u) / (1 - exp(-Delta)), so a narrow spike is resolved. It needs
# (w - mu) / beta above -700 so that s(w) does not underflow.
gumbel_reference <- function(q, pwindow, mu, beta, cdf, positive = TRUE) {
  log_sw <- -(pwindow - mu) / beta
  delta <- exp(log_sw) * expm1(pwindow / beta)
  top <- min(delta, 60)
  mass <- -expm1(-delta)
  vapply(q, function(qq) {
    integrand <- function(u) {
      cdf((qq - pwindow) + beta * log1p(exp(log(u) - log_sw))) * exp(-u)
    }
    breaks <- c(0, top, top * c(1e-6, 1e-4, 1e-2, 0.1, 0.3, 0.6))
    if (positive && qq > 0 && qq < pwindow) {
      u_kink <- exp(log_sw) * expm1((pwindow - qq) / beta)
      if (u_kink > 0 && u_kink < top) breaks <- c(breaks, u_kink)
    }
    breaks <- sort(unique(breaks))
    sum(vapply(seq_len(length(breaks) - 1L), function(i) {
      stats::integrate(
        integrand, breaks[i], breaks[i + 1L],
        rel.tol = 5e-14, abs.tol = 0, subdivisions = 5000L,
        stop.on.error = FALSE
      )$value
    }, numeric(1))) / mass
  }, numeric(1))
}

# Largest error of a CDF, relative where the reference is at least `floor`
gumbel_error <- function(actual, expected, floor = 1e-12) {
  max(abs(actual - expected) / pmax(expected, floor))
}

# Delay families, with the observations q at which to compare the CDFs
gumbel_families <- function() {
  list(
    list(
      label = "exponential rate 1", pdist = pexp, args = list(rate = 1),
      q = c(0.005, 0.02, 0.05, 0.1, 0.5, 1, 2.5, 8), positive = TRUE
    ),
    list(
      label = "gamma shape 3 rate 1", pdist = pgamma,
      args = list(shape = 3, rate = 1),
      q = c(0.01, 0.03, 0.05, 0.1, 0.5, 1, 2.5, 8), positive = TRUE
    ),
    list(
      label = "gamma shape 0.6 rate 2", pdist = pgamma,
      args = list(shape = 0.6, rate = 2),
      q = c(0.002, 0.02, 0.1, 0.5, 1, 2.5), positive = TRUE
    ),
    list(
      label = "normal mean 3", pdist = pnorm,
      args = list(mean = 3, sd = 2),
      q = c(-5, -1, 0.3, 1, 2.5, 3, 5, 8, 12), positive = FALSE
    ),
    list(
      label = "normal mean -1", pdist = pnorm,
      args = list(mean = -1, sd = 3),
      q = c(-8, -1, 0.3, 1, 2.5, 4, 8), positive = FALSE
    )
  )
}

# The families with a series solution
gumbel_series_families <- function() gumbel_families()[4:5]

gumbel_cdf <- function(family) {
  function(x) do.call(family$pdist, c(list(x), family$args))
}

gumbel_object <- function(family, mu, beta) {
  do.call(
    new_pcens,
    c(
      list(
        pdist = family$pdist, dprimary = dtgumbel,
        primary_args = list(mu = mu, beta = beta)
      ),
      family$args
    )
  )
}

gumbel_label <- function(family, pwindow, mu, beta) {
  sprintf(
    "%s, pwindow = %g, mu = %g, beta = %g",
    family$label, pwindow, mu, beta
  )
}

# Windows where the density is a narrow spike, mu / beta of 15 to 50, and
# delays long relative to the window so that the series is not available
gumbel_spike_settings <- function() {
  list(
    c(mu = 1.5, beta = 0.1, w = 1),
    c(mu = 2, beta = 0.1, w = 1),
    c(mu = 1, beta = 0.02, w = 0.3),
    c(mu = 0.5, beta = 0.025, w = 1),
    c(mu = 1.5, beta = 0.05, w = 7),
    c(mu = 3, beta = 0.06, w = 2)
  )
}

gumbel_spike_families <- function() {
  list(
    list(
      label = "normal mean 3 sd 2", pdist = pnorm,
      args = list(mean = 3, sd = 2),
      q = c(-3, 0.5, 1, 3, 5, 8, 20), positive = FALSE
    ),
    list(
      label = "exponential rate 1", pdist = pexp, args = list(rate = 1),
      q = c(0.05, 0.5, 1, 3, 8, 20), positive = TRUE
    ),
    list(
      label = "gamma shape 3 rate 1", pdist = pgamma,
      args = list(shape = 3, rate = 1),
      q = c(0.5, 1, 3, 8, 20), positive = TRUE
    ),
    list(
      label = "lognormal meanlog 1 sdlog 0.5", pdist = plnorm,
      args = list(meanlog = 1, sdlog = 0.5),
      q = c(0.5, 1, 3, 8, 20), positive = TRUE
    ),
    list(
      label = "Weibull shape 2 scale 3", pdist = pweibull,
      args = list(shape = 2, scale = 3),
      q = c(0.5, 1, 3, 8, 20), positive = TRUE
    )
  )
}
