# Helpers for the truncated Gumbel primary event window tests.
#
# The reference values integrate the delay CDF against the window density
# with tight tolerances, unlike `pcens_cdf.default()` which uses the default
# `stats::integrate()` tolerances.

# Reference primary event censored CDF for any delay CDF `cdf(x)`.
#
# The integral is taken in u = s(z) - s(w), with s(z) = exp(-(z - mu) / beta),
# where the window density is exp(-u) / (1 - exp(-Delta)) on [0, Delta] and
# z = mu - beta log(s(w) + u). The weight is smooth, so a narrow spike of
# the density in z, for a large mu / beta, is resolved. The package
# integrates in u only for mu at or above the window end and in z otherwise,
# so the two share the delay argument only. The integral is truncated at
# u = 60, which leaves a mass of 1e-26, and split at the kink where the
# delay CDF leaves zero and on a ladder of scales. It needs
# (w - mu) / beta above -700 so that s(w) does not underflow.
gumbel_reference <- function(q, pwindow, mu, beta, cdf, positive = TRUE) {
  log_sw <- -(pwindow - mu) / beta
  delta <- exp(log_sw) * expm1(pwindow / beta)
  top <- min(delta, 60)
  mass <- -expm1(-delta)
  vapply(q, function(qq) {
    # The delay at q - z = (q - w) + beta log(1 + u / s(w)), which keeps a
    # small difference when u / s(w) is tiny
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
# and absolute, relative to `floor`, below it.
gumbel_error <- function(actual, expected, floor = 1e-12) {
  max(abs(actual - expected) / pmax(expected, floor))
}

# Delay families with a transform for positive tilts. The exponential and
# gamma rates are large as the series needs the rate above N / beta, so
# these are short delays in the units of the window.
gumbel_families <- function() {
  list(
    list(
      label = "exponential rate 60", pdist = pexp, args = list(rate = 60),
      q = c(0.005, 0.02, 0.05, 0.1, 0.5, 1, 2.5), positive = TRUE
    ),
    list(
      label = "gamma shape 3 rate 80", pdist = pgamma,
      args = list(shape = 3, rate = 80),
      q = c(0.01, 0.03, 0.05, 0.1, 0.5, 1, 2.5), positive = TRUE
    ),
    list(
      label = "gamma shape 0.6 rate 100", pdist = pgamma,
      args = list(shape = 0.6, rate = 100),
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

# Label for a failed expectation.
gumbel_label <- function(family, pwindow, mu, beta) {
  sprintf(
    "%s, pwindow = %g, mu = %g, beta = %g",
    family$label, pwindow, mu, beta
  )
}

# Windows where the density is a narrow spike, mu / beta of 15 to 50, at the
# upper edge of the window, inside it, and at a very narrow window. The
# delays are long relative to the window so that the series is not
# available for the exponential and the gamma.
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
    )
  )
}
