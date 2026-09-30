# Helpers for the truncated Gumbel primary event window tests.
#
# The reference values integrate the delay CDF against the window density
# with tight tolerances, unlike `pcens_cdf.default()` which uses the default
# `stats::integrate()` tolerances.

# Reference primary event censored CDF for any delay CDF `cdf(x)`. The
# integral is split at the kink where the delay CDF leaves zero.
gumbel_reference <- function(q, pwindow, mu, beta, cdf, positive = TRUE) {
  vapply(q, function(qq) {
    integrand <- function(z) {
      cdf(qq - z) * dtgumbel(z, 0, pwindow, mu, beta)
    }
    breaks <- c(0, if (positive && qq > 0 && qq < pwindow) qq, pwindow)
    sum(vapply(seq_len(length(breaks) - 1L), function(i) {
      stats::integrate(
        integrand, breaks[i], breaks[i + 1L],
        rel.tol = 1e-13, abs.tol = 0, subdivisions = 2000L
      )$value
    }, numeric(1)))
  }, numeric(1))
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
