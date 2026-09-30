# Helpers for the truncated logistic primary event window tests.
#
# The reference values integrate the delay CDF against the window density
# with tight tolerances. Unlike `pcens_cdf.default()`, which uses the default
# `stats::integrate()` tolerances, it is accurate to about 1e-12 relative,
# including deep in the lower tail. The integral is split at the kink where
# the CDF of a delay on the non-negative reals leaves zero, and at the
# location of the window.
tlogis_reference <- function(q, pwindow, location, scale, cdf,
                             positive = TRUE) {
  vapply(q, function(qq) {
    integrand <- function(p) {
      cdf(qq - p) * dtlogis(p, 0, pwindow, location, scale)
    }
    breaks <- sort(unique(c(
      0, pwindow,
      if (positive && qq > 0 && qq < pwindow) qq,
      if (location > 0 && location < pwindow) location
    )))
    sum(vapply(seq_len(length(breaks) - 1L), function(i) {
      stats::integrate(
        integrand, breaks[i], breaks[i + 1L],
        rel.tol = 1e-12, abs.tol = 0, subdivisions = 2000L
      )$value
    }, numeric(1)))
  }, numeric(1))
}

# Delay CDF of a family from `exptilt_families()` as a function of the
# quantile alone.
tlogis_object <- function(family, location, scale) {
  do.call(
    new_pcens,
    c(
      list(
        pdist = family$pdist, dprimary = dtlogis,
        primary_args = list(location = location, scale = scale)
      ),
      family$args
    )
  )
}

# Label for a failed expectation.
tlogis_label <- function(family, pwindow, location, scale) {
  sprintf(
    "%s, pwindow = %g, location = %g, scale = %g",
    family$label, pwindow, location, scale
  )
}
