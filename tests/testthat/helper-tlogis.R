# Reference CDF by integrating the delay CDF against the window density to a
# relative 1e-12, split at the delay kink and around the location.
tlogis_reference <- function(q, pwindow, location, scale, cdf,
                             positive = TRUE) {
  vapply(q, function(qq) {
    integrand <- function(p) {
      cdf(qq - p) * dtlogis(p, 0, pwindow, location, scale)
    }
    centre <- min(max(location, 0), pwindow)
    spike <- centre + scale * c(-40, -20, -10, -5, -2, 0, 2, 5, 10, 20, 40)
    breaks <- sort(unique(c(
      0, pwindow,
      if (positive && qq > 0 && qq < pwindow) qq,
      spike[spike > 0 & spike < pwindow]
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

# Log CDF of a gamma delay with a truncated logistic primary, by integration
# with the tight tolerances of `tlogis_reference()`. With `max_delay` it is
# the sum of the log PMFs of the integer delays 0 to `max_delay` that the
# vectorised Stan PMF returns, with no truncation.
tlogis_gamma_lcdf_reference <- function(shape, rate, d, pwindow, location,
                                        scale, max_delay = NULL) {
  pdist <- function(x) stats::pgamma(x, shape, rate)
  if (is.null(max_delay)) {
    return(log(tlogis_reference(d, pwindow, location, scale, pdist)))
  }
  cdfs <- tlogis_reference(
    0:(max_delay + 1), pwindow, location, scale, pdist
  )
  sum(log(diff(cdfs)))
}

# Gradient of tlogis_gamma_lcdf_reference() in the shape by central
# differences in the log shape with two Richardson extrapolations.
tlogis_shape_grad_ref <- function(shape, rate, d, pwindow,
                                  location, scale,
                                  max_delay = NULL) {
  central <- function(h) {
    (tlogis_gamma_lcdf_reference(
      shape * (1 + h), rate, d, pwindow, location, scale, max_delay
    ) - tlogis_gamma_lcdf_reference(
      shape * (1 - h), rate, d, pwindow, location, scale, max_delay
    )) / (2 * shape * h)
  }
  h <- 1e-2
  r <- c(central(h), central(h / 2), central(h / 4))
  e <- (4 * r[2:3] - r[1:2]) / 3
  (16 * e[2] - e[1]) / 15
}
