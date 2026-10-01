# Helpers for the Weibull and generalised gamma tilt tests. Both are in the
# Stacy family with CDF P(k, (x / scale)^shape), with k = 1 for the Weibull.

exptilt_stacy_families <- function() {
  stacy_family <- function(label, dist, args, shape, scale, k) {
    list(
      label = label, pdist = dist$p, rdist = dist$r, ddist = dist$d,
      args = args, rate = Inf, positive = TRUE,
      stacy = list(shape = shape, scale = scale, k = k)
    )
  }
  weibull <- list(p = pweibull, r = rweibull, d = dweibull)
  families <- list(
    stacy_family("weibull 1.5 5", weibull, list(shape = 1.5, scale = 5),
      1.5, 5, 1),
    stacy_family("weibull 0.7 5", weibull, list(shape = 0.7, scale = 5),
      0.7, 5, 1),
    stacy_family("weibull 3 2", weibull, list(shape = 3, scale = 2), 3, 2, 1)
  )
  if (requireNamespace("flexsurv", quietly = TRUE)) {
    gengamma <- list(
      p = flexsurv::pgengamma.orig, r = flexsurv::rgengamma.orig,
      d = flexsurv::dgengamma.orig
    )
    families <- c(families, list(
      stacy_family("gengamma 1.3 4 2.5", gengamma,
        list(shape = 1.3, scale = 4, k = 2.5), 1.3, 4, 2.5),
      stacy_family("gengamma 0.8 3 0.6", gengamma,
        list(shape = 0.8, scale = 3, k = 0.6), 0.8, 3, 0.6),
      stacy_family("gengamma 3 1 0.2", gengamma,
        list(shape = 3, scale = 1, k = 0.2), 3, 1, 0.2)
    ))
  }
  families
}

# Reference transform int_0^t exp(xi u) f(u) du
stacy_transform_reference <- function(family, t, xi) {
  vapply(t, function(tt) {
    stats::integrate(
      function(u) exp(xi * u) * do.call(family$ddist, c(list(u), family$args)),
      0, tt,
      rel.tol = 1e-12, abs.tol = 0, subdivisions = 2000L
    )$value
  }, numeric(1))
}

# Reference PMF from the delay survival function, which does not cancel in
# the upper tail where the CDF is within rounding of 1
exptilt_pmf_reference <- function(family, x, pwindow, rho) {
  survival <- function(u) {
    do.call(family$pdist, c(list(u), family$args, list(lower.tail = FALSE)))
  }
  vapply(x, function(xx) {
    stats::integrate(
      function(z) {
        (survival(xx - z) - survival(xx + 1 - z)) *
          exptilt_window_density(z, pwindow, rho)
      },
      0, pwindow,
      rel.tol = 1e-12, abs.tol = 0, subdivisions = 5000L,
      stop.on.error = FALSE
    )$value
  }, numeric(1))
}
