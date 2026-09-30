# Reference for the uniform primary analytical CDFs, integrating the delay CDF
# over the primary window at a tighter tolerance than `pcens_cdf.default()`.
# The upper limit is min(pwindow, q) as the delay CDF is 0 beyond it.
unif_reference <- function(pdist, q, pwindow, ...) {
  vapply(
    q,
    function(d) {
      if (d <= 0) {
        return(0)
      }
      upper <- min(pwindow, d)
      stats::integrate(
        function(p) pdist(d - p, ...) / pwindow,
        lower = 0, upper = upper, rel.tol = 1e-12, abs.tol = 0,
        subdivisions = 1000L
      )$value
    },
    numeric(1)
  )
}

# Relative agreement with an absolute floor
expect_close <- function(object, expected, rtol = 1e-6, atol = 1e-13,
                         info = NULL) {
  testthat::expect_true(
    all(abs(object - expected) <= atol + rtol * abs(expected)),
    info = info
  )
}

unif_cases <- function() {
  cases <- list(
    list(
      name = "gamma", pdist = pgamma, rdist = rgamma,
      grid = list(
        list(shape = 0.5, scale = 2), list(shape = 1, scale = 1),
        list(shape = 3, scale = 2), list(shape = 20, scale = 0.5),
        list(shape = 0.1, scale = 10)
      )
    ),
    list(
      name = "lognormal", pdist = plnorm, rdist = rlnorm,
      grid = list(
        list(meanlog = 0, sdlog = 1), list(meanlog = 1.5, sdlog = 0.6),
        list(meanlog = 2, sdlog = 0.2), list(meanlog = -1, sdlog = 1.5),
        list(meanlog = 1, sdlog = 2)
      )
    ),
    list(
      name = "weibull", pdist = pweibull, rdist = rweibull,
      grid = list(
        list(shape = 0.7, scale = 3), list(shape = 1, scale = 2),
        list(shape = 1.6, scale = 6), list(shape = 4, scale = 5),
        list(shape = 0.3, scale = 1)
      )
    )
  )
  if (requireNamespace("flexsurv", quietly = TRUE)) {
    cases[[4]] <- list(
      name = "gengamma", pdist = flexsurv::pgengamma.orig,
      rdist = flexsurv::rgengamma.orig,
      grid = list(
        list(shape = 1.5, scale = 4, k = 1.2),
        list(shape = 0.8, scale = 2, k = 0.5),
        list(shape = 2.5, scale = 6, k = 3),
        list(shape = 1, scale = 2, k = 2)
      )
    )
  }
  cases
}
