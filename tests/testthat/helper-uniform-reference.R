# Reference for the uniform primary analytical CDFs, integrating the delay CDF
# over the primary window at a tighter tolerance than `pcens_cdf.default()`.
# The upper limit is min(pwindow, q - lower) as the delay CDF is 0 beyond it.
unif_reference <- function(pdist, q, pwindow, ..., lower = 0) {
  vapply(
    q,
    function(d) {
      if (d <= lower) {
        return(0)
      }
      upper <- min(pwindow, d - lower)
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

# Lower end of the delay support of a case
case_lower <- function(case) if (is.null(case$lower)) 0 else case$lower

# The delays and parameters on which every uniform primary analytical CDF is
# checked
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
    ),
    list(
      name = "exponential", pdist = pexp, rdist = rexp,
      grid = list(
        list(rate = 0.4), list(rate = 1), list(rate = 0.1), list(rate = 3)
      )
    ),
    list(
      name = "beta", pdist = pbeta, rdist = rbeta,
      grid = list(
        list(shape1 = 2, shape2 = 3), list(shape1 = 0.5, shape2 = 0.5),
        list(shape1 = 1, shape2 = 4), list(shape1 = 5, shape2 = 1.5)
      )
    ),
    list(
      name = "chi-square", pdist = pchisq, rdist = rchisq,
      grid = list(list(df = 4), list(df = 1), list(df = 0.5), list(df = 10))
    ),
    list(
      name = "normal", pdist = pnorm, rdist = rnorm, lower = -Inf,
      grid = list(
        list(mean = 0, sd = 1), list(mean = 3, sd = 2),
        list(mean = 10, sd = 0.3), list(mean = -2, sd = 3)
      )
    )
  )
  if (requireNamespace("flexsurv", quietly = TRUE)) {
    cases[[length(cases) + 1L]] <- list(
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
  if (requireNamespace("extraDistr", quietly = TRUE)) {
    cases[[length(cases) + 1L]] <- list(
      name = "pareto", pdist = extraDistr::ppareto,
      rdist = extraDistr::rpareto,
      grid = list(
        list(a = 2.5, b = 0.5), list(a = 1, b = 1), list(a = 0.7, b = 0.2),
        list(a = 5, b = 2)
      )
    )
    cases[[length(cases) + 1L]] <- list(
      name = "inverse gamma", pdist = extraDistr::pinvgamma,
      rdist = extraDistr::rinvgamma,
      grid = list(
        list(alpha = 3, beta = 2), list(alpha = 1.5, beta = 1),
        list(alpha = 2, beta = 3), list(alpha = 10, beta = 20),
        list(alpha = 1.1, beta = 0.5)
      )
    )
    cases[[length(cases) + 1L]] <- list(
      name = "inverse chi-square", pdist = extraDistr::pinvchisq,
      rdist = extraDistr::rinvchisq,
      grid = list(
        list(nu = 5), list(nu = 3), list(nu = 4, tau = 2),
        list(nu = 10, tau = 2), list(nu = 2.5, tau = 0.3)
      )
    )
  }
  cases
}
