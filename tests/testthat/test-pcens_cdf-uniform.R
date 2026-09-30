# The analytical uniform primary CDFs for the gamma, lognormal, Weibull and
# generalised gamma delays share `.pcens_cdf_uniform()`. It combines terms
# G(t) = t F(t) - E tilde F(t) at both ends of the primary window.
# These tests check the solutions against numerical integration at a tight
# tolerance, including tails, `q < pwindow`, `q` near 0 and edge parameters.

# Numerical reference: integrate the delay CDF over the primary window with
# a much tighter tolerance than `pcens_cdf.default()` uses. The upper limit
# is min(pwindow, q) because the delay CDF is 0 beyond it, which would put a
# kink inside the interval.
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

# Relative agreement with an absolute floor for values that are zero to
# working precision.
expect_close <- function(object, expected, rtol = 1e-6, atol = 1e-13,
                         info = NULL) {
  expect_true(
    all(abs(object - expected) <= atol + rtol * abs(expected)),
    info = info
  )
}

unif_cases <- function() {
  cases <- list(
    list(
      name = "gamma", pdist = pgamma,
      grid = list(
        list(shape = 0.5, scale = 2), list(shape = 1, scale = 1),
        list(shape = 3, scale = 2), list(shape = 20, scale = 0.5),
        list(shape = 0.1, scale = 10)
      )
    ),
    list(
      name = "lognormal", pdist = plnorm,
      grid = list(
        list(meanlog = 0, sdlog = 1), list(meanlog = 1.5, sdlog = 0.6),
        list(meanlog = 2, sdlog = 0.2), list(meanlog = -1, sdlog = 1.5),
        list(meanlog = 1, sdlog = 2)
      )
    ),
    list(
      name = "weibull", pdist = pweibull,
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

test_that("uniform primary analytical CDFs match tight numerical
  integration across parameters, windows and tails", {
  q <- c(1e-8, 1e-3, 0.05, 0.3, 0.75, 1, 1.5, 2.5, 6, 12, 40)
  for (case in unif_cases()) {
    for (args in case$grid) {
      for (pwindow in c(0.5, 1, 3)) {
        obj <- do.call(
          new_pcens, c(list(case$pdist, dunif), args)
        )
        info <- paste(
          case$name, toString(names(args)), toString(unlist(args)),
          "pwindow", pwindow
        )
        expect_close(
          pcens_cdf(obj, q, pwindow),
          do.call(unif_reference, c(list(case$pdist, q, pwindow), args)),
          info = info
        )
      }
    }
  }
})

test_that("uniform primary analytical CDFs are 0 at and below 0 and tend to
  1", {
  for (case in unif_cases()) {
    obj <- do.call(
      new_pcens, c(list(case$pdist, dunif), case$grid[[1]])
    )
    expect_identical(pcens_cdf(obj, c(-5, -0.1, 0), 1), c(0, 0, 0))
    expect_equal(pcens_cdf(obj, 1e4, 1), 1, tolerance = 1e-10)
  }
})

test_that("uniform primary analytical CDFs handle empty, missing and
  non-positive q", {
  for (case in unif_cases()) {
    obj <- do.call(
      new_pcens, c(list(case$pdist, dunif), case$grid[[1]])
    )
    expect_identical(pcens_cdf(obj, numeric(0), 1), numeric(0))
    expect_error(pcens_cdf(obj, c(1, NA), 1), "missing")
    expect_error(pcens_cdf(obj, NaN, 1), "missing")
    expect_identical(pcens_cdf(obj, -Inf, 1), 0)
  }
})

test_that("uniform primary analytical CDFs recycle a vector pwindow
  element-wise, including when some q are not positive", {
  q <- c(-1, 0.4, 2, 5, 0, 9)
  pwindow <- c(1, 2, 0.5, 3, 1, 4)
  for (case in unif_cases()) {
    obj <- do.call(
      new_pcens, c(list(case$pdist, dunif), case$grid[[3]])
    )
    expected <- vapply(
      seq_along(q), function(i) pcens_cdf(obj, q[[i]], pwindow[[i]]),
      numeric(1)
    )
    expect_equal(
      pcens_cdf(obj, q, pwindow), expected, tolerance = 1e-14,
      info = case$name
    )
  }
})

test_that(".pcens_cdf_uniform combines the terms at both ends of the window
  for a new delay family", {
  # An exponential delay is not a package analytical solution. Its CDF is
  # F(t) = 1 - exp(-rate t) and its partial expectation, the integral of
  # x f(x) from 0 to t, is (1 - exp(-rate t) (1 + rate t)) / rate. The
  # terms are G(t) = t F(t) - partial expectation.
  rate <- 0.4
  terms_fn <- function(t) {
    t * -expm1(-rate * t) - (1 - exp(-rate * t) * (1 + rate * t)) / rate
  }
  q <- c(-1, 0, 0.2, 0.9, 1, 3, 10, 40)
  for (pwindow in c(0.5, 1, 4)) {
    expect_close(
      .pcens_cdf_uniform(q, pwindow, terms_fn),
      unif_reference(pexp, q, pwindow, rate = rate)
    )
  }
})
