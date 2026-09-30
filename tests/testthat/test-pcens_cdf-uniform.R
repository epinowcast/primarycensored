# The analytical uniform primary CDFs for the gamma, lognormal, Weibull and
# generalised gamma delays share `.pcens_cdf_uniform()`. It combines terms
# G(t) = t F(t) - E tilde F(t) at both ends of the primary window.
# These tests check the solutions against numerical integration at a tight
# tolerance, including tails, `q < pwindow`, `q` near 0 and edge parameters.

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
      pcens_cdf(obj, q, pwindow), expected,
      tolerance = 1e-14,
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
