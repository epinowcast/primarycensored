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

test_that("uniform primary analytical CDFs are 1 at Inf and 0 at -Inf", {
  for (case in unif_cases()) {
    for (args in case$grid) {
      obj <- do.call(new_pcens, c(list(case$pdist, dunif), args))
      expect_identical(
        pcens_cdf(obj, c(-Inf, Inf), 1), c(0, 1), info = case$name
      )
      # Infinite delays sit alongside finite ones and a vector pwindow
      expect_equal(
        pcens_cdf(obj, c(Inf, 2, -Inf, Inf), c(1, 1, 2, 0.5))[c(1, 3, 4)],
        c(1, 0, 1),
        info = case$name
      )
    }
  }
})

test_that("uniform primary analytical CDFs are 1 far into the upper tail,
  where the lower tail form cancels", {
  q <- c(1e6, 1e9, 1e12, 1e14, 1e15, 1e16, 1e17)
  for (case in unif_cases()) {
    for (pwindow in c(0.1, 1, 5)) {
      obj <- do.call(
        new_pcens, c(list(case$pdist, dunif), case$grid[[1]])
      )
      expect_equal(
        pcens_cdf(obj, q, pwindow), rep(1, length(q)),
        tolerance = 1e-9,
        info = paste(case$name, "pwindow", pwindow)
      )
    }
  }
})

test_that("uniform primary analytical CDFs match numerical integration
  across the switch to the upper tail form and with narrow windows", {
  # Windows start above and below the delay mean for each case.
  q <- c(0.5, 2, 5, 8, 15, 30, 60, 150)
  for (case in unif_cases()) {
    for (args in case$grid[1:3]) {
      for (pwindow in c(1e-3, 0.05, 2)) {
        obj <- do.call(new_pcens, c(list(case$pdist, dunif), args))
        expect_close(
          pcens_cdf(obj, q, pwindow),
          do.call(unif_reference, c(list(case$pdist, q, pwindow), args)),
          info = paste(
            case$name, toString(unlist(args)), "pwindow", pwindow
          )
        )
      }
    }
  }
})

test_that("uniform primary analytical CDFs check pwindow", {
  for (case in unif_cases()) {
    obj <- do.call(
      new_pcens, c(list(case$pdist, dunif), case$grid[[3]])
    )
    expect_error(pcens_cdf(obj, c(1, 2), NA_real_), "pwindow")
    expect_error(pcens_cdf(obj, c(1, 2), c(1, NA)), "pwindow")
    expect_error(pcens_cdf(obj, c(1, 2), -1), "pwindow")
    # A vector pwindow is recycled against q, as a longer pwindow was before
    expect_equal(
      pcens_cdf(obj, 2, c(1, 2)),
      vapply(c(1, 2), function(w) pcens_cdf(obj, 2, w), numeric(1)),
      tolerance = 1e-14
    )
    expect_identical(pcens_cdf(obj, numeric(0), 1), numeric(0))
  }
})

test_that("uniform primary analytical CDFs give the delay CDF for the zero
  width elements of a vector pwindow", {
  q <- c(2, 2, 5, -1, 0)
  pwindow <- c(0, 1, 0, 0, 0)
  for (case in unif_cases()) {
    args <- case$grid[[3]]
    obj <- do.call(new_pcens, c(list(case$pdist, dunif), args))
    expected <- c(
      do.call(case$pdist, c(list(2), args)),
      pcens_cdf(obj, 2, 1),
      do.call(case$pdist, c(list(5), args)),
      0, 0
    )
    expect_equal(
      pcens_cdf(obj, q, pwindow), expected, tolerance = 1e-14,
      info = case$name
    )
  }
})

test_that("uniform primary analytical CDFs with a vector pwindow match
  element-wise calls when some windows use the upper tail form", {
  q <- c(0.5, 1e6, 3, 1e9, 60, 25, Inf, -2)
  pwindow <- c(1, 0.1, 2, 5, 0.01, 1, 1, 1)
  for (case in unif_cases()) {
    obj <- do.call(
      new_pcens, c(list(case$pdist, dunif), case$grid[[1]])
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
    expect_equal(expected[c(2, 4)], c(1, 1), tolerance = 1e-9)
  }
})
