# The analytical uniform primary CDFs for the gamma, lognormal, Weibull and
# generalised gamma delays share `.pcens_cdf_uniform()`.
# These tests check them against numerical integration at a tight tolerance.
# The agreement with Stan is in test-stan-primarycensored_analytical_cdf.R.

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
  # Exponential delay, with partial expectation (1 - exp(-rt) (1 + rt)) / r
  rate <- 0.4
  terms_fn <- function(t) {
    t * -expm1(-rate * t) - (1 - exp(-rate * t) * (1 + rate * t)) / rate
  }
  upper_fn <- function(t) -exp(-rate * t) / rate
  q <- c(-1, 0, 0.2, 0.9, 1, 3, 10, 40)
  for (pwindow in c(0.5, 1, 4)) {
    expect_close(
      .pcens_cdf_uniform(
        q, pwindow, terms_fn, upper_fn, 1 / rate,
        function(t) pexp(t, rate)
      ),
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

test_that("uniform primary analytical CDFs match numerical integration for
  windows many orders of magnitude narrower than the delay", {
  cases <- list(
    list(pgamma, list(shape = 500, scale = 100), 49967, 1e-6),
    list(pgamma, list(shape = 500, scale = 100), c(4.9e4, 5.1e4), 1e-5),
    list(pweibull, list(shape = 0.1, scale = 1), 4.29e6, 1e-6),
    list(pweibull, list(shape = 0.1, scale = 100), 4.29e8, 1e-3),
    list(pweibull, list(shape = 0.1, scale = 100), 4.29e8, 1e-6),
    list(pweibull, list(shape = 1.5, scale = 5), c(3, 8, 20), 1e-8),
    list(plnorm, list(meanlog = 2, sdlog = 2), c(5, 4e3, 1e5), 1e-7),
    list(plnorm, list(meanlog = 0, sdlog = 1), c(1, 3, 40), 1e-9)
  )
  for (case in cases) {
    obj <- do.call(new_pcens, c(list(case[[1]], dunif), case[[2]]))
    expect_close(
      pcens_cdf(obj, case[[3]], case[[4]]),
      do.call(unif_reference, c(list(case[[1]], case[[3]], case[[4]]),
                                case[[2]])),
      info = paste(toString(unlist(case[[2]])), "pwindow", case[[4]])
    )
  }
})

test_that("uniform primary analytical CDFs handle a vector pwindow with both
  narrow and wide windows", {
  q <- c(5, 4e3, 1e5, 3)
  pwindow <- c(1e-7, 1, 1e-2, 1e-9)
  for (case in unif_cases()) {
    obj <- do.call(
      new_pcens, c(list(case$pdist, dunif), case$grid[[3]])
    )
    expect_equal(
      pcens_cdf(obj, q, pwindow),
      vapply(
        seq_along(q), function(i) pcens_cdf(obj, q[[i]], pwindow[[i]]),
        numeric(1)
      ),
      tolerance = 1e-14,
      info = case$name
    )
  }
})

test_that(".check_pwindow errors for missing q and invalid pwindow", {
  expect_null(.check_pwindow(c(1, 2), c(0, 1)))
  expect_error(.check_pwindow(NA_real_, 1), "q must not")
  expect_error(.check_pwindow(1, numeric(0)), "pwindow")
  expect_error(.check_pwindow(1, NA_real_), "pwindow")
  expect_error(.check_pwindow(1, -1), "pwindow")
})

test_that("uniform primary analytical CDFs warn when q and a vector pwindow
  have lengths that do not recycle evenly", {
  for (case in unif_cases()) {
    obj <- do.call(
      new_pcens, c(list(case$pdist, dunif), case$grid[[3]])
    )
    expect_warning(
      pcens_cdf(obj, c(1, 2, 3), c(1, 2)), "multiple"
    )
    expect_warning(
      pcens_cdf(obj, c(1, 2), c(1, 2, 3)), "multiple"
    )
    expect_no_warning(pcens_cdf(obj, c(1, 2, 3, 4), c(1, 2)))
    expect_no_warning(pcens_cdf(obj, 1, c(1, 2, 3)))
    expect_no_warning(pcens_cdf(obj, c(1, 2, 3), 1))
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
    expect_equal(
      pcens_cdf(obj, 2, c(1, 2)),
      vapply(c(1, 2), function(w) pcens_cdf(obj, 2, w), numeric(1)),
      tolerance = 1e-14
    )
    expect_identical(pcens_cdf(obj, numeric(0), 1), numeric(0))
    expect_identical(pcens_cdf(obj, numeric(0), c(1, 2)), numeric(0))
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

test_that("uniform primary analytical CDFs match the empirical CDF of
  rprimarycensored samples", {
  set.seed(378)
  n <- 2e5
  for (case in unif_cases()) {
    for (pwindow in c(0.5, 3)) {
      args <- case$grid[[3]]
      obj <- do.call(new_pcens, c(list(case$pdist, dunif), args))
      draws <- do.call(
        rprimarycensored,
        c(list(n, case$rdist, pwindow = pwindow, swindow = 0), args)
      )
      q <- stats::quantile(draws, c(0.02, 0.1, 0.5, 0.9, 0.98))
      # Within five standard errors of the empirical CDF
      expect_close(
        unname(pcens_cdf(obj, q, pwindow)),
        vapply(q, function(x) mean(draws <= x), numeric(1)),
        rtol = 0, atol = 6e-3,
        info = paste(case$name, "pwindow", pwindow)
      )
    }
  }
})
