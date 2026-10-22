# The analytical uniform primary CDFs share `.pcens_cdf_uniform()`.
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
  if (requireNamespace("flexsurv", quietly = TRUE)) {
    # A small Q gives a large shape, where the terms at q < 0 are undefined
    obj <- new_pcens(
      flexsurv::pgengamma, dunif, mu = 1, sigma = 0.6, Q = 1e-3
    )
    expect_identical(pcens_cdf(obj, c(-2, -1e-12, 0, 0), 1), rep(0, 4))
  }
})

test_that("uniform primary analytical CDFs handle empty input and error for
  missing q and invalid pwindow", {
  for (case in unif_cases()) {
    obj <- do.call(
      new_pcens, c(list(case$pdist, dunif), case$grid[[1]])
    )
    expect_identical(pcens_cdf(obj, numeric(0), 1), numeric(0))
    expect_identical(pcens_cdf(obj, numeric(0), c(1, 2)), numeric(0))
    expect_error(pcens_cdf(obj, c(1, NA), 1), "missing")
    expect_error(pcens_cdf(obj, NaN, 1), "missing")
    expect_error(pcens_cdf(obj, c(1, 2), NA_real_), "pwindow")
    expect_error(pcens_cdf(obj, c(1, 2), c(1, NA)), "pwindow")
    expect_error(pcens_cdf(obj, c(1, 2), -1), "pwindow")
  }
})

test_that("uniform primary analytical CDFs with a vector pwindow match
  element-wise calls across non-positive, narrow and upper tail cases", {
  q <- c(-1, 0.4, 2, 5, 0, 9, 5, 4e3, 1e5, 3, 1e6, 1e9, 60, 25, Inf)
  pwindow <- c(
    1, 2, 0.5, 3, 1, 4, 1e-7, 1, 1e-2, 1e-9, 0.1, 5, 0.01, 1, 1
  )
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
    expect_equal(
      pcens_cdf(obj, 2, c(1, 2)),
      vapply(c(1, 2), function(w) pcens_cdf(obj, 2, w), numeric(1)),
      tolerance = 1e-14
    )
  }
})

test_that(".pcens_cdf_uniform combines the terms at both ends of the window
  for a new delay family", {
  # Exponential delay, with partial expectation (1 - exp(-rt) (1 + rt)) / r
  rate <- 0.4
  terms_fn <- function(t, p) {
    t * -expm1(-p$rate * t) - p$mean * (1 - exp(-p$rate * t) * (1 + p$rate * t))
  }
  upper_fn <- function(t, p) -p$mean * exp(-p$rate * t)
  delay_cdf <- function(t, p) pexp(t, p$rate)
  q <- c(-1, 0, 0.2, 0.9, 1, 3, 10, 40)
  for (pwindow in c(0.5, 1, 4)) {
    expect_close(
      .pcens_cdf_uniform(
        q, pwindow, terms_fn, upper_fn, delay_cdf,
        list(rate = rate, mean = 1 / rate)
      ),
      unif_reference(pexp, q, pwindow, rate = rate)
    )
  }
})

test_that(".uniform_endpoint_diff evaluates each distinct window end once
  when enough lower ends match a delay", {
  counter <- new.env()
  fn <- function(t, p) {
    counter$n <- counter$n + length(t)
    t^2
  }
  count_evaluations <- function(q, pwindow, share = TRUE) {
    counter$n <- 0L
    lo <- pmax.int(q - pwindow, 0)
    expect_identical(
      .uniform_endpoint_diff(fn, q, lo, list(), share), q^2 - lo^2
    )
    counter$n
  }
  distinct_ends <- function(q, pwindow) {
    length(unique(c(q, pmax.int(q - pwindow, 0))))
  }
  q <- 0:80
  expect_identical(count_evaluations(q, 2), distinct_ends(q, 2))
  q_half <- seq(0, 40, by = 0.5)
  expect_identical(count_evaluations(q_half, 1.5), distinct_ends(q_half, 1.5))
  q_mixed <- c(seq(0.37, 9.9, length.out = 30), 0:50)
  expect_identical(count_evaluations(q_mixed, 1), distinct_ends(q_mixed, 1))
  # Not shared without the flag, for short inputs, or when few lower ends
  # match a delay
  expect_identical(count_evaluations(q, 2, share = FALSE), 2L * length(q))
  expect_identical(count_evaluations(0:30, 1), 62L)
  expect_identical(count_evaluations(q, 1.5), 2L * length(q))
  q_none <- seq(1.37, 40, length.out = 60)
  expect_identical(count_evaluations(q_none, 1), 2L * length(q_none))
  q_zeros <- sort(c(0, 0, 0, seq(1.37, 40, length.out = 200)))
  expect_identical(count_evaluations(q_zeros, 1), 2L * length(q_zeros))
})

test_that("uniform primary analytical CDFs that share terms match single
  point evaluation", {
  grids <- list(
    list(q = 0:100, pwindow = 1),
    list(q = 0:100, pwindow = 3),
    list(q = seq(0, 40, by = 0.5), pwindow = 1.5),
    list(q = seq(0, 20, by = 0.25), pwindow = 0.25),
    list(q = c(seq(0.37, 9.9, length.out = 30), 0:50), pwindow = 1),
    list(q = c(0:100, 2e3, 5e3), pwindow = 1)
  )
  for (case in unif_cases()) {
    for (args in case$grid[1:2]) {
      obj <- do.call(new_pcens, c(list(case$pdist, dunif), args))
      for (grid in grids) {
        expect_equal(
          pcens_cdf(obj, grid$q, grid$pwindow),
          vapply(
            grid$q, function(x) pcens_cdf(obj, x, grid$pwindow), numeric(1)
          ),
          tolerance = 1e-14,
          info = paste(case$name, toString(unlist(args)))
        )
      }
    }
  }
})

test_that("pprimarycensored and dprimarycensored agree on a lattice where
  uniform primary terms are shared", {
  x <- 0:80
  for (case in unif_cases()) {
    args <- c(list(pdist = case$pdist, pwindow = 2), case$grid[[3]])
    p <- do.call(pprimarycensored, c(list(q = x), args))
    expect_equal(
      p,
      vapply(
        x, function(i) do.call(pprimarycensored, c(list(q = i), args)),
        numeric(1)
      ),
      tolerance = 1e-14,
      info = case$name
    )
    d <- do.call(dprimarycensored, c(list(x = x, D = 81), args))
    p_lo <- do.call(pprimarycensored, c(list(q = x, L = 0, D = 81), args))
    p_up <- do.call(
      pprimarycensored, c(list(q = x + 1, L = 0, D = 81), args)
    )
    expect_equal(d, p_up - p_lo, tolerance = 1e-12, info = case$name)
  }
})

test_that("uniform primary analytical CDFs evaluate vector delay parameters
  per point", {
  q <- c(0:10, 2e3, 1e7)
  shape <- seq(1, 3, length.out = 13)
  obj <- new_pcens(pgamma, dunif, shape = shape, scale = 1)
  expected <- vapply(seq_along(q), function(i) {
    pcens_cdf(new_pcens(pgamma, dunif, shape = shape[i], scale = 1), q[i], 1)
  }, numeric(1))
  expect_no_warning(pcens_cdf(obj, q, 1))
  expect_equal(pcens_cdf(obj, q, 1), expected, tolerance = 1e-14)
  expect_equal(
    pcens_cdf(obj, 2, 1),
    vapply(shape, function(s) {
      pcens_cdf(new_pcens(pgamma, dunif, shape = s, scale = 1), 2, 1)
    }, numeric(1)),
    tolerance = 1e-14
  )
  expect_warning(pcens_cdf(obj, 1:5, 1), "multiple")
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

test_that("uniform primary analytical CDFs match numerical integration where
  the upper tail form is used and the tail is not negligible", {
  q <- c(30, 60, 120, 240)
  for (case in unif_cases()) {
    for (args in case$grid) {
      obj <- do.call(new_pcens, c(list(case$pdist, dunif), args))
      expect_close(
        pcens_cdf(obj, q, 0.01),
        do.call(unif_reference, c(list(case$pdist, q, 0.01), args)),
        info = paste(case$name, toString(unlist(args)))
      )
    }
  }
})

test_that("uniform primary analytical CDFs agree on either side of the
  switch to the narrow window form", {
  d <- 5
  for (case in unif_cases()) {
    for (args in case$grid[1:3]) {
      obj <- do.call(new_pcens, c(list(case$pdist, dunif), args))
      for (ratio in c(0.999, 1.001)) {
        w <- d / (ratio * .narrow_window_ratio)
        narrow <- pcens_cdf(obj, d, w)
        wide <- with_mocked_bindings(
          pcens_cdf(obj, d, w),
          .narrow_window_ratio = Inf
        )
        expect_close(
          narrow, wide, rtol = 1e-8,
          info = paste(case$name, toString(unlist(args)), "ratio", ratio)
        )
      }
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
