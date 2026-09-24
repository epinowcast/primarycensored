# The fast paths added for #347 must give the same results as the paths
# they skip. Each test switches a fast path off with a mock, or compares
# against per-observation calls to the exported functions.

equiv_delays <- list(
  list(pdist = pgamma, args = list(shape = 3, scale = 2)),
  list(pdist = plnorm, args = list(meanlog = 1, sdlog = 0.5)),
  list(pdist = pweibull, args = list(shape = 1.5, scale = 3))
)

equiv_primaries <- list(
  list(dprimary = dunif, args = list()),
  list(dprimary = dexpgrowth, args = list(r = 0.2))
)

# Call an exported function for a delay and primary
call_pcens <- function(fn, x, delay, primary, ...) {
  suppressMessages(do.call(
    fn,
    c(
      list(
        x, delay$pdist, ...,
        dprimary = primary$dprimary, primary_args = primary$args
      ),
      delay$args
    )
  ))
}

# Evaluate `code` with names found only by .extract_function_name()
with_slow_names <- function(code) {
  testthat::local_mocked_bindings(.registry_name = function(func) NULL)
  code
}

test_that("registry name lookups match .extract_function_name()", {
  for (delay in equiv_delays) {
    for (primary in equiv_primaries) {
      for (bounds in list(c(-Inf, Inf), c(-Inf, 15), c(1, 15), c(1, Inf))) {
        for (pwindow in c(1, 2)) {
          x <- seq(max(0, bounds[1]), min(14.5, bounds[2] - 0.5), by = 0.5)
          args <- list(pwindow = pwindow, L = bounds[1], D = bounds[2])
          for (fn in list(dprimarycensored, pprimarycensored)) {
            fast <- do.call(call_pcens, c(list(fn, x, delay, primary), args))
            slow <- with_slow_names(
              do.call(call_pcens, c(list(fn, x, delay, primary), args))
            )
            expect_identical(fast, slow)
          }
        }
      }
      new_args <- c(
        list(delay$pdist, primary$dprimary, primary_args = primary$args),
        delay$args
      )
      expect_identical(
        do.call(new_pcens, new_args),
        with_slow_names(do.call(new_pcens, new_args))
      )
    }
  }
})

test_that("pcens_pmf gives the same results for sorted and unsorted x", {
  obj <- new_pcens(pgamma, dunif, shape = 3, scale = 2)
  x <- seq(0, 12, by = 0.5)
  for (swindow in c(0.3, 1, 2.5)) {
    for (D in c(Inf, 13)) {
      sorted <- suppressMessages(pcens_pmf(obj, x, 1, swindow, D = D))
      unsorted <- suppressMessages(pcens_pmf(obj, rev(x), 1, swindow, D = D))
      expect_identical(rev(unsorted), sorted)
    }
  }
})

equiv_fit_params <- data.frame(
  swindow = c(1, 1, 2, 2, 1, 1),
  pwindow = c(1, 1, 1, 2, 2, 2),
  L = c(-Inf, -Inf, -Inf, 0, 0, 0),
  D = c(Inf, Inf, 15, 15, 20, 20)
)

test_that(".dpcens and .ppcens match per-observation calls", {
  x <- c(1, 4, 2, 6, 3, 8)
  p <- equiv_fit_params
  for (delay in equiv_delays) {
    for (primary in equiv_primaries) {
      expected_d <- vapply(seq_along(x), function(i) {
        call_pcens(
          dprimarycensored, x[i], delay, primary,
          pwindow = p$pwindow[i], swindow = p$swindow[i], L = p$L[i],
          D = p$D[i]
        )
      }, numeric(1))
      expected_p <- vapply(seq_along(x), function(i) {
        call_pcens(
          pprimarycensored, x[i], delay, primary,
          pwindow = p$pwindow[i], L = p$L[i], D = p$D[i]
        )
      }, numeric(1))
      fit_args <- c(
        list(x, p, delay$pdist, primary$dprimary, primary$args), delay$args
      )
      expect_identical(do.call(.dpcens, fit_args), expected_d)
      expect_identical(do.call(.ppcens, fit_args), expected_p)
    }
  }
})

test_that(".dpcens and .ppcens give NaN when any input errors", {
  # fitdistrplus probes the fitted functions with inputs like these
  for (x in list(c(0, 1, NA), c(0, 1, Inf, NaN, -1))) {
    for (fn in list(.dpcens, .ppcens)) {
      expect_identical(
        fn(x, equiv_fit_params, pgamma, dunif, list(), shape = 2, scale = 1),
        rep(NaN, length(x))
      )
    }
  }
  expect_length(
    .ppcens(
      numeric(0), equiv_fit_params, pgamma, dunif, list(),
      shape = 2, scale = 1
    ),
    0L
  )
})

test_that(".dpcens and .ppcens treat NULL primary_args as an empty list", {
  x <- c(1, 4, 2, 6, 3, 8)
  for (fn in list(.dpcens, .ppcens)) {
    expect_identical(
      fn(x, equiv_fit_params, pgamma, dunif, NULL, shape = 3, scale = 2),
      fn(x, equiv_fit_params, pgamma, dunif, list(), shape = 3, scale = 2)
    )
  }
})

test_that("a custom primary without a CDF matches a new_pcens object", {
  dprim <- function(x, min, max) dunif(x, min, max)
  obj <- new_pcens(pgamma, dprim, shape = 3, scale = 2)
  expect_null(obj$pprimary)
  x <- 0:10
  expect_identical(
    dprimarycensored(x, pgamma, dprimary = dprim, shape = 3, scale = 2),
    pcens_pmf(obj, x, 1)
  )
  expect_identical(
    pprimarycensored(x, pgamma, dprimary = dprim, shape = 3, scale = 2),
    pcens_cdf(obj, x, 1)
  )
})

test_that("fitdistdoublecens gives the same fit when rebuilding per call", {
  skip_if_not_installed("fitdistrplus")
  withr::local_seed(101)
  n <- 100
  pwindows <- rep_len(c(1, 2), n)
  Ds <- rep(c(15, Inf), each = n / 2)
  delays <- vapply(seq_len(n), function(i) {
    rprimarycensored(
      1, rgamma,
      shape = 3, scale = 2, pwindow = pwindows[i], swindow = 1, D = Ds[i]
    )
  }, numeric(1))
  censdata <- data.frame(
    left = delays, right = delays + 1, pwindow = pwindows, D = Ds
  )
  fit <- function() {
    fitdistdoublecens(
      censdata, "gamma",
      start = list(shape = 1, scale = 1), truncation_check_multiplier = NULL
    )
  }
  cached <- fit()

  # Drop the cache so every evaluation builds a new object and groups
  fit_state <- .fit_pcens_state
  testthat::local_mocked_bindings(
    .fit_pcens_state = function(cache, ...) fit_state(NULL, ...)
  )
  rebuilt <- fit()

  expect_identical(cached$estimate, rebuilt$estimate)
  expect_identical(cached$loglik, rebuilt$loglik)
  expect_identical(cached$vcov, rebuilt$vcov)
})
