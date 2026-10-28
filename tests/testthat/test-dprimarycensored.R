test_that("dprimarycensored sums to 1 for discrete values", {
  pwindow <- 1
  D <- 10
  pmf <- dpcens(
    0:(D - 1), plnorm, pwindow,
    D = D, meanlog = 1, sdlog = 1
  )
  expect_equal(sum(pmf), 1, tolerance = 1e-6)
})

test_that("dprimarycensored handles log probabilities", {
  pwindow <- 1
  D <- 10
  pmf <- dpcens(
    0:(D - 1), plnorm, pwindow,
    D = D, meanlog = 1, sdlog = 1
  )
  log_pmf <- dpcens(
    0:(D - 1), plnorm, pwindow,
    D = D, meanlog = 1, sdlog = 1, log = TRUE
  )
  expect_equal(exp(log_pmf), pmf, tolerance = 1e-6)
})

test_that("dprimarycensored handles non-finite D", {
  pwindow <- 1
  D <- Inf
  pmf <- dpcens(
    0:100, plnorm, pwindow,
    D = D, meanlog = 1, sdlog = 1
  )
  expect_lt(sum(pmf), 1)
  expect_equal(sum(pmf), 1, tolerance = 0.01)
})

test_that("dprimarycensored matches difference of pprimarycensored", {
  x <- c(1, 2, 3)
  pwindow <- 5
  swindow <- 0.5
  D <- 10

  pmf <- dpcens(
    x, plnorm,
    pwindow = pwindow, swindow = swindow, D = D,
    meanlog = 0, sdlog = 1
  )
  cdf_diff <- sapply(x, function(xi) {
    ppcens(
      xi + swindow, plnorm,
      pwindow = pwindow, D = D,
      meanlog = 0, sdlog = 1
    ) -
      ppcens(
        xi, plnorm,
        pwindow = pwindow, D = D,
        meanlog = 0, sdlog = 1
      )
  })

  expect_equal(pmf, cdf_diff, tolerance = 1e-6)
})

test_that(
  "dprimarycensored clips upper endpoint at D when right > D and left < D",
  {
    pwindow <- 1
    swindow <- 1
    D <- 10
    # When the lower edge is strictly below D, the clipped likelihood
    # equals (F_cens(D) - F_cens(D - 0.5)) / F_cens(D), using the
    # un-normalised censored CDF (computed via D = Inf).
    res <- dpcens(
      D - 0.5, plnorm,
      pwindow = pwindow, swindow = swindow, D = D,
      meanlog = 0, sdlog = 1
    )
    raw_F_D <- pprimarycensored(
      D, plnorm,
      pwindow = pwindow, L = -Inf, D = Inf, meanlog = 0, sdlog = 1
    )
    raw_F_lower <- pprimarycensored(
      D - 0.5, plnorm,
      pwindow = pwindow, L = -Inf, D = Inf, meanlog = 0, sdlog = 1
    )
    expected <- (raw_F_D - raw_F_lower) / raw_F_D
    expect_equal(res, expected, tolerance = 1e-8)
  }
)

test_that(
  "dprimarycensored errors when any x >= D",
  {
    pwindow <- 1
    swindow <- 1
    D <- 10
    # x == D: degenerate observation; rejected.
    expect_error(
      dpcens(
        D, plnorm,
        pwindow = pwindow, swindow = swindow, D = D,
        meanlog = 0, sdlog = 1
      ),
      "Upper truncation point is greater than D"
    )
    # x > D: cannot be observed under truncation; rejected.
    expect_error(
      dpcens(
        D + 1, plnorm,
        pwindow = pwindow, swindow = swindow, D = D,
        meanlog = 0, sdlog = 1
      ),
      "Upper truncation point is greater than D"
    )
    # Vectorised: error if any element violates the constraint.
    expect_error(
      dpcens(
        c(1, D + 1), plnorm,
        pwindow = pwindow, swindow = swindow, D = D,
        meanlog = 0, sdlog = 1
      ),
      "Upper truncation point is greater than D"
    )
  }
)

test_that("dprimarycensored errors for x below L = 0", {
  d <- -1
  pwindow <- 1
  swindow <- 0.5
  D <- 10

  expect_error(
    dpcens(
      d, plnorm,
      pwindow = pwindow, swindow = swindow, L = 0, D = D,
      meanlog = 0, sdlog = 1
    ),
    "values of x are below L"
  )
  expect_error(
    dpcens(
      c(8, d), plnorm,
      pwindow = pwindow, swindow = swindow, L = 0, D = D,
      meanlog = 0, sdlog = 1
    ),
    "values of x are below L"
  )
})

test_that("dprimarycensored returns non-negative values", {
  # Test case from issue #238: lognormal with specific parameters

  # that previously produced negative values due to floating-point precision
  pmf <- dpcens(
    x = seq(0, 29), plnorm, pwindow = 1, swindow = 1, D = 30,
    meanlog = 0.55, sdlog = 0.27
  )
  expect_true(all(pmf >= 0), info = "PMF should never be negative")

  # Also test with infinite D

  pmf_inf <- dpcens(
    x = seq(0, 29), plnorm, pwindow = 1, swindow = 1, D = Inf,
    meanlog = 0.55, sdlog = 0.27
  )
  expect_true(
    all(pmf_inf >= 0),
    info = "PMF with D=Inf should never be negative"
  )

  # Test with other distributions
  pmf_gamma <- dpcens(
    x = seq(0, 29), pgamma, pwindow = 1, swindow = 1, D = 30,
    shape = 2, rate = 0.5
  )
  expect_true(all(pmf_gamma >= 0), info = "Gamma PMF should never be negative")

  # Test with exponential growth primary distribution
  pmf_expgrowth <- dpcens(
    x = seq(0, 29), plnorm, pwindow = 1, swindow = 1, D = 30,
    dprimary = dexpgrowth,
    primary_args = list(r = 0.2),
    meanlog = 0.55, sdlog = 0.27
  )
  expect_true(
    all(pmf_expgrowth >= 0),
    info = "PMF with expgrowth primary should never be negative"
  )

  # Test with negative exponential growth primary distribution
  pmf_expgrowth_neg <- dpcens(
    x = seq(0, 29), plnorm, pwindow = 1, swindow = 1, D = 30,
    dprimary = dexpgrowth,
    primary_args = list(r = -0.2),
    meanlog = 0.55, sdlog = 0.27
  )
  expect_true(
    all(pmf_expgrowth_neg >= 0),
    info = "PMF with negative expgrowth primary should never be negative"
  )
})

test_that("dprimarycensored sums to 1 over [L, D) with L > 0", {
  pwindow <- 1
  D <- 10
  L <- 2
  pmf <- dpcens(
    L:(D - 1), plnorm, pwindow,
    D = D, L = L, meanlog = 1, sdlog = 1
  )
  expect_equal(sum(pmf), 1, tolerance = 1e-6)
})

test_that("dprimarycensored errors for x < L", {
  pwindow <- 1
  D <- 10
  L <- 2

  expect_error(
    dpcens(
      0:(L - 1), plnorm, pwindow,
      D = D, L = L, meanlog = 1, sdlog = 1
    ),
    "values of x are below L"
  )
})

test_that("dprimarycensored errors when L >= D", {
  expect_error(
    dpcens(5, plnorm, pwindow = 1, D = 10, L = 10, meanlog = 1, sdlog = 1),
    "L must be less than D"
  )
  expect_error(
    dpcens(5, plnorm, pwindow = 1, D = 10, L = 15, meanlog = 1, sdlog = 1),
    "L must be less than D"
  )
})

test_that("dprimarycensored default matches L = 0 for positive support", {
  pwindow <- 1
  D <- 10
  pmf_default <- dpcens(
    0:(D - 1), plnorm, pwindow,
    D = D, meanlog = 1, sdlog = 1
  )
  pmf_explicit <- dpcens(
    0:(D - 1), plnorm, pwindow,
    D = D, L = 0, meanlog = 1, sdlog = 1
  )
  expect_equal(pmf_default, pmf_explicit, tolerance = 1e-12)
})

test_that("dprimarycensored is consistent with pprimarycensored for L > 0", {
  pwindow <- 1
  D <- 10
  L <- 2
  x <- L:(D - 2)
  pmf <- dpcens(x, plnorm, pwindow, D = D, L = L, meanlog = 1, sdlog = 1)
  cdf_diff <- sapply(x, function(xi) {
    ppcens(
      xi + 1, plnorm, pwindow,
      D = D, L = L,
      meanlog = 1, sdlog = 1
    ) -
      ppcens(
        xi, plnorm, pwindow,
        D = D, L = L,
        meanlog = 1, sdlog = 1
      )
  })
  expect_equal(pmf, cdf_diff, tolerance = 1e-6)
})

test_that("dprimarycensored handles log probabilities with L > 0", {
  pwindow <- 1
  D <- 10
  L <- 2
  pmf <- dpcens(
    L:(D - 1), plnorm, pwindow,
    D = D, L = L, meanlog = 1, sdlog = 1
  )
  log_pmf <- dpcens(
    L:(D - 1), plnorm, pwindow,
    D = D, L = L, meanlog = 1, sdlog = 1, log = TRUE
  )
  expect_equal(exp(log_pmf), pmf, tolerance = 1e-6)
})

test_that("dprimarycensored works with L > 0 and D = Inf", {
  pwindow <- 1
  L <- 2
  pmf <- dpcens(
    L:10, plnorm, pwindow,
    D = Inf, L = L, meanlog = 1, sdlog = 1
  )
  expect_true(all(pmf >= 0))
  expect_lt(sum(pmf), 1) # Should be less than 1 since we're not summing all
})

test_that("dprimarycensored handles signed-support delays with negative L", {
  L <- -1
  D <- 2
  swindow <- 0.5
  x <- seq(L, D - swindow, by = swindow)
  pmf <- dpcens(
    x, pnorm,
    pwindow = 1, swindow = swindow, L = L, D = D, mean = 0, sd = 1
  )
  expect_true(all(pmf >= 0))
  expect_equal(sum(pmf), 1, tolerance = 1e-6)
  cdf_diff <- ppcens(
    x + swindow, pnorm,
    pwindow = 1, L = L, D = D, mean = 0, sd = 1
  ) - ppcens(
    x, pnorm,
    pwindow = 1, L = L, D = D, mean = 0, sd = 1
  )
  expect_equal(pmf, cdf_diff, tolerance = 1e-6)
})

test_that("dprimarycensored handles L = -Inf for signed-support delays", {
  # Exercises the `is.infinite(L) && is.infinite(D)` fast path that skips
  # truncation normalisation entirely.
  x <- -2:2
  pmf <- dpcens(
    x, pnorm,
    pwindow = 1, swindow = 1, L = -Inf, D = Inf, mean = 0, sd = 1
  )
  expect_true(all(pmf >= 0))
  expect_length(pmf, length(x))
})

test_that("dprimarycensored handles L not in unique_points", {
  # Use non-integer L that won't be in the computed unique_points
  # This exercises the branch where L needs to be computed separately
  pwindow <- 1
  swindow <- 1
  L <- 1.5
  D <- 10
  x <- 2:(D - 1)

  pmf <- dpcens(
    x, plnorm, pwindow,
    swindow = swindow, L = L, D = D, meanlog = 1, sdlog = 1
  )

  # Verify correctness by comparing to CDF differences
  cdf_upper <- ppcens(
    x + swindow, plnorm, pwindow,
    L = L, D = D, meanlog = 1, sdlog = 1
  )
  cdf_lower <- ppcens(
    x, plnorm, pwindow,
    L = L, D = D, meanlog = 1, sdlog = 1
  )
  expected_pmf <- cdf_upper - cdf_lower

  expect_true(all(pmf >= 0))
  expect_equal(pmf, expected_pmf, tolerance = 1e-10)
})

test_that("dprimarycensored at the d=1 step boundary returns pmf[1]", {
  # Regression for the lower-support boundary case under primary
  # censoring with pwindow = 1. F_cens(1) is zero because F_step is
  # zero on [0, 1) and the integral over a single jump point has zero
  # measure, so the unnormalised PMF at d = 1 reduces to F_cens(2),
  # which equals pmf[1] under the right-continuous step convention.
  boundaries <- 0:3
  pmf_vals <- c(0.25, 0.5, 0.25)
  pm <- dprimarycensored(
    x = 1, pdiscretestep,
    pwindow = 1, swindow = 1,
    boundaries = boundaries, pmf = pmf_vals
  )
  expect_equal(pm, pmf_vals[1], tolerance = 1e-10)
})

# Per-row reference: one scalar call for every observation
per_row <- function(f, x, ...) {
  settings <- lapply(list(...), rep_len, length(x))
  vapply(seq_along(x), function(i) {
    suppressMessages(do.call(f, c(list(x[i]), lapply(settings, `[[`, i))))
  }, numeric(1))
}

test_that("dprimarycensored with vector settings matches per-row calls", {
  set.seed(5)
  n <- 200
  x <- sample(0:12, n, replace = TRUE)
  pwindow <- sample(c(0, 1, 2), n, replace = TRUE)
  swindow <- sample(c(0, 0.5, 1, 2, Inf), n, replace = TRUE)
  L <- sample(c(-Inf, 0, -0.5), n, replace = TRUE)
  D <- sample(c(Inf, 25, 22.5, 40), n, replace = TRUE)
  cases <- list(
    function(x, ...) {
      dprimarycensored(x, pgamma, ..., shape = 2.5, scale = 1.5)
    },
    function(x, ...) {
      dprimarycensored(x, plnorm, ..., meanlog = 1.2, sdlog = 0.6)
    },
    function(x, ...) {
      dprimarycensored(
        x, pweibull, ...,
        dprimary = dexpgrowth, primary_args = list(r = 0.2),
        shape = 1.8, scale = 4
      )
    }
  )
  for (f in cases) {
    expect_equal(
      suppressMessages(f(x, pwindow, swindow, L, D)),
      per_row(f, x, pwindow, swindow, L, D),
      tolerance = 1e-10
    )
  }
  # A single x is recycled against vector settings
  expect_equal(
    cases[[1]](3, pwindow = c(1, 2), swindow = c(0, 1), D = c(10, Inf)),
    per_row(cases[[1]], c(3, 3), c(1, 2), c(0, 1), -Inf, c(10, Inf)),
    tolerance = 1e-10
  )
  # Only densities and no truncation, so no CDF points
  expect_equal(
    cases[[2]](1:4, pwindow = c(1, 2), swindow = 0),
    per_row(cases[[2]], 1:4, c(1, 2), 0),
    tolerance = 1e-10
  )
})

test_that("dprimarycensored evaluates each CDF point once per pwindow", {
  calls <- new.env(parent = emptyenv())
  calls$q <- list()
  local_mocked_bindings(
    pcens_cdf = function(object, q, pwindow, ...) {
      calls$q[[length(calls$q) + 1L]] <- q
      pgamma(q, shape = 2, scale = 1)
    }
  )
  set.seed(4)
  n <- 600
  suppressMessages(dprimarycensored(
    sample(0:20, n, replace = TRUE), pgamma,
    pwindow = sample.int(2, n, replace = TRUE),
    swindow = sample(c(0.5, 1, 2), n, replace = TRUE),
    L = sample(c(-Inf, 0), n, replace = TRUE),
    D = sample(c(Inf, 25, 40), n, replace = TRUE),
    shape = 2, scale = 1
  ))
  expect_length(calls$q, 2)
  for (q in calls$q) {
    expect_identical(anyDuplicated(q), 0L)
  }
})

test_that("dprimarycensored reports the row of invalid vector settings", {
  dens <- function(x, ...) dprimarycensored(x, pgamma, ..., shape = 2, rate = 1)
  expect_error(
    dens(c(1, 5, 3), D = c(10, 4, 10)),
    "Maximum x is 5 and D is 4 (row 2)",
    fixed = TRUE
  )
  expect_error(
    dens(c(1, 5, 3), L = c(0, 0, 4)),
    "Minimum x is 3 and L is 4 (row 3)",
    fixed = TRUE
  )
  expect_error(
    dens(c(1, 2, 3), L = c(0, 0, 4), D = c(5, 6, 4)),
    "First invalid row: 3 (L = 4, D = 4)",
    fixed = TRUE
  )
  expect_message(
    dens(c(1, 3.5, 3), D = c(10, 4, 10)),
    "It is 4.5 and D is 4 (row 2)",
    fixed = TRUE
  )
})
