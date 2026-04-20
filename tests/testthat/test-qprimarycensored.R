test_that("qprimarycensored handles boundary cases correctly", {
  # Default `L = -Inf` so the `p = 0` quantile is `-Inf`.
  q <- qpcens(c(0, 1), plnorm, pwindow = 1, meanlog = 1, sdlog = 1)
  expect_identical(q[1], -Inf)
  expect_true(is.na(q[2]))

  # Explicit `L = 0` reproduces the old default behaviour.
  q_l0 <- qpcens(c(0, 1), plnorm, pwindow = 1, L = 0, meanlog = 1, sdlog = 1)
  expect_identical(q_l0[1], 0)
  expect_true(is.na(q_l0[2]))
})

test_that("qprimarycensored returns monotonically increasing values", {
  probs <- seq(0.1, 0.9, by = 0.1)
  q <- qpcens(probs, plnorm, pwindow = 1, meanlog = 1, sdlog = 1)
  expect_true(all(diff(q) >= 0))
})

test_that("qprimarycensored and pprimarycensored are consistent", {
  probs <- c(0.25, 0.5, 0.75)
  q <- qpcens(probs, plnorm, pwindow = 1, meanlog = 2, sdlog = 1)
  p <- ppcens(q, plnorm, pwindow = 1, meanlog = 2, sdlog = 1)
  expect_equal(p, probs, tolerance = 1e-4)
})

test_that("qprimarycensored works with custom primary distributions", {
  probs <- c(0.25, 0.5, 0.75)
  pwindow <- 5

  q_uniform <- qpcens(
    probs,
    plnorm,
    pwindow,
    meanlog = 2,
    sdlog = 1
  )

  q_expgrowth <- qpcens(
    probs,
    plnorm,
    pwindow,
    dprimary = dexpgrowth,
    dprimary_args = list(r = 0.2),
    meanlog = 2,
    sdlog = 1
  )

  expect_false(all(q_uniform == q_expgrowth))
})

test_that("qprimarycensored handles truncation correctly", {
  probs <- c(0.25, 0.5, 0.75)
  D <- 10

  # Get quantiles with and without truncation
  q_untruncated <- qpcens(
    probs,
    plnorm,
    pwindow = 1,
    meanlog = 2,
    sdlog = 1
  )

  q_truncated <- qpcens(
    probs,
    plnorm,
    pwindow = 1,
    D = D,
    meanlog = 2,
    sdlog = 1
  )

  # Truncated quantiles should be less than or equal to D
  expect_true(all(q_truncated <= D))

  # Truncated quantiles should be less than or equal to untruncated
  expect_true(all(q_truncated <= q_untruncated))

  # Check consistency with pprimarycensored
  p_truncated <- ppcens(
    q_truncated,
    plnorm,
    pwindow = 1,
    D = D,
    meanlog = 2,
    sdlog = 1
  )
  expect_equal(p_truncated, probs, tolerance = 1e-4)
})


test_that("qprimarycensored handles left truncation boundary cases", {
  L <- 2
  q <- qpcens(
    c(0, 1), plnorm,
    pwindow = 1, D = 10, L = L, meanlog = 1, sdlog = 1
  )
  expect_identical(q[1], L)
  expect_true(is.na(q[2]))
})

test_that("qprimarycensored returns values >= L with L > 0", {
  probs <- seq(0.1, 0.9, by = 0.1)
  L <- 2
  D <- 10
  q <- qpcens(probs, plnorm, pwindow = 1, D = D, L = L, meanlog = 1, sdlog = 1)
  expect_true(all(q >= L))
  expect_true(all(q <= D))
})

test_that("qprimarycensored and pprimarycensored are consistent with L > 0", {
  probs <- c(0.25, 0.5, 0.75)
  L <- 2
  D <- 10
  q <- qpcens(probs, plnorm, pwindow = 1, D = D, L = L, meanlog = 1, sdlog = 1)
  p <- ppcens(q, plnorm, pwindow = 1, D = D, L = L, meanlog = 1, sdlog = 1)
  expect_equal(p, probs, tolerance = 1e-4)
})

test_that("qprimarycensored errors when L >= D", {
  expect_error(
    qpcens(0.5, plnorm, pwindow = 1, D = 10, L = 10, meanlog = 1, sdlog = 1),
    "L must be less than D"
  )
  expect_error(
    qpcens(0.5, plnorm, pwindow = 1, D = 10, L = 15, meanlog = 1, sdlog = 1),
    "L must be less than D"
  )
})

test_that("qprimarycensored inverts ppcens for signed-support delays", {
  probs <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  L <- -3
  D <- 3
  q <- qpcens(
    probs, pnorm,
    pwindow = 1, L = L, D = D, mean = 0, sd = 1
  )
  expect_true(all(q >= L & q <= D))
  back <- ppcens(
    q, pnorm,
    pwindow = 1, L = L, D = D, mean = 0, sd = 1
  )
  expect_equal(back, probs, tolerance = 1e-4)
})

test_that("qprimarycensored inverts pprimarycensored with L = -Inf", {
  probs <- c(0.25, 0.5, 0.75)
  q <- qpcens(
    probs, pnorm,
    pwindow = 1, L = -Inf, D = 3, mean = 0, sd = 1
  )
  expect_true(all(is.finite(q)))
  expect_true(all(q < 3))
  back <- ppcens(
    q, pnorm,
    pwindow = 1, L = -Inf, D = 3, mean = 0, sd = 1
  )
  expect_equal(back, probs, tolerance = 1e-4)
})

test_that(
  "qprimarycensored default matches L = 0 for positive-support delays",
  {
    # For positive-support delays `F_cens(0) = 0`, so the renormalisation
    # applied by an explicit `L = 0` is a numerical no-op relative to the
    # new `L = -Inf` default.
    probs <- c(0.25, 0.5, 0.75)
    q_default <- qpcens(
      probs, plnorm,
      pwindow = 1, D = 10, meanlog = 1, sdlog = 1
    )
    q_explicit <- qpcens(
      probs, plnorm,
      pwindow = 1, D = 10, L = 0, meanlog = 1, sdlog = 1
    )
    expect_equal(q_default, q_explicit, tolerance = 1e-6)
  }
)
