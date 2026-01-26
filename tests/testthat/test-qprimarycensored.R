test_that("qprimarycensored handles boundary cases correctly", {
  q <- qpcens(c(0, 1), plnorm, pwindow = 1, meanlog = 1, sdlog = 1)
  expect_identical(q[1], 0)
  expect_true(is.na(q[2]))
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

# Left truncation (L parameter) tests

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

test_that("qprimarycensored errors when L < 0", {
  expect_error(
    qpcens(0.5, plnorm, pwindow = 1, D = 10, L = -1, meanlog = 1, sdlog = 1),
    "L must be non-negative"
  )
})

test_that("qprimarycensored with L = 0 matches default behaviour", {
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
})
