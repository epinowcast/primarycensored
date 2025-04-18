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
