test_that("pcens_quantile handles boundary cases correctly", {
  pdist <- pgamma
  dprimary <- dunif
  obj <- new_pcens(
    pdist,
    dprimary,
    list(),
    shape = 2,
    rate = 1
  )

  expect_identical(pcens_quantile(obj, p = 0, pwindow = 1), 0)
  expect_true(is.na(pcens_quantile(obj, p = 1, pwindow = 1)))
})

test_that("pcens_quantile and pcens_cdf are consistent", {
  pdist <- pgamma
  dprimary <- dunif

  shapes <- c(2, 4)
  rates <- c(2, 3)
  pwindows <- 1
  probs <- c(0.25, 0.5, 0.75)

  for (shape in shapes) {
    for (rate in rates) {
      for (pwindow in pwindows) {
        obj <- new_pcens(
          pdist,
          dprimary,
          list(),
          shape = shape,
          rate = rate
        )

        # Get quantiles
        quantiles <- pcens_quantile(obj, p = probs, pwindow = pwindow)

        # Check CDF at those quantiles matches original probabilities
        cdfs <- pcens_cdf(obj, q = quantiles, pwindow = pwindow)

        expect_equal(
          cdfs,
          probs,
          tolerance = 1e-4,
          info = sprintf(
            "Mismatch for shape = %s, rate = %s, pwindow = %s",
            shape,
            rate,
            pwindow
          )
        )
      }
    }
  }
})

test_that("pcens_quantile returns monotonically increasing values", {
  pdist <- pgamma
  dprimary <- dunif
  obj <- new_pcens(
    pdist,
    dprimary,
    list(),
    shape = 2,
    rate = 1
  )

  probs <- seq(0.1, 0.9, by = 0.1)
  quantiles <- pcens_quantile(obj, p = probs, pwindow = 1)

  expect_gt(min(diff(quantiles)), 0)
})

test_that("pcens_quantile works with different initial values", {
  pdist <- pgamma
  dprimary <- dunif
  obj <- new_pcens(
    pdist,
    dprimary,
    list(),
    shape = 2,
    rate = 1
  )

  p <- 0.5
  pwindow <- 1

  result1 <- pcens_quantile(obj, p = p, pwindow = pwindow, init = 5)
  result2 <- pcens_quantile(obj, p = p, pwindow = pwindow, init = 10)

  expect_equal(result1, result2, tolerance = 1e-4)
})

test_that("pcens_quantile respects tolerance parameter", {
  pdist <- pgamma
  dprimary <- dunif
  obj <- new_pcens(
    pdist,
    dprimary,
    list(),
    shape = 2,
    rate = 1
  )

  p <- 0.5
  pwindow <- 1

  # Compute CDF at quantile for different tolerances
  q_loose <- pcens_quantile(obj, p = p, pwindow = pwindow, tol = 1e-4)
  cdf_loose <- pcens_cdf(obj, q = q_loose, pwindow = pwindow)

  q_tight <- pcens_quantile(obj, p = p, pwindow = pwindow, tol = 1e-8)
  cdf_tight <- pcens_cdf(obj, q = q_tight, pwindow = pwindow)

  # Tight tolerance should give closer match to target probability
  expect_lte(abs(cdf_tight - p), abs(cdf_loose - p))
})
