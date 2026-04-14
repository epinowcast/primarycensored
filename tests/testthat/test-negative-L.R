test_that("pprimarycensored accepts negative L", {
  expect_no_error(
    ppcens(0, pnorm, pwindow = 1, L = -5, D = 10, mean = 0, sd = 1)
  )
})

test_that("pprimarycensored matches full-R convolution for pnorm", {
  q_vals <- seq(-3, 3, by = 0.5)
  reference <- vapply(q_vals, function(q) {
    stats::integrate(
      function(p) {
        stats::pnorm(q - p, mean = 0, sd = 1) * stats::dunif(p, 0, 1)
      },
      0, 1
    )$value
  }, numeric(1))
  result <- ppcens(
    q_vals, pnorm,
    pwindow = 1, L = -Inf, D = Inf, mean = 0, sd = 1
  )
  expect_equal(result, reference, tolerance = 1e-6)
})

test_that("pprimarycensored with L = -3, D = 3, pnorm is a proper CDF", {
  q_vals <- seq(-3, 3, by = 0.25)
  cdf <- ppcens(
    q_vals, pnorm,
    pwindow = 1, L = -3, D = 3, mean = 0, sd = 1
  )
  f_L <- ppcens(-3, pnorm, pwindow = 1, L = -3, D = 3, mean = 0, sd = 1)
  f_D <- ppcens(3, pnorm, pwindow = 1, L = -3, D = 3, mean = 0, sd = 1)
  expect_equal(f_L, 0, tolerance = 1e-6)
  expect_equal(f_D, 1, tolerance = 1e-6)
  expect_true(all(diff(cdf) >= 0))
  expect_true(all(cdf >= 0 & cdf <= 1))
})

test_that("pprimarycensored returns proper truncated value for pnorm, L = 0", {
  raw_0 <- stats::integrate(
    function(p) stats::pnorm(0 - p), 0, 1
  )$value
  raw_q <- stats::integrate(
    function(p) stats::pnorm(0.5 - p), 0, 1
  )$value
  reference <- (raw_q - raw_0) / (1 - raw_0)
  result <- ppcens(
    0.5, pnorm,
    pwindow = 1, L = 0, D = Inf, mean = 0, sd = 1
  )
  expect_equal(result, reference, tolerance = 1e-6)
  f_at_0 <- ppcens(
    0, pnorm,
    pwindow = 1, L = 0, D = Inf, mean = 0, sd = 1
  )
  expect_equal(f_at_0, 0, tolerance = 1e-6)
})

test_that("dprimarycensored handles x with 0 and negative values for pnorm", {
  x <- -2:2
  swindow <- 1
  d_vals <- dpcens(
    x, pnorm,
    pwindow = 1, swindow = swindow, L = -3, D = 3, mean = 0, sd = 1
  )
  expect_true(all(d_vals >= 0))
  upper <- ppcens(
    max(x) + swindow, pnorm,
    pwindow = 1, L = -3, D = 3, mean = 0, sd = 1
  )
  lower <- ppcens(
    min(x), pnorm,
    pwindow = 1, L = -3, D = 3, mean = 0, sd = 1
  )
  expect_equal(sum(d_vals), upper - lower, tolerance = 1e-6)
})

test_that("dprimarycensored matches diff of pprimarycensored with negative L", {
  x <- seq(-2, 2, by = 1)
  d_pcens <- dpcens(
    x, pnorm,
    pwindow = 1, swindow = 1, L = -3, D = 3, mean = 0, sd = 1
  )
  upper <- ppcens(
    x + 1, pnorm,
    pwindow = 1, L = -3, D = 3, mean = 0, sd = 1
  )
  lower <- ppcens(
    x, pnorm,
    pwindow = 1, L = -3, D = 3, mean = 0, sd = 1
  )
  norm_upper <- ppcens(
    3, pnorm,
    pwindow = 1, L = -3, D = 3, mean = 0, sd = 1
  )
  norm_lower <- ppcens(
    -3, pnorm,
    pwindow = 1, L = -3, D = 3, mean = 0, sd = 1
  )
  reference <- (upper - lower) / (norm_upper - norm_lower)
  expect_equal(d_pcens, reference, tolerance = 1e-6)
})

test_that("qprimarycensored inverts pprimarycensored with negative L", {
  probs <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  q <- qpcens(
    probs, pnorm,
    pwindow = 1, L = -3, D = 3, mean = 0, sd = 1
  )
  expect_true(all(q >= -3 & q <= 3))
  back <- ppcens(
    q, pnorm,
    pwindow = 1, L = -3, D = 3, mean = 0, sd = 1
  )
  expect_equal(back, probs, tolerance = 1e-4)
})

test_that("rprimarycensored produces samples in [L, D) for negative L", {
  set.seed(42)
  samples <- rpcens(
    n = 5000, rnorm,
    pwindow = 1, swindow = 1, L = -3, D = 3, mean = 0, sd = 1
  )
  expect_length(samples, 5000)
  expect_true(all(samples >= -3 & samples < 3))
  expect_lt(abs(mean(samples)), 0.2)
  expect_lt(abs(stats::sd(samples) - 1), 0.2)
})

test_that("pprimarycensored with L = 0 matches reference for plnorm", {
  q_vals <- seq(0.1, 5, by = 0.5)
  reference <- vapply(q_vals, function(q) {
    stats::integrate(
      function(p) {
        stats::plnorm(q - p, meanlog = 0, sdlog = 1) *
          stats::dunif(p, 0, 1)
      },
      0, 1
    )$value
  }, numeric(1))
  result <- ppcens(
    q_vals, plnorm,
    pwindow = 1, L = 0, D = Inf, meanlog = 0, sdlog = 1
  )
  expect_equal(result, reference, tolerance = 1e-6)
})

test_that("pprimarycensored still errors when L >= D", {
  expect_error(
    ppcens(
      5, pnorm,
      pwindow = 1, L = 10, D = 10, mean = 0, sd = 1
    ),
    "L must be less than D"
  )
  expect_error(
    ppcens(
      5, pnorm,
      pwindow = 1, L = 15, D = 10, mean = 0, sd = 1
    ),
    "L must be less than D"
  )
})
