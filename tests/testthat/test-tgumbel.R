# Parameter sets covering a peak inside the window, a window far above the
# peak (declining density), far below it (growing density), a non-zero min
# and a wide scale.
tgumbel_cases <- list(
  list(min = 0, max = 1, mu = 0.5, beta = 0.3, label = "peak inside"),
  list(min = 0, max = 2, mu = -0.5, beta = 0.2, label = "declining"),
  list(min = 0, max = 2, mu = 1.5, beta = 0.1, label = "growing"),
  list(min = 0, max = 1, mu = 0, beta = 1, label = "wide scale"),
  list(min = 5, max = 7, mu = 6, beta = 0.5, label = "non-zero min"),
  list(min = -2, max = 1, mu = -0.5, beta = 0.4, label = "negative min")
)

for (tc in tgumbel_cases) {
  test_that(paste("dtgumbel integrates to 1:", tc$label), {
    integral <- integrate(
      function(x) dtgumbel(x, tc$min, tc$max, tc$mu, tc$beta),
      tc$min, tc$max,
      rel.tol = 1e-10
    )$value
    expect_equal(integral, 1, tolerance = 1e-9)
  })

  test_that(paste("ptgumbel matches the integral of dtgumbel:", tc$label), {
    x <- seq(tc$min, tc$max, length.out = 7)[2:6]
    from_integral <- vapply(x, function(xx) {
      integrate(
        function(t) dtgumbel(t, tc$min, tc$max, tc$mu, tc$beta),
        tc$min, xx,
        rel.tol = 1e-12
      )$value
    }, numeric(1))
    expect_equal(
      ptgumbel(x, tc$min, tc$max, tc$mu, tc$beta), from_integral,
      tolerance = 1e-9
    )
  })

  test_that(paste("ptgumbel tails and log scale are consistent:", tc$label), {
    x <- seq(tc$min, tc$max, length.out = 9)
    lower <- ptgumbel(x, tc$min, tc$max, tc$mu, tc$beta)
    upper <- ptgumbel(
      x, tc$min, tc$max, tc$mu, tc$beta,
      lower.tail = FALSE
    )
    expect_equal(lower + upper, rep(1, length(x)), tolerance = 1e-12)
    # Compare on the log scale only where the probability does not underflow
    keep <- lower > 1e-300 & upper > 1e-300
    expect_equal(
      exp(ptgumbel(x, tc$min, tc$max, tc$mu, tc$beta, log.p = TRUE))[keep],
      lower[keep],
      tolerance = 1e-12
    )
    expect_equal(
      exp(ptgumbel(
        x, tc$min, tc$max, tc$mu, tc$beta,
        lower.tail = FALSE, log.p = TRUE
      ))[keep],
      upper[keep],
      tolerance = 1e-12
    )
    expect_identical(lower[1], 0)
    expect_identical(lower[length(x)], 1)
  })

  test_that(paste("dtgumbel is the derivative of ptgumbel:", tc$label), {
    x <- seq(tc$min, tc$max, length.out = 7)[2:6]
    h <- 1e-6 * (tc$max - tc$min)
    numeric_derivative <- (
      ptgumbel(x + h, tc$min, tc$max, tc$mu, tc$beta) -
        ptgumbel(x - h, tc$min, tc$max, tc$mu, tc$beta)
    ) / (2 * h)
    expect_equal(
      dtgumbel(x, tc$min, tc$max, tc$mu, tc$beta), numeric_derivative,
      tolerance = 1e-6
    )
    dens <- dtgumbel(x, tc$min, tc$max, tc$mu, tc$beta)
    keep <- dens > 1e-300
    expect_equal(
      exp(dtgumbel(x, tc$min, tc$max, tc$mu, tc$beta, log = TRUE))[keep],
      dens[keep],
      tolerance = 1e-12
    )
  })

  test_that(paste("rtgumbel samples match ptgumbel:", tc$label), {
    set.seed(123)
    samples <- rtgumbel(20000, tc$min, tc$max, tc$mu, tc$beta)
    expect_true(all(samples >= tc$min))
    expect_true(all(samples <= tc$max))
    ks <- suppressWarnings(stats::ks.test(
      samples,
      function(x) ptgumbel(x, tc$min, tc$max, tc$mu, tc$beta)
    ))
    expect_gt(ks$p.value, 0.001)
  })
}

test_that("dtgumbel is zero and ptgumbel saturates outside the window", {
  expect_identical(dtgumbel(c(-0.1, 1.1), 0, 1, 0.5, 0.3), c(0, 0))
  expect_identical(ptgumbel(c(-0.1, 1.1), 0, 1, 0.5, 0.3), c(0, 1))
  expect_identical(
    ptgumbel(c(-0.1, 1.1), 0, 1, 0.5, 0.3, lower.tail = FALSE), c(1, 0)
  )
  expect_identical(
    dtgumbel(c(-0.1, 1.1), 0, 1, 0.5, 0.3, log = TRUE), c(-Inf, -Inf)
  )
  expect_identical(dtgumbel(NA_real_, 0, 1, 0.5, 0.3), NA_real_)
  expect_identical(ptgumbel(NA_real_, 0, 1, 0.5, 0.3), NA_real_)
})

test_that("dtgumbel tends to the exponentially tilted window", {
  # With mu far below the window, G(z) - G(0) is proportional to
  # exp(-z / beta) - 1, the exponentially tilted CDF with r = -1 / beta
  x <- seq(0, 1, length.out = 6)
  expect_equal(
    dtgumbel(x, 0, 1, mu = -8, beta = 0.5),
    dexpgrowth(x, 0, 1, r = -2),
    tolerance = 1e-6
  )
  expect_equal(
    ptgumbel(x, 0, 1, mu = -8, beta = 0.5),
    pexpgrowth(x, 0, 1, r = -2),
    tolerance = 1e-6
  )
})

test_that("ptgumbel tends to the uniform window for a wide scale", {
  x <- seq(0, 1, length.out = 6)
  expect_equal(
    ptgumbel(x, 0, 1, mu = 0.5, beta = 1e4), x,
    tolerance = 1e-3
  )
})

test_that("tgumbel functions stay finite at extreme parameters", {
  x <- seq(0, 1, length.out = 6)
  for (params in list(
    c(mu = 30, beta = 0.05), c(mu = -30, beta = 0.05),
    c(mu = 0.5, beta = 1e-3), c(mu = 0.5, beta = 1e3)
  )) {
    dens <- dtgumbel(x, 0, 1, params[["mu"]], params[["beta"]])
    cdf <- ptgumbel(x, 0, 1, params[["mu"]], params[["beta"]])
    samples <- rtgumbel(50, 0, 1, params[["mu"]], params[["beta"]])
    expect_true(all(is.finite(dens)), info = toString(params))
    expect_true(all(is.finite(cdf)), info = toString(params))
    expect_true(all(diff(cdf) >= 0), info = toString(params))
    expect_true(all(is.finite(samples)), info = toString(params))
  }
})

test_that("tgumbel functions reject invalid parameters", {
  expect_error(dtgumbel(0.5, 0, 1, mu = 0, beta = 0), "beta")
  expect_error(ptgumbel(0.5, 0, 1, mu = 0, beta = -1), "beta")
  expect_error(rtgumbel(2, 0, 1, mu = Inf, beta = 1), "mu")
  expect_error(dtgumbel(0.5, 1, 1, mu = 0, beta = 1), "min")
  expect_error(dtgumbel(0.5, 0, 1, mu = NA_real_, beta = 1), "mu")
})

test_that("tgumbel is registered as a primary distribution", {
  expect_identical(pcd_stan_dist_id("tgumbel", type = "primary"), 4L)
  expect_identical(pcd_stan_dist_id("truncated gumbel", "primary"), 4L)
  expect_identical(pcd_dist_name("tgumbel", type = "primary"), "dtgumbel")
  expect_identical(attr(dtgumbel, "name"), "dtgumbel")
  expect_identical(attr(ptgumbel, "name"), "ptgumbel")
  obj <- new_pcens(
    pnorm, dtgumbel,
    primary_args = list(mu = 0.5, beta = 0.3), mean = 3, sd = 2
  )
  expect_identical(class(obj)[1], "pcens_pnorm_dtgumbel")
  expect_identical(obj$pprimary, ptgumbel)
})

test_that("a truncated Gumbel primary works through the public functions", {
  set.seed(1)
  draws <- rprimarycensored(
    2000, rnorm,
    rprimary = rtgumbel,
    rprimary_args = list(mu = 0.5, beta = 0.3),
    pwindow = 2, swindow = 0, mean = 5, sd = 1
  )
  expect_length(draws, 2000)
  ks <- suppressWarnings(stats::ks.test(
    draws,
    function(x) {
      pprimarycensored(
        x, pnorm,
        pwindow = 2, dprimary = dtgumbel,
        primary_args = list(mu = 0.5, beta = 0.3), mean = 5, sd = 1
      )
    }
  ))
  expect_gt(ks$p.value, 0.001)
})
