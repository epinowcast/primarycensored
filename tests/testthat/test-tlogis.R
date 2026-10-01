# Window, location and scale sets covering the midpoint before, inside and
# after the window, a sharp and a flat transition, and a non-zero min
tlogis_cases <- list(
  list(min = 0, max = 1, location = -0.5, scale = 0.2, label = "m below"),
  list(min = 0, max = 1, location = 0, scale = 0.1, label = "m at min"),
  list(min = 0, max = 1, location = 0.5, scale = 0.2, label = "m inside"),
  list(min = 0, max = 2, location = 1, scale = 1, label = "m mid, wide"),
  list(min = 0, max = 1, location = 1.5, scale = 0.1, label = "m above"),
  list(min = 0, max = 2, location = 40, scale = 5, label = "flat, far above"),
  list(min = 10, max = 20, location = 14, scale = 2, label = "non-zero min")
)

# Reference density from the definition, for parameters where the direct
# form is well conditioned
tlogis_naive_density <- function(x, tc) {
  mass <- plogis(tc$max, tc$location, tc$scale) -
    plogis(tc$min, tc$location, tc$scale)
  dlogis(x, tc$location, tc$scale) / mass
}

for (tc in tlogis_cases) {
  test_that(paste("dtlogis integrates to 1:", tc$label), {
    integral <- integrate(
      function(x) dtlogis(x, tc$min, tc$max, tc$location, tc$scale),
      tc$min, tc$max,
      rel.tol = 1e-10
    )$value
    expect_equal(integral, 1, tolerance = 1e-8)
  })

  test_that(paste("dtlogis matches the definition:", tc$label), {
    x <- seq(tc$min, tc$max, length.out = 9)
    expect_equal(
      dtlogis(x, tc$min, tc$max, tc$location, tc$scale),
      tlogis_naive_density(x, tc),
      tolerance = 1e-10
    )
  })

  test_that(paste("ptlogis matches the integral of dtlogis:", tc$label), {
    x <- seq(tc$min, tc$max, length.out = 7)[2:6]
    for (xi in x) {
      expect_equal(
        ptlogis(xi, tc$min, tc$max, tc$location, tc$scale),
        integrate(
          function(t) dtlogis(t, tc$min, tc$max, tc$location, tc$scale),
          tc$min, xi,
          rel.tol = 1e-10
        )$value,
        tolerance = 1e-8
      )
    }
  })

  test_that(paste("ptlogis boundary values and tails:", tc$label), {
    expect_identical(ptlogis(tc$min, tc$min, tc$max, tc$location, tc$scale), 0)
    expect_identical(ptlogis(tc$max, tc$min, tc$max, tc$location, tc$scale), 1)
    expect_identical(
      ptlogis(tc$min - 1, tc$min, tc$max, tc$location, tc$scale), 0
    )
    expect_identical(
      ptlogis(tc$max + 1, tc$min, tc$max, tc$location, tc$scale), 1
    )
    x <- seq(tc$min, tc$max, length.out = 7)[2:6]
    lower <- ptlogis(x, tc$min, tc$max, tc$location, tc$scale)
    upper <- ptlogis(
      x, tc$min, tc$max, tc$location, tc$scale,
      lower.tail = FALSE
    )
    expect_equal(lower + upper, rep(1, length(x)), tolerance = 1e-12)
    expect_equal(
      ptlogis(x, tc$min, tc$max, tc$location, tc$scale, log.p = TRUE),
      log(lower),
      tolerance = 1e-12
    )
    expect_equal(
      ptlogis(
        x, tc$min, tc$max, tc$location, tc$scale,
        lower.tail = FALSE, log.p = TRUE
      ),
      log(upper),
      tolerance = 1e-12
    )
  })

  test_that(paste("rtlogis samples match ptlogis:", tc$label), {
    set.seed(20260930)
    n <- 20000
    samples <- rtlogis(n, tc$min, tc$max, tc$location, tc$scale)
    expect_length(samples, n)
    expect_true(all(samples >= tc$min & samples <= tc$max))
    ks <- suppressWarnings(ks.test(
      samples,
      function(q) ptlogis(q, tc$min, tc$max, tc$location, tc$scale)
    ))
    expect_gt(ks$p.value, 1e-3)
    # The sample mean is within five standard errors of the integral
    mean_expected <- integrate(
      function(x) x * dtlogis(x, tc$min, tc$max, tc$location, tc$scale),
      tc$min, tc$max,
      rel.tol = 1e-10
    )$value
    sd_expected <- sqrt(integrate(
      function(x) {
        (x - mean_expected)^2 *
          dtlogis(x, tc$min, tc$max, tc$location, tc$scale)
      },
      tc$min, tc$max,
      rel.tol = 1e-10
    )$value)
    expect_lt(abs(mean(samples) - mean_expected), 5 * sd_expected / sqrt(n))
  })

  test_that(paste("dtlogis is zero outside the window:", tc$label), {
    expect_identical(
      dtlogis(
        c(tc$min - 0.1, tc$max + 0.1), tc$min, tc$max, tc$location, tc$scale
      ),
      c(0, 0)
    )
    expect_identical(
      dtlogis(
        c(tc$min - 0.1, tc$max + 0.1), tc$min, tc$max, tc$location, tc$scale,
        log = TRUE
      ),
      c(-Inf, -Inf)
    )
  })
}

test_that("the tlogis functions are stable far from the location", {
  # 1 - plogis() underflows here, the window mass is about exp(-300)
  loc <- -6
  sc <- 0.02
  expect_true(is.finite(dtlogis(0.3, 0, 1, loc, sc, log = TRUE)))
  expect_equal(
    integrate(
      function(x) dtlogis(x, 0, 1, loc, sc), 0, 1,
      rel.tol = 1e-10
    )$value,
    1,
    tolerance = 1e-8
  )
  # The window density is that of an exponential truncated at 1, rate 1 / sc
  expect_equal(
    dtlogis(0.3, 0, 1, loc, sc),
    dexp(0.3, 1 / sc) / pexp(1, 1 / sc),
    tolerance = 1e-8
  )
  expect_equal(
    ptlogis(0.3, 0, 1, loc, sc),
    pexp(0.3, 1 / sc) / pexp(1, 1 / sc),
    tolerance = 1e-8
  )
  set.seed(1)
  samples <- rtlogis(1000, 0, 1, loc, sc)
  expect_true(all(is.finite(samples) & samples >= 0 & samples <= 1))
  # The mirrored window, far below the location
  expect_equal(
    dtlogis(0.7, 0, 1, 7, sc),
    dexp(0.3, 1 / sc) / pexp(1, 1 / sc),
    tolerance = 1e-8
  )
})

test_that("tlogis approaches the uniform for a large scale", {
  expect_equal(
    dtlogis(c(0.2, 0.8), 0, 1, 0.5, 1e5), c(1, 1),
    tolerance = 1e-5
  )
  expect_equal(ptlogis(0.3, 0, 1, 0.5, 1e5), 0.3, tolerance = 1e-5)
})

test_that("tlogis functions recycle arguments and keep names", {
  expect_length(dtlogis(c(0.1, 0.5, 0.9), 0, 1, c(0, 0.5, 1), 0.3), 3)
  expect_length(ptlogis(c(0.1, 0.5, 0.9), 0, 1, c(0, 0.5, 1), 0.3), 3)
  expect_identical(attr(dtlogis, "name"), "dtlogis")
  expect_identical(attr(ptlogis, "name"), "ptlogis")
  expect_length(rtlogis(0, 0, 1, 0.5, 1), 0L)
})

test_that("tlogis functions check their arguments", {
  expect_error(dtlogis(0.5, 0, 1, 0.5, 0), "scale")
  expect_error(dtlogis(0.5, 0, 1, 0.5, -1), "scale")
  expect_error(ptlogis(0.5, 0, 1, 0.5, 0), "scale")
  expect_error(rtlogis(3, 0, 1, 0.5, 0), "scale")
  expect_error(dtlogis(0.5, 1, 1, 0.5, 1), "min")
  expect_error(ptlogis(0.5, 2, 1, 0.5, 1), "min")
  expect_error(rtlogis(3, 1, 0, 0.5, 1), "min")
  expect_error(dtlogis(0.5, 0, 1, Inf, 1), "location")
})

test_that("tlogis is registered as primary distribution 3", {
  expect_identical(pcd_dist_name("tlogis", type = "primary"), "dtlogis")
  expect_identical(
    pcd_dist_name("truncated logistic", type = "primary"), "dtlogis"
  )
  expect_identical(pcd_stan_dist_id("tlogis", type = "primary"), 3L)
  expect_identical(
    pcd_stan_dist_id("truncated logistic", type = "primary"), 3L
  )
  registry <- primarycensored::pcd_primary_distributions
  row <- registry[registry$name == "tlogis", ]
  expect_identical(row$pprimary, "ptlogis")
})

test_that("tlogis works as the primary of pcens objects", {
  obj <- new_pcens(
    pdist = pgamma, dprimary = dtlogis,
    primary_args = list(location = 0.5, scale = 0.2),
    shape = 2, rate = 1
  )
  expect_s3_class(obj, "pcens_pgamma_dtlogis")
  expect_identical(obj$pprimary, ptlogis)
  expect_identical(
    attr(obj$pprimary, "name"), "ptlogis"
  )
  expect_error(
    new_pcens(
      pdist = pgamma, dprimary = dtlogis, pprimary = pexpgrowth,
      primary_args = list(location = 0.5, scale = 0.2),
      shape = 2, rate = 1
    ),
    "different distributions"
  )
  expect_null(check_dprimary(
    dtlogis,
    pwindow = 2, dprimary_args = list(location = 1, scale = 0.3)
  ))
})

test_that("rprimarycensored and pprimarycensored run with tlogis", {
  set.seed(3)
  samples <- rprimarycensored(
    2000, rgamma, rprimary = rtlogis, pwindow = 2, swindow = 0,
    rprimary_args = list(location = 0.5, scale = 0.3),
    shape = 3, rate = 1
  )
  expect_length(samples, 2000)
  empirical <- ecdf(samples)
  q <- c(1, 3, 6)
  theoretical <- pprimarycensored(
    q, pgamma, dprimary = dtlogis, pwindow = 2,
    primary_args = list(location = 0.5, scale = 0.3),
    shape = 3, rate = 1
  )
  expect_equal(empirical(q), theoretical, tolerance = 0.04)
})

test_that("non-parametric delays use the truncated logistic CDF", {
  boundaries <- c(0, 1, 3, 6, 10)
  pmf <- c(0.2, 0.3, 0.35, 0.15)
  obj <- new_pcens(
    pdist = pdiscretestep, dprimary = dtlogis,
    primary_args = list(location = 0.5, scale = 0.3),
    boundaries = boundaries, pmf = pmf
  )
  q <- c(0.5, 1.5, 3, 5, 8, 12)
  expect_equal(
    pcens_cdf(obj, q, 2),
    pcens_cdf(obj, q, 2, use_numeric = TRUE),
    tolerance = 1e-6
  )
})

narrow_delays <- list(
  list(
    label = "exponential rate 1.5", pdist = pexp, args = list(rate = 1.5),
    positive = TRUE
  ),
  list(
    label = "gamma shape 2", pdist = pgamma,
    args = list(shape = 2, rate = 1), positive = TRUE
  ),
  list(
    label = "lognormal", pdist = plnorm,
    args = list(meanlog = 1, sdlog = 0.5), positive = TRUE
  ),
  list(
    label = "weibull", pdist = pweibull,
    args = list(shape = 2, scale = 2), positive = TRUE
  ),
  list(
    label = "normal", pdist = pnorm, args = list(mean = 3, sd = 2),
    positive = FALSE
  )
)

# Relative difference, absolute below 1e-10 where `stats::integrate()` has
# its absolute tolerance
narrow_diff <- function(actual, expected) {
  max(abs(actual - expected) / pmax(expected, 1e-10))
}

test_that("the numerical CDF resolves a narrow truncated logistic primary", {
  pwindow <- 2
  q <- c(0.5, 1, 3, 10)
  for (delay in narrow_delays) {
    cdf <- function(x) do.call(delay$pdist, c(list(x), delay$args))
    for (location in c(-0.5, 0.7, 2.5)) {
      for (scale in c(0.02, 0.005, 0.001)) {
        obj <- tlogis_object(delay, location, scale)
        expected <- tlogis_reference(
          q, pwindow, location, scale, cdf, delay$positive
        )
        label <- tlogis_label(delay, pwindow, location, scale)
        expect_lt(
          narrow_diff(pcens_cdf(obj, q, pwindow), expected), 1e-6,
          label = label
        )
        expect_lt(
          narrow_diff(
            pcens_cdf(obj, q, pwindow, use_numeric = TRUE), expected
          ), 1e-6,
          label = paste("numeric:", label)
        )
      }
    }
  }
})

test_that("the numerical CDF resolves a named wrapper of dtlogis", {
  # A named wrapper gets the tlogis classes and so the same break points
  wrapper <- function(x, min = 0, max = 1, location = 0, scale = 1,
                      log = FALSE) {
    dtlogis(x, min, max, location, scale, log)
  }
  named <- add_name_attribute(wrapper, "dtlogis")
  pwindow <- 2
  q <- c(0.72, 1, 1.5, 3)
  cdf <- function(x) pgamma(x, shape = 2, rate = 1)
  for (scale in c(1e-3, 1e-4)) {
    obj <- new_pcens(
      pdist = pgamma, dprimary = named,
      primary_args = list(location = 0.7, scale = scale),
      shape = 2, rate = 1
    )
    expect_s3_class(obj, "pcens_pgamma_dtlogis")
    expected <- tlogis_reference(q, pwindow, 0.7, scale, cdf)
    expect_lt(
      narrow_diff(pcens_cdf(obj, q, pwindow), expected), 1e-6,
      label = paste("scale", scale)
    )
    expect_lt(
      narrow_diff(pcens_cdf(obj, q, pwindow, use_numeric = TRUE), expected),
      1e-6,
      label = paste("numeric, scale", scale)
    )
  }
})
