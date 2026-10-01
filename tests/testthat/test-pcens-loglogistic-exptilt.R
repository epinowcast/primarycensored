# Log-logistic delays with an exponentially tilted primary. The transform is
# the series F(t) sum_n (xi t)^n / n! r_{n/shape}(A), used where |xi| t is at
# most 10, and the numerical method is used beyond.

families <- loglogistic_families()
pwindows <- c(0.5, 1, 2, 7)
rhos <- c(
  -1, -0.5, -0.05, -3e-4, -1e-4, -1e-5, -1e-8, 1e-8, 1e-5, 1e-4, 3e-4, 0.05,
  0.5, 1
)
series_limit <- 10

loglogistic_object <- function(family, rho) {
  exptilt_object(family, rho)
}

test_that("log-logistic delays dispatch to the tilted primary method", {
  obj <- loglogistic_object(families[[2]], 0.2)
  expect_s3_class(obj, "pcens_pllogis_dexpgrowth")
  expect_false(
    is.null(utils::getS3method(
      "pcens_cdf", "pcens_pllogis_dexpgrowth",
      optional = TRUE
    ))
  )
  expect_identical(.pcens_tilt_lower(obj), 0)
  expect_true(.pcens_tilt_available(obj, -0.3))
  # The transform over a finite range exists for any tilt
  expect_true(.pcens_tilt_available(obj, 5))
})

test_that("tilt transforms match numerical integration of the density", {
  ts <- c(1e-3, 0.4, 1, 3.5, 8, 12)
  xis <- c(0, -1, -0.25, -0.05, 0.05, 0.25, 0.6)
  for (family in families) {
    obj <- loglogistic_object(family, 0.2)
    shape <- family$args$shape
    scale <- family$args$scale
    for (xi in xis) {
      keep <- abs(xi) * ts <= series_limit
      t <- ts[keep]
      lower <- .pcens_tilt_transform(obj, t, xi)
      expected <- vapply(
        t, loglogistic_transform_ref, numeric(1), xi, shape, scale
      )
      expect_equal(
        exp(lower), expected,
        tolerance = 1e-9,
        info = paste(family$label, "xi", xi)
      )
      upper <- .pcens_tilt_transform(obj, t, xi, upper = TRUE)
      if (xi == 0) {
        expect_equal(
          exp(lower) + exp(upper), rep(1, length(t)),
          tolerance = 1e-14
        )
      } else {
        # The upper transform has no closed form and is not available
        expect_true(all(is.nan(upper)))
      }
    }
  }
})

test_that("tilt transforms are accurate deep in the lower tail", {
  obj <- loglogistic_object(families[[4]], 0.2)
  t <- c(1e-6, 1e-3, 0.05)
  lower <- .pcens_tilt_transform(obj, t, -0.3)
  expected <- vapply(
    t, loglogistic_transform_ref, numeric(1), -0.3, 4.5, 8
  )
  # log T follows 4.5 log(t / 8) far below the median
  expect_equal(lower, log(expected), tolerance = 1e-10)
  expect_equal(
    .pcens_tilt_transform(obj, 1e-300, -0.3),
    4.5 * (log(1e-300) - log(8)),
    tolerance = 1e-12
  )
})

test_that("tilt transforms are 0 below the support", {
  obj <- loglogistic_object(families[[2]], 0.2)
  expect_identical(.pcens_tilt_transform(obj, c(-2, 0), -0.3), c(-Inf, -Inf))
  expect_identical(.pcens_tilt_transform(obj, c(-2, 0), 0), c(-Inf, -Inf))
  expect_identical(
    .pcens_tilt_transform(obj, c(-2, 0), 0, upper = TRUE), c(0, 0)
  )
  expect_true(all(is.nan(
    .pcens_tilt_transform(obj, c(-2, 0), -0.3, upper = TRUE)
  )))
})

test_that("the transform is NaN beyond its limit", {
  obj <- loglogistic_object(families[[2]], 0.2)
  # The limit is |xi| t of 10, so 20 for xi of -0.5
  inside <- .pcens_tilt_transform(obj, c(19, 20), -0.5)
  expect_true(all(is.finite(inside)))
  expect_true(all(is.nan(.pcens_tilt_transform(obj, c(20.5, 30), -0.5))))
  expect_true(all(is.finite(.pcens_tilt_transform(obj, c(30, 1e3), 0))))
  # A shape below 0.2 has no series
  slow <- exptilt_object(
    list(pdist = pllogis_test, args = list(shape = 0.15, scale = 2)), 0.2
  )
  expect_true(all(is.nan(.pcens_tilt_transform(slow, c(1, 3), -0.5))))
})

test_that("tilt moments match numerical integration", {
  ts <- c(1e-3, 0.4, 1, 3.5, 12, 60)
  for (family in families) {
    obj <- loglogistic_object(family, 0.2)
    shape <- family$args$shape
    scale <- family$args$scale
    moments <- .pcens_tilt_moments(obj, c(-1, 0, ts))
    expect_identical(unname(moments[1:2, ]), matrix(-Inf, 2, 3))
    for (i in seq_along(ts)) {
      g <- vapply(1:3, function(k) {
        stats::integrate(
          function(u) (ts[i] - u)^k * dllogis_test(u, shape, scale),
          0, ts[i],
          rel.tol = 1e-12, abs.tol = 0, subdivisions = 2000L
        )$value
      }, numeric(1))
      expect_equal(
        exp(unname(moments[i + 2, ])), g,
        tolerance = 1e-8,
        info = paste(family$label, "t", ts[i])
      )
    }
  }
})

test_that("the analytic CDF matches a reference integral", {
  worst <- 0
  for (family in families) {
    for (pwindow in pwindows) {
      q <- sort(c(
        1e-6, 1e-3, 0.3 * pwindow, pwindow - 1e-3, pwindow, pwindow + 1e-3,
        2, 3, 6, 12, 25
      ))
      for (rho in rhos) {
        # Beyond the limit the numerical method is used, tested below
        inside <- q[abs(rho) * q <= series_limit]
        obj <- loglogistic_object(family, rho)
        expected <- exptilt_reference(
          inside, pwindow, rho, exptilt_cdf(family)
        )
        actual <- pcens_cdf(obj, inside, pwindow)
        worst <- max(worst, max_rel_diff(actual, expected))
        expect_lt(
          max_rel_diff(actual, expected), 1e-7,
          label = exptilt_label(family, pwindow, rho)
        )
      }
    }
  }
  expect_lt(worst, 1e-7)
})

test_that("the analytic CDF agrees with use_numeric = TRUE", {
  q <- c(0.05, 0.5, 1.5, 3, 6, 12, 20)
  for (family in families) {
    for (pwindow in c(1, 2, 7)) {
      for (rho in c(-1, -0.5, -1e-8, 1e-8, 0.3, 0.5, 1)) {
        obj <- loglogistic_object(family, rho)
        expect_equal(
          pcens_cdf(obj, q, pwindow),
          pcens_cdf(obj, q, pwindow, use_numeric = TRUE),
          tolerance = 1e-6,
          info = exptilt_label(family, pwindow, rho)
        )
      }
    }
  }
})

test_that("points beyond the limit use the numerical method", {
  obj <- loglogistic_object(families[[3]], 0.8)
  q <- c(2, 6, 12.4, 12.6, 30, 80)
  # The limit is 10 / 0.8 = 12.5
  result <- pcens_cdf(obj, q, 2)
  numeric <- pcens_cdf(obj, q, 2, use_numeric = TRUE)
  beyond <- q > 12.5
  expect_equal(result[beyond], numeric[beyond], tolerance = 1e-6)
  expect_equal(result[!beyond], numeric[!beyond], tolerance = 1e-6)
  expected <- exptilt_reference(q, 2, 0.8, exptilt_cdf(families[[3]]))
  expect_lt(max_rel_diff(result, expected), 1e-5)
  # A shape below 0.01 is numerical
  tiny <- list(pdist = pllogis_test, args = list(shape = 0.005, scale = 2))
  expect_false(.pcens_tilt_available(loglogistic_object(tiny, 0.3), -0.3))
  expect_identical(
    pcens_cdf(loglogistic_object(tiny, 0.3), c(1, 4), 2),
    pcens_cdf(loglogistic_object(tiny, 0.3), c(1, 4), 2, use_numeric = TRUE)
  )
  # A shape below 0.2 is numerical throughout
  slow <- exptilt_object(
    list(pdist = pllogis_test, args = list(shape = 0.15, scale = 2)), 0.3
  )
  expect_equal(
    pcens_cdf(slow, c(1, 4), 2),
    pcens_cdf(slow, c(1, 4), 2, use_numeric = TRUE),
    tolerance = 1e-6
  )
})

test_that("the CDF is accurate relative to its smaller tail near the limit", {
  for (case in conditioning_cases) {
    obj <- exptilt_object(
      list(pdist = pllogis_test, args = list(shape = case[1], scale = case[2])),
      case[3]
    )
    q <- c(2, 4, 6, 8, 9, 9.5, 9.9, 10) / abs(case[3])
    actual <- pcens_cdf(obj, q, case[4])
    for (i in seq_along(q)) {
      ref <- loglogistic_censored_reference(
        q[i], case[1], case[2], case[3], case[4]
      )
      # Below 1e-9 the reference cannot resolve the smaller tail
      if (min(ref) < 1e-9) next
      expect_lt(
        abs(actual[i] - ref[["cdf"]]) / min(ref), 1e-7,
        label = paste(
          "shape", case[1], "scale", case[2], "rho", case[3], "w", case[4],
          "q", q[i]
        )
      )
    }
  }
})

test_that("the numerical fallback resolves sharp shapes and tails", {
  # Shape 50 makes the delay CDF a near step at the scale, which the
  # default integration misses
  sharp <- function(scale, rho) {
    exptilt_object(
      list(pdist = pllogis_test, args = list(shape = 50, scale = scale)), rho
    )
  }
  ref <- loglogistic_censored_reference(5, 50, 1, 4, 7)
  expect_equal(
    pcens_cdf(sharp(1, 4), 5, 7), ref[["cdf"]],
    tolerance = 1e-8
  )
  ref <- loglogistic_censored_reference(10, 50, 5, -2, 7)
  expect_lt(
    abs((1 - pcens_cdf(sharp(5, -2), 10, 7)) - ref[["survival"]]) /
      ref[["survival"]],
    1e-6
  )
})

test_that("use_numeric = TRUE uses the default method for a tilt", {
  obj <- loglogistic_object(families[[2]], 0.5)
  expect_identical(
    pcens_cdf(obj, c(1, 4), 2, use_numeric = TRUE),
    pcens_cdf.default(obj, c(1, 4), 2)
  )
})

test_that("the analytic CDF is continuous in the tilt through zero", {
  q <- c(0.5, 2, 6, 15)
  for (family in families) {
    for (pwindow in c(1, 4)) {
      at_zero <- pcens_cdf(
        loglogistic_object(family, 0), q, pwindow
      )
      for (rho in c(-1e-6, 1e-6, -1e-3, 1e-3)) {
        expect_equal(
          pcens_cdf(loglogistic_object(family, rho), q, pwindow), at_zero,
          tolerance = 10 * abs(rho),
          info = exptilt_label(family, pwindow, rho)
        )
      }
    }
  }
})

test_that("the analytic CDF has no jump where the small tilt form ends", {
  q <- c(0.5, 2, 6, 15)
  for (family in families) {
    for (pwindow in c(0.5, 2)) {
      rho <- 1e-4 / pwindow
      below <- pcens_cdf(
        loglogistic_object(family, rho * (1 - 1e-9)), q, pwindow
      )
      above <- pcens_cdf(
        loglogistic_object(family, rho * (1 + 1e-9)), q, pwindow
      )
      expect_equal(below, above, tolerance = 1e-8)
    }
  }
})

test_that("the analytic CDF is accurate for q near zero", {
  for (family in families) {
    cdf <- exptilt_cdf(family)
    for (rho in c(-0.5, -1e-3, 1e-3, 0.5)) {
      for (pwindow in c(1, 5)) {
        q <- c(1e-10, 1e-7, 1e-4, 1e-2)
        expected <- exptilt_reference(q, pwindow, rho, cdf)
        actual <- pcens_cdf(loglogistic_object(family, rho), q, pwindow)
        expect_lt(
          max_rel_diff(actual, expected), 1e-7,
          label = exptilt_label(family, pwindow, rho)
        )
      }
    }
  }
})

test_that("the analytic CDF handles boundary values of q", {
  obj <- loglogistic_object(families[[3]], 0.3)
  expect_identical(pcens_cdf(obj, c(-Inf, -1, 0), 2), c(0, 0, 0))
  expect_identical(pcens_cdf(obj, Inf, 2), 1)
  expect_identical(pcens_cdf(obj, numeric(0), 2), numeric(0))
  expect_true(is.na(pcens_cdf(obj, NA_real_, 2)))
})

test_that("the analytic CDF is a CDF", {
  q <- seq(0, 40, by = 0.25)
  for (family in families) {
    for (rho in c(-0.3, 0.05, 0.3)) {
      result <- pcens_cdf(loglogistic_object(family, rho), q, 3)
      expect_true(all(result >= 0 & result <= 1))
      # Numerical fallback values are within 1e-6, so allow that in the
      # difference
      expect_true(all(diff(result) > -1e-6))
    }
  }
})

test_that("the pmf matches differences of the reference CDF", {
  for (family in families) {
    cdf <- exptilt_cdf(family)
    for (rho in c(-0.3, 1e-8, 0.3)) {
      obj <- loglogistic_object(family, rho)
      x <- 0:11
      inside <- abs(rho) * max(x + 1) <= series_limit
      if (!inside) next
      ref_cdf <- exptilt_reference(0:12, 3, rho, cdf)
      expect_equal(
        pcens_pmf(obj, x, pwindow = 3), diff(ref_cdf),
        tolerance = 1e-8
      )
    }
  }
})

test_that("the uniform CDF is the tilted CDF at zero tilt", {
  q <- c(0.5, 2, 6, 15)
  for (family in families) {
    expect_equal(
      pcens_cdf(uniform_object(family), q, 3),
      pcens_cdf(loglogistic_object(family, 0), q, 3),
      tolerance = 1e-9
    )
  }
})

test_that("the numerical fallback keeps the lower tail under a strong tilt", {
  # The primary mass sits near the window end, so the CDF is tiny for
  # q < pwindow even where the delay CDF at the window midpoint is above 0.5
  cases <- list(
    c(40, 0.25, 3, 20, 10.5),
    c(41.2, 0.94, 2.9, 27.6, 17.7),
    c(1.66, 2.06, 2.52, 27.5, 16.6),
    c(6.6, 0.25, 2.24, 26, 16.2)
  )
  for (case in cases) {
    obj <- exptilt_object(
      list(pdist = pllogis_test, args = list(shape = case[1], scale = case[2])),
      case[3]
    )
    ref <- loglogistic_censored_reference(
      case[5], case[1], case[2], case[3], case[4]
    )
    # The fallback is used as rho * q is beyond the series limit
    expect_gt(case[3] * case[5], series_limit)
    expect_lt(
      abs(pcens_cdf(obj, case[5], case[4]) - ref[["cdf"]]) / min(ref), 1e-6,
      label = paste(case, collapse = " ")
    )
  }
})

test_that("the small tilt form is accurate where its difference cancels", {
  cases <- list(
    c(115.5, 16.9, 0.00179, 1e-5, 19.8),
    c(0.5, 5, 1, 5e-5, 1e5),
    c(8.72, 0.2153, 0.0355, 0.00125, 1.72),
    c(0.01, 1, 1, -1e-5, 1e8)
  )
  for (case in cases) {
    obj <- exptilt_object(
      list(pdist = pllogis_test, args = list(shape = case[1], scale = case[2])),
      case[4]
    )
    # The small tilt form applies as |rho| w is below 1e-4
    expect_lt(abs(case[4]) * case[3], 1e-4)
    ref <- loglogistic_censored_reference(
      case[5], case[1], case[2], case[4], case[3]
    )
    actual <- pcens_cdf(obj, case[5], case[3])
    expect_lt(
      abs(actual - ref[["cdf"]]) / max(min(ref), 1e-9), 1e-6,
      label = paste(case, collapse = " ")
    )
  }
})

test_that("the tilted analytic CDFs match rprimarycensored samples", {
  withr::local_seed(2)
  n <- 20000
  pwindow <- 2
  for (family in families[c(1, 3, 5)]) {
    for (rho in c(-0.5, 0.4)) {
      samples <- loglogistic_samples(n, family, pwindow, rho)
      qs <- unname(stats::quantile(samples, c(0.05, 0.25, 0.5, 0.75, 0.95)))
      expect_lt(
        max(abs(
          vapply(qs, function(q) mean(samples <= q), numeric(1)) -
            pcens_cdf(loglogistic_object(family, rho), qs, pwindow)
        )),
        0.015,
        label = exptilt_label(family, pwindow, rho)
      )
    }
  }
})
