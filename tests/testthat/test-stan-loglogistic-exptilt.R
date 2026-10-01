skip_on_cran()

# Log-logistic delay (dist_id 31, params = [scale, shape]) with a tilted
# primary (primary_id 2, params = r).

ll_limit <- 10


test_that("dispatch checks include the tilted log-logistic solutions", {
  expect_identical(check_for_exptilt(31L, 2L), 1L)
  expect_identical(check_for_analytical_vectorized(31L, 2L, 3), 1L)
  expect_identical(check_for_analytical_vectorized(31L, 2L, 0.5), 0L)
  # The transform exists for any tilt, the series region depends on the point
  for (xi in c(-5, -0.5, 0, 0.5, 5)) {
    expect_identical(check_for_tilt_transform(31L, xi, c(3, 2)), 1L)
  }
  for (rho in c(-0.7, 0.3, 0.7)) {
    expect_identical(
      check_for_analytical_params(31L, c(3, 2), 2L, rho), 1L
    )
  }
  # A tiny shape has no solution at all
  expect_identical(check_for_tilt_transform(31L, -0.5, c(3, 0.005)), 0L)
  expect_identical(check_for_tilt_transform(31L, -0.5, c(3, 0.01)), 1L)
  expect_identical(
    check_for_analytical_params(31L, c(3, 0.005), 2L, 0.3), 0L
  )
})

test_that("Stan log-logistic tilt transforms match the R transforms", {
  ts <- c(1e-3, 0.4, 1, 3.5, 8, 12)
  for (case in ll_cases) {
    obj <- exptilt_object(ll_case_family(case), 0.1)
    for (xi in c(0, -1, -0.25, -0.05, 0.05, 0.25, 0.6)) {
      t <- ts[abs(xi) * ts <= ll_limit]
      pairs <- vapply(
        t, log_tilt_transform_pair, numeric(2), 31L, xi, case$params
      )
      lower <- pairs[1, ]
      upper <- pairs[2, ]
      expect_equal(
        lower, .pcens_tilt_transform(obj, t, xi),
        tolerance = 1e-10,
        info = ll_case_label(case, xi = xi)
      )
      if (xi == 0) {
        expect_equal(
          upper, .pcens_tilt_transform(obj, t, xi, upper = TRUE),
          tolerance = 1e-12
        )
      } else {
        expect_true(all(is.nan(upper)))
      }
    }
  }
})

test_that("Stan log-logistic transforms are nan beyond the series limit", {
  params <- c(3, 2)
  expect_true(all(is.finite(log_tilt_transform_pair(10, 31L, 1, params)[1])))
  expect_true(is.nan(log_tilt_transform_pair(10.01, 31L, 1, params)[1]))
  expect_true(is.nan(log_tilt_transform_pair(10.01, 31L, -1, params)[1]))
  expect_true(is.finite(log_tilt_transform_pair(1e6, 31L, 0, params)[1]))
  expect_identical(log_tilt_transform_pair(-5, 31L, -1, params)[1], -Inf)
  # A shape below 0.2 has no series
  expect_true(is.nan(log_tilt_transform_pair(1, 31L, 1, c(3, 0.19))[1]))
  expect_true(is.finite(log_tilt_transform_pair(1, 31L, 1, c(3, 0.2))[1]))
  expect_identical(loglogistic_tilt_in_range(10, 1, 2), 1L)
  expect_identical(loglogistic_tilt_in_range(10.01, -1, 2), 0L)
  expect_identical(loglogistic_tilt_in_range(1e6, 0, 2), 1L)
  expect_identical(loglogistic_tilt_in_range(1, 1, 0.19), 0L)
})

test_that("Stan log-logistic transforms are 0 or nan below the support", {
  params <- c(3, 2)
  expect_identical(log_tilt_transform_pair(0, 31L, -0.3, params)[1], -Inf)
  expect_identical(log_tilt_transform_pair(-2, 31L, 0.3, params)[1], -Inf)
  expect_identical(log_tilt_transform_pair(-2, 31L, 0, params)[2], 0)
  expect_true(is.nan(log_tilt_transform_pair(-2, 31L, -0.3, params)[2]))
})

test_that("Stan log-logistic tilt moments match the R moments", {
  ts <- c(1e-3, 0.4, 1, 3.5, 12, 60, 1e4)
  for (case in ll_cases) {
    obj <- exptilt_object(ll_case_family(case), 0.1)
    expected <- .pcens_tilt_moments(obj, ts)
    actual <- t(vapply(
      ts, primarycensored_tilt_moments, numeric(3), 31L, case$params
    ))
    expect_equal(
      unname(actual), unname(expected),
      tolerance = 1e-9,
      info = ll_case_label(case)
    )
  }
  expect_identical(
    primarycensored_tilt_moments(-1, 31L, c(3, 2)), rep(-Inf, 3)
  )
})

test_that("primarycensored_exptilt_lcdf matches the R implementation", {
  d <- c(1e-4, 0.3, 1, 2.5, 6, 15, 30)
  for (case in ll_cases) {
    for (pwindow in c(0.5, 2, 7)) {
      for (rho in c(-0.2, -1e-6, 1e-6, 0.3)) {
        obj <- exptilt_object(ll_case_family(case), rho)
        inside <- d[abs(rho) * d <= ll_limit]
        expect_equal(
          exp(vapply(
            inside, primarycensored_exptilt_lcdf, numeric(1),
            31L, case$params, pwindow, rho
          )),
          pcens_cdf(obj, inside, pwindow),
          tolerance = 1e-9,
          info = ll_case_label(case, pwindow = pwindow, r = rho)
        )
      }
    }
  }
})

test_that("points beyond the series limit use a tight numerical CDF", {
  params <- c(5, 2)
  rho <- 0.8
  ode_model <- ll_ode_model()
  for (d in c(13, 20, 50)) {
    numeric_lcdf <- loglogistic_numeric_lcdf(d, params, 2, 2L, rho)
    expect_identical(
      primarycensored_exptilt_lcdf(d, 31L, params, 2, rho), numeric_lcdf
    )
    expect_identical(
      primarycensored_lcdf(d, 31L, params, 2, 0, Inf, 2L, rho), numeric_lcdf
    )
    expect_equal(
      exp(numeric_lcdf), ll_ode_cdf(ode_model, params, d, 2, 2L, rho),
      tolerance = 1e-4
    )
    ref <- loglogistic_censored_reference(d, 2, 5, rho, 2)
    smaller <- if (ref[["cdf"]] > 0.5) {
      -expm1(numeric_lcdf)
    } else {
      exp(numeric_lcdf)
    }
    expect_lt(abs(smaller - min(ref)) / min(ref), 1e-8)
  }
  # A small shape is always numerical, and the tight CDF is accurate
  slow <- c(3, 0.15)
  ref <- loglogistic_censored_reference(2, 0.15, 3, 0.3, 2)
  expect_equal(
    exp(primarycensored_exptilt_lcdf(2, 31L, slow, 2, 0.3)),
    ref[["cdf"]],
    tolerance = 1e-6
  )
  # A tiny shape uses the ODE for both primaries
  tiny <- c(3, 0.005)
  expect_identical(
    primarycensored_lcdf(2, 31L, tiny, 2, 0, Inf, 2L, 0.3),
    log(ll_ode_cdf(ode_model, tiny, 2, 2, 2L, 0.3))
  )
  expect_identical(
    primarycensored_lcdf(2, 31L, tiny, 2, 0, Inf, 1L, numeric(0)),
    log(ll_ode_cdf(ode_model, tiny, 2, 2, 1L, numeric(0)))
  )
  # The small tilt forms do not need the series
  expect_false(
    identical(
      primarycensored_exptilt_lcdf(30, 31L, params, 2, 1e-5),
      log(ll_ode_cdf(ode_model, params, 30, 2, 2L, 1e-5))
    )
  )
})

test_that("the numerical fallback tail rule uses the same primary mean", {
  for (rho in c(-3, -0.2, 0, 1e-4, 0.5, 4)) {
    for (pwindow in c(0.5, 2)) {
      expect_equal(
        .loglogistic_primary_mean(rho, pwindow),
        loglogistic_primary_mean(
          if (rho == 0) 1L else 2L, if (rho == 0) numeric(0) else rho,
          pwindow
        ),
        tolerance = 1e-13
      )
    }
  }
})

test_that("the numerical CDF resolves a deep lower tail of a sharp shape", {
  # Delays either side of the series limit |rho| d = 10
  params <- c(20, 40)
  for (d in c(4, 5, 6, 7, 8)) {
    ref <- loglogistic_censored_reference(d, 40, 20, 2.5, 1)
    expect_lt(
      abs(expm1(primarycensored_lcdf(d, 31L, params, 1, 0, Inf, 2L, 2.5) -
        log(ref[["cdf"]]))),
      1e-7,
      label = paste("d", d)
    )
  }
})

test_that("the tilted CDF is accurate relative to its smaller tail", {
  for (case in conditioning_cases) {
    params <- c(case[2], case[1])
    q <- c(2, 4, 6, 8, 9, 9.5, 9.9, 10) / abs(case[3])
    for (i in seq_along(q)) {
      ref <- loglogistic_censored_reference(
        q[i], case[1], case[2], case[3], case[4]
      )
      # Below 1e-6 the tolerance of the numerical CDF cannot resolve the
      # smaller tail
      if (min(ref) < 1e-6) next
      actual <- exp(primarycensored_exptilt_lcdf(
        q[i], 31L, params, case[4], case[3]
      ))
      expect_lt(
        abs(actual - ref[["cdf"]]) / min(ref), 1e-7,
        label = paste(
          "shape", case[1], "scale", case[2], "rho", case[3], "w", case[4],
          "q", q[i]
        )
      )
    }
  }
})

test_that("a sharp shape beyond the series limit gives a valid PMF", {
  # The numerical CDF can exceed 1 by its tolerance and be out of order
  params <- c(6, 30)
  for (d in 12:25) {
    expect_lte(primarycensored_lcdf(d, 31L, params, 1, 0, Inf, 2L, 1), 0)
    expect_lte(primarycensored_exptilt_lcdf(d, 31L, params, 1, 1), 0)
  }
  for (big_d in c(Inf, 41)) {
    lpmf <- primarycensored_sone_lpmf_vectorized(
      40L, 0, big_d, 31L, params, 1, 2L, 1
    )
    expect_false(anyNA(lpmf))
    expect_true(all(lpmf <= 0))
    expect_equal(sum(exp(lpmf)), if (is.infinite(big_d)) {
      exp(primarycensored_lcdf(41, 31L, params, 1, 0, Inf, 2L, 1))
    } else {
      1
    }, tolerance = 1e-6)
  }
})

test_that("primarycensored_lcdf and primarycensored_cdf use the tilted
  solution and agree with the ODE path", {
  d <- c(0.2, 1, 2.5, 6, 15)
  ode_model <- ll_ode_model()
  for (case in ll_cases) {
    cdf <- ll_case_cdf(case)
    for (pwindow in c(1, 3)) {
      for (rho in c(-1, -0.5, -1e-8, 1e-8, 0.5)) {
        inside <- d[abs(rho) * d <= ll_limit]
        info <- ll_case_label(case, pwindow = pwindow, r = rho)
        expect_identical(
          check_for_analytical_params(31L, case$params, 2L, rho), 1L
        )
        expected <- exptilt_reference(inside, pwindow, rho, cdf)
        lcdf <- vapply(
          inside, primarycensored_lcdf, numeric(1),
          31L, case$params, pwindow, 0, Inf, 2L, rho
        )
        expect_lt(max_rel_diff(exp(lcdf), expected), 1e-7, label = info)
        plain <- vapply(
          inside, primarycensored_cdf, numeric(1),
          31L, case$params, pwindow, 0, Inf, 2L, rho
        )
        expect_lt(max_rel_diff(plain, expected), 1e-7, label = info)
        # The ODE path has absolute and relative tolerances of 1e-6, and
        # about 1e-4 for a shape below 1 where the density is singular at 0
        ode <- ll_ode_cdf(ode_model, case$params, inside, pwindow, 2L, rho)
        expect_lt(max(abs(plain - ode)), 1e-4, label = info)
      }
    }
  }
})

test_that("the vectorised tilted CDF matches the per delay CDF", {
  for (case in ll_cases) {
    for (pwindow in c(1, 3, 7)) {
      for (rho in c(-0.4, -1e-6, 1e-6, 0.05, 0.3, 1)) {
        n <- 25L
        vec <- primarycensored_analytical_lcdf_vectorized(
          1L, n, 31L, case$params, pwindow, 2L, rho
        )
        single <- vapply(
          1:n, primarycensored_exptilt_lcdf, numeric(1),
          31L, case$params, pwindow, rho
        )
        expect_equal(
          vec, single,
          tolerance = 1e-9,
          info = ll_case_label(case, pwindow = pwindow, r = rho)
        )
      }
    }
  }
})

test_that("primarycensored_lcdf_vectorized shares terms for the tilted
  primary", {
  params <- c(5, 2)
  vec <- primarycensored_lcdf_vectorized(1L, 20L, 31L, params, 3, 2L, 0.3)
  single <- vapply(
    1:20, primarycensored_lcdf, numeric(1),
    31L, params, 3, 0, Inf, 2L, 0.3
  )
  expect_equal(vec, single, tolerance = 1e-9)
})

test_that("the tilted vectorised PMF matches the per delay PMF with
  truncation", {
  params <- c(5, 2)
  for (rho in c(0.3, -0.3)) {
    for (bounds in list(c(0, Inf), c(0, 25), c(2, 25))) {
      vec <- primarycensored_sone_lpmf_vectorized(
        15L, bounds[1], bounds[2], 31L, params, 3, 2L, rho
      )
      single <- vapply(1:15, function(d) {
        primarycensored_lpmf(
          d - 1L, 31L, params, 3, d, bounds[1], bounds[2], 2L, rho
        )
      }, numeric(1))
      expect_equal(vec[1:15], single, tolerance = 1e-8)
    }
  }
})

test_that("log-logistic tilted log CDFs have finite gradients matching
  finite differences", {
  model <- ll_gradient_model()
  # Direct, small window and small delay forms; d below and above pwindow;
  # the tilted delay, the uniform primary and a point beyond the series limit
  points <- list(
    list(d = 0.3, pwindow = 2, rho = 0.4),
    list(d = 1, pwindow = 1, rho = -0.2),
    list(d = 2.5, pwindow = 2, rho = 0.4),
    list(d = 6, pwindow = 3, rho = -0.15),
    list(d = 14, pwindow = 3, rho = 0.5),
    list(d = 2.5, pwindow = 2, rho = 1e-6),
    list(d = 12, pwindow = 7, rho = -1e-5),
    list(d = 0.0001, pwindow = 2, rho = 0.4),
    list(d = 4, pwindow = 2, rho = 0)
  )
  for (case in ll_cases) {
    for (point in points) {
      expect_ll_gradient(model, case, point, 2L)
    }
  }
})

test_that("the tilted vectorised log PMF has finite gradients matching
  finite differences", {
  model <- ll_gradient_model()
  points <- list(
    list(d = 12, pwindow = 3, rho = 0.4, primary_id = 2L),
    list(d = 12, pwindow = 3, rho = -0.15, primary_id = 2L),
    list(d = 6, pwindow = 2, rho = 1e-6, primary_id = 2L)
  )
  for (case in ll_cases[c(2, 3, 4)]) {
    for (point in points) {
      expect_ll_gradient(model, case, point, point$primary_id, TRUE)
    }
  }
})

test_that("the vectorised log PMF is finite for shapes below the series
  minimum where (t / scale)^shape is beyond the tail start", {
  model <- ll_gradient_model()
  # The ODE path is used for a shape below 0.2, so no tail series pivots
  # are needed, and these overflow for large powers n / shape
  cases <- list(
    list(params = c(1e-5, 0.08)),
    list(params = c(1e-9, 0.05)),
    list(params = c(1e-6, 0.15))
  )
  points <- list(
    list(d = 30, pwindow = 1, rho = 0.3, primary_id = 2L),
    list(d = 30, pwindow = 1, rho = -0.3, primary_id = 2L),
    list(d = 31, pwindow = 1, rho = 0.3, primary_id = 2L),
    list(d = 30, pwindow = 2, rho = 0, primary_id = 1L)
  )
  for (case in cases) {
    for (point in points) {
      for (vectorised in c(FALSE, TRUE)) {
        expect_ll_gradient(
          model, case, point, point$primary_id, vectorised, tolerance = 1e-3
        )
      }
    }
  }
})

test_that("the Stan tilted CDF resolves the upper tail", {
  for (case in upper_tail_cases) {
    ref <- loglogistic_censored_reference(
      case[5], case[1], case[2], case[4], case[3]
    )
    lcdf <- primarycensored_lcdf(
      case[5], 31L, c(case[2], case[1]), case[3], 0, Inf, 2L, case[4]
    )
    # The survival from the CDF is resolved to 1e-16 in absolute terms
    expect_lt(
      abs(-expm1(lcdf) - ref[["survival"]]) / ref[["survival"]], 1e-6,
      label = paste(case, collapse = " ")
    )
  }
})

test_that("the Stan tilted PMF in the upper tail is accurate and valid", {
  for (case in upper_tail_pmf_cases) {
    params <- c(case[2], case[1])
    expected <- loglogistic_pmf_at(case)
    lpmf <- primarycensored_lpmf(
      as.integer(case[5]), 31L, params, case[3], case[5] + 1, 0, Inf, 2L,
      case[4]
    )
    expect_false(is.nan(lpmf), label = paste(case, collapse = " "))
    expect_lt(
      abs(exp(lpmf) - expected) / expected, 1e-4,
      label = paste(case, collapse = " ")
    )
    # The CDF does not decrease in the delay
    lcdf <- vapply(
      case[5] + 0:3, primarycensored_lcdf, numeric(1),
      31L, params, case[3], 0, Inf, 2L, case[4]
    )
    expect_true(all(diff(lcdf) >= 0), label = paste(case, collapse = " "))
  }
})

test_that("the Stan tilted numerical CDF does not reach the step limit at
  a large delay", {
  case <- c(1.57, 1.2, 0.31, -2.16, 3e4)
  expect_no_error(
    primarycensored_lcdf(
      case[5], 31L, c(case[2], case[1]), case[3], 0, Inf, 2L, case[4]
    )
  )
})

test_that("the Stan small tilt form is accurate where its difference cancels", {
  cases <- list(
    c(115.5, 16.9, 0.00179, 1e-5, 19.8),
    c(0.5, 5, 1, 5e-5, 1e5),
    c(8.72, 0.2153, 0.0355, 0.00125, 1.72),
    c(0.01, 1, 1, -1e-5, 1e8)
  )
  for (case in cases) {
    expect_lt(abs(case[4]) * case[3], 1e-2)
    ref <- loglogistic_censored_reference(
      case[5], case[1], case[2], case[4], case[3]
    )
    lcdf <- primarycensored_lcdf(
      case[5], 31L, c(case[2], case[1]), case[3], 0, Inf, 2L, case[4]
    )
    expect_lt(
      abs(exp(lcdf) - ref[["cdf"]]) / max(min(ref), 1e-9), 1e-6,
      label = paste(case, collapse = " ")
    )
    expect_lt(
      abs(-expm1(lcdf) - ref[["survival"]]) / max(ref[["survival"]], 1e-9),
      1e-6,
      label = paste(case, collapse = " ")
    )
  }
})

test_that("the tilted fallback has gradients matching differences", {
  model <- ll_gradient_model()
  # Tilted points beyond the series limit, in the upper tail and in the
  # lower tail
  tilted <- list(
    list(case = ll_cases[[3]], d = 30, pwindow = 2, rho = 0.5),
    list(case = ll_cases[[2]], d = 60, pwindow = 3, rho = -0.4),
    list(case = ll_cases[[4]], d = 4, pwindow = 10, rho = 3),
    list(case = ll_cases[[1]], d = 500, pwindow = 1, rho = 0.2),
    list(case = ll_cases[[3]], d = 3, pwindow = 10, rho = 2),
    list(case = list(params = c(20, 40)), d = 5, pwindow = 1, rho = 2.5)
  )
  for (point in tilted) {
    expect_ll_gradient(model, point$case, point, 2L, tolerance = 1e-3)
  }
  expect_ll_gradient(
    model, list(params = c(20, 40)), list(d = 15, pwindow = 1, rho = 2.5),
    2L, vectorised = TRUE, tolerance = 1e-2
  )
})

test_that("Stan tilted CDFs match rprimarycensored samples", {
  withr::local_seed(3)
  n <- 20000
  pwindow <- 2
  for (case in ll_cases[c(1, 3, 5)]) {
    for (rho in c(-0.5, 0.4)) {
      samples <- loglogistic_samples(n, ll_case_family(case), pwindow, rho)
      qs <- unname(stats::quantile(samples, c(0.05, 0.25, 0.5, 0.75, 0.95)))
      stan_cdf <- exp(vapply(
        qs, primarycensored_lcdf, numeric(1), 31L, case$params, pwindow, 0,
        Inf, 2L, rho
      ))
      expect_lt(
        max(abs(vapply(qs, function(q) mean(samples <= q), numeric(1)) -
          stan_cdf)),
        0.015,
        label = ll_case_label(case, pwindow = pwindow, r = rho)
      )
    }
  }
})
