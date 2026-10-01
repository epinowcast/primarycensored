skip_on_cran()

# Stan solutions for exponential, gamma and normal delays with a truncated
# logistic primary (primary_id 3 with primary_params = c(location, scale)).
# These tests check the Stan series against the R implementation and a
# reference integral, the dispatch and fallback to the ODE path, the shared
# endpoint vectorised form, and gradients.

# Delays below, inside and above the window. Stan returns -Inf for probabilities
# below the smallest double, so values that underflow are not compared.
tlogis_stan_delays <- function(case, pwindow) {
  scale <- if (case$params[length(case$params)] > 100) 1e-3 else 1
  sort(c(
    1e-6, 1e-3, 0.3 * pwindow, pwindow - 1e-3, pwindow, pwindow + 1e-3,
    2, 3, 6, 12, 25, if (case$dist_id == 18L) c(-10, -3, -0.5)
  ) * scale)
}

# Whether a case uses the analytic path for a location, scale and window
tlogis_case_analytic <- function(case, location, scale, pwindow) {
  as.logical(check_for_tlogis_params(
    case$dist_id, case$params, c(location, scale), pwindow
  ))
}

test_that("check_for_tlogis is structural and for primary_id 3", {
  for (dist_id in c(2L, 4L, 18L)) {
    expect_identical(check_for_tlogis(dist_id, 3L), 1L)
    expect_identical(check_for_tlogis(dist_id, 2L), 0L)
    expect_identical(check_for_tlogis(dist_id, 1L), 0L)
    expect_identical(check_for_analytical(dist_id, 3L), 1L)
  }
  for (dist_id in c(1L, 3L, 5L, 9L)) {
    expect_identical(check_for_tlogis(dist_id, 3L), 0L)
    expect_identical(check_for_analytical(dist_id, 3L), 0L)
  }
  # The non-parametric delays are analytic for this primary
  for (dist_id in c(26L, 27L, 28L)) {
    expect_identical(check_for_analytical(dist_id, 3L), 1L)
  }
  expect_identical(check_for_analytical(2L, 1L), 1L)
  expect_identical(check_for_analytical(4L, 2L), 1L)
  expect_identical(check_for_analytical(2L, 2L), 1L)
})

test_that("the series rule matches R", {
  ranges <- list(
    c(0, 1), c(0.2, 1), c(0.9, 1), c(0, 0.6), c(0.1, 0.3), c(1e-3, 1e-2)
  )
  for (tol in c(1e-8, 1e-12)) {
    for (range in ranges) {
      r_terms <- .tlogis_series_terms(
        log(range[1]), log(range[2]), log(tol)
      )
      stan_terms <- tlogis_series_terms(log(range[1]), log(range[2]), log(tol))
      expect_identical(
        stan_terms, as.integer(unname(r_terms)),
        info = paste(tol, toString(range))
      )
    }
  }
  # No rule gives the marker -1
  expect_identical(
    tlogis_series_terms(log(1e-300), 0, log(1e-300)), c(-1L, -1L)
  )
  expect_identical(tlogis_series_terms(-Inf, 0, -Inf), c(-1L, -1L))
})

test_that("the weights match R", {
  for (args in list(c(3L, 5L), c(0L, 7L), c(9L, 17L), c(4L, 0L), c(0L, 1L))) {
    expect_equal(
      tlogis_weights(args[1], args[2]),
      .tlogis_weights(args[1], args[2]),
      tolerance = 1e-14
    )
  }
})

test_that("the plan matches R", {
  cases <- list(
    c(-0.5, 0.2, 1), c(0, 0.1, 1), c(0.5, 0.2, 1), c(1, 1, 2),
    c(1.5, 0.1, 1), c(40, 5, 2), c(-6, 0.02, 0.5)
  )
  for (cs in cases) {
    location <- cs[1]
    scale <- cs[2]
    pwindow <- cs[3]
    log_mass <- .tlogis_log_diff(0, pwindow, location, scale)
    plan <- tlogis_plan(location, scale, pwindow, log_mass)
    obj <- tlogis_object(
      list(pdist = pnorm, args = list(mean = 3, sd = 2)), location, scale
    )
    r_plan <- .tlogis_plan(obj, pwindow)
    r_form <- if (is.null(r_plan$pos)) {
      0L
    } else {
      c(A = 1L, P = 2L)[[r_plan$pos$form]]
    }
    expect_identical(plan[1], 1L, info = toString(cs))
    expect_identical(plan[2], r_form, info = toString(cs))
    if (!is.null(r_plan$pos)) {
      expect_identical(plan[3:4], as.integer(c(r_plan$pos$n0, r_plan$pos$M)))
    }
    if (!is.null(r_plan$neg)) {
      expect_identical(plan[5:6], as.integer(c(r_plan$neg$n0, r_plan$neg$M)))
    } else {
      expect_identical(plan[5:6], c(0L, 0L))
    }
  }
})

test_that("check_for_tlogis_params needs the tilts of the series", {
  # A location after the window needs only negative tilts
  expect_identical(
    check_for_tlogis_params(2L, c(2.5, 0.4), c(3, 0.5), 2), 1L
  )
  expect_identical(
    check_for_tlogis_params(4L, 0.3, c(3, 0.5), 2), 1L
  )
  # Positive tilts of the series reach the rate of these delays
  expect_identical(
    check_for_tlogis_params(2L, c(2.5, 0.4), c(0.5, 0.2), 2), 0L
  )
  expect_identical(
    check_for_tlogis_params(4L, 0.3, c(-0.5, 1), 2), 0L
  )
  expect_identical(
    check_for_tlogis_params(2L, c(2.5, 1000), c(0.5, 1), 2), 1L
  )
  # The normal delay has no restriction
  expect_identical(
    check_for_tlogis_params(18L, c(3, 2), c(0.5, 0.05), 2), 1L
  )
  # Invalid windows, scales and locations
  expect_identical(check_for_tlogis_params(18L, c(3, 2), c(0.5, 0.2), 0), 0L)
  expect_identical(check_for_tlogis_params(18L, c(3, 2), c(0.5, 0.2), -1), 0L)
  expect_identical(
    check_for_tlogis_params(18L, c(3, 2), c(0.5, 0.2), Inf), 0L
  )
  expect_identical(check_for_tlogis_params(18L, c(3, 2), c(0.5, 0), 2), 0L)
  expect_identical(check_for_tlogis_params(18L, c(3, 2), c(0.5, -1), 2), 0L)
  expect_identical(
    check_for_tlogis_params(18L, c(3, 2), c(Inf, 0.2), 2), 0L
  )
  # check_for_analytical_params adds it to the structural checks
  expect_identical(
    check_for_analytical_params(2L, c(2.5, 0.4), 3L, c(3, 0.5), 2), 1L
  )
  expect_identical(
    check_for_analytical_params(2L, c(2.5, 0.4), 3L, c(0.5, 0.2), 2), 0L
  )
  expect_identical(
    check_for_analytical_params(3L, c(2.5, 0.4), 3L, c(3, 0.5), 2), 0L
  )
  expect_identical(
    check_for_analytical_params(2L, c(2.5, 0.4), 1L, numeric(0), 2), 1L
  )
  expect_identical(
    check_for_analytical_params(2L, c(2.5, 0.4), 2L, -0.4, 2), 0L
  )
  expect_identical(
    check_for_analytical_params(2L, c(2.5, 0.4), 2L, 0.5, 2), 1L
  )
})

test_that("primarycensored_tlogis_lcdf matches a reference integral", {
  locations <- c(-0.5, 0, 0.5, 1, 1.5)
  scales <- c(0.1, 0.2, 1)
  for (case in tlogis_stan_cases) {
    positive <- case$dist_id != 18L
    cdf <- tlogis_case_cdf(case)
    for (pwindow in c(1, 2)) {
      d <- tlogis_stan_delays(case, pwindow)
      for (location in locations) {
        for (scale in scales) {
          if (!tlogis_case_analytic(case, location, scale, pwindow)) {
            next
          }
          expected <- tlogis_reference(
            d, pwindow, location, scale, cdf, positive
          )
          actual <- exp(vapply(
            d, primarycensored_tlogis_lcdf, numeric(1),
            case$dist_id, case$params, pwindow, location, scale
          ))
          expect_lt(
            max_rel_diff(actual, expected), 1e-8,
            label = tlogis_case_label(
              case,
              pwindow = pwindow, location = location, scale = scale
            )
          )
        }
      }
    }
  }
})

test_that("primarycensored_tlogis_lcdf matches the R implementation", {
  for (case in tlogis_stan_cases) {
    for (pwindow in c(0.5, 2, 7)) {
      d <- tlogis_stan_delays(case, pwindow)
      for (location in c(-2, 0, 0.3, 1, pwindow, pwindow + 1.5)) {
        for (scale in c(0.15, 0.7, 4)) {
          if (!tlogis_case_analytic(case, location, scale, pwindow)) {
            next
          }
          obj <- tlogis_case_obj(case, location, scale)
          expect_equal(
            exp(vapply(
              d, primarycensored_tlogis_lcdf, numeric(1),
              case$dist_id, case$params, pwindow, location, scale
            )),
            pcens_cdf(obj, d, pwindow),
            tolerance = 1e-9,
            info = tlogis_case_label(
              case,
              pwindow = pwindow, location = location, scale = scale
            )
          )
        }
      }
    }
  }
})

test_that("primarycensored_tlogis_lcdf matches Stan random draws", {
  set.seed(20260930)
  n <- 5000
  pwindow <- 2
  probs <- seq(0.1, 0.9, by = 0.2)
  # A large rate keeps the analytic path for a location before the window
  fast_exponential <- list(
    dist_id = 4L, params = 200, pdist = pexp, args = list(rate = 200)
  )
  cases <- c(
    tlogis_stan_cases[c(1, 4, 6, 7)], list(fast_exponential)
  )
  n_analytic <- 0L
  for (case in cases) {
    rdist <- switch(as.character(case$dist_id),
      "4" = function(n) rexp(n, case$args$rate),
      "2" = function(n) rgamma(n, case$args$shape, case$args$rate),
      "18" = function(n) rnorm(n, case$args$mean, case$args$sd)
    )
    for (window in list(c(-0.5, 0.3), c(1, 0.4), c(4, 0.5))) {
      if (!tlogis_case_analytic(case, window[1], window[2], pwindow)) {
        next
      }
      n_analytic <- n_analytic + 1L
      primary <- vapply(
        seq_len(n),
        function(i) tlogis_rng(0, pwindow, window[1], window[2]),
        numeric(1)
      )
      q <- unname(quantile(primary + rdist(n), probs))
      cdf <- exp(vapply(
        q, primarycensored_tlogis_lcdf, numeric(1),
        case$dist_id, case$params, pwindow, window[1], window[2]
      ))
      expect_lt(
        max(abs(cdf - probs)), 0.03,
        label = tlogis_case_label(
          case,
          pwindow = pwindow, location = window[1], scale = window[2]
        )
      )
    }
  }
  # Every location for the large rates, and the location after the window
  # for the other delays
  expect_identical(n_analytic, 11L)
})

test_that("the small delay form is used and continuous", {
  # q / scale below 1e-4 for a delay on the non-negative reals
  for (case in tlogis_stan_cases[c(2, 4, 6)]) {
    obj <- tlogis_case_obj(case, 3, 2)
    d <- c(1e-9, 1e-7, 1e-5, 1.9999e-4, 2.0001e-4, 1e-3)
    expect_equal(
      exp(vapply(
        d, primarycensored_tlogis_lcdf, numeric(1),
        case$dist_id, case$params, 2, 3, 2
      )),
      pcens_cdf(obj, d, 2),
      tolerance = 1e-8,
      info = tlogis_case_label(case)
    )
  }
})

test_that("primarycensored_lcdf and primarycensored_cdf use the analytical
  solution and agree with the ODE path", {
  d <- c(0.2, 1, 2.5, 6, 15)
  for (case in tlogis_stan_cases) {
    lower <- tlogis_case_lower(case)
    cdf <- tlogis_case_cdf(case)
    for (pwindow in c(1, 3)) {
      for (window in list(c(-0.5, 0.5), c(0.5, 0.2), c(pwindow + 1, 1))) {
        location <- window[1]
        scale <- window[2]
        if (!tlogis_case_analytic(case, location, scale, pwindow)) {
          next
        }
        dd <- d * if (case$params[length(case$params)] > 100) 1e-3 else 1
        info <- tlogis_case_label(
          case,
          pwindow = pwindow, location = location, scale = scale
        )
        expect_identical(
          check_for_analytical_params(
            case$dist_id, case$params, 3L, window, pwindow
          ), 1L
        )
        expected <- tlogis_reference(
          dd, pwindow, location, scale, cdf, case$dist_id != 18L
        )
        lcdf <- vapply(
          dd, primarycensored_lcdf, numeric(1),
          case$dist_id, case$params, pwindow, lower, Inf, 3L, window
        )
        expect_lt(max_rel_diff(exp(lcdf), expected), 1e-8, label = info)
        plain <- vapply(
          dd, primarycensored_cdf, numeric(1),
          case$dist_id, case$params, pwindow, lower, Inf, 3L, window
        )
        expect_lt(max_rel_diff(plain, expected), 1e-8, label = info)
        # The ODE path is solved to a relative tolerance of 1e-9 and an
        # absolute tolerance of 1e-10, and less for a shape below 1 where
        # the density is singular at 0
        ode <- vapply(
          dd, primarycensored_numeric_cdf, numeric(1),
          case$dist_id, case$params, pwindow, 3L, window
        )
        expect_lt(max(abs(plain - ode)), 1e-6, label = info)
      }
    }
  }
})

test_that("inadmissible series use the ODE path", {
  cases <- list(
    list(dist_id = 4L, params = 0.3, window = c(0.5, 0.2)),
    list(dist_id = 4L, params = 0.3, window = c(-1, 1)),
    list(dist_id = 2L, params = c(2.5, 0.4), window = c(0.5, 0.2)),
    list(dist_id = 2L, params = c(2.5, 0.4), window = c(1, 1))
  )
  for (case in cases) {
    expect_identical(
      check_for_analytical_params(
        case$dist_id, case$params, 3L, case$window, 2
      ), 0L
    )
    for (d in c(0.5, 2, 5, 10)) {
      ode <- primarycensored_numeric_cdf(
        d, case$dist_id, case$params, 2, 3L, case$window
      )
      expect_identical(
        primarycensored_cdf(
          d, case$dist_id, case$params, 2, 0, Inf, 3L, case$window
        ),
        ode
      )
      expect_identical(
        primarycensored_lcdf(
          d, case$dist_id, case$params, 2, 0, Inf, 3L, case$window
        ),
        log(ode)
      )
    }
  }
})

test_that("the analytical function rejects an inadmissible series", {
  expect_error(
    primarycensored_analytical_lcdf(
      2, 2L, c(2.5, 0.4), 2, 0, Inf, 3L, c(0.5, 0.2)
    ),
    "truncated logistic"
  )
  # Also for a delay that would use the small delay form, and a bad scale
  expect_error(
    primarycensored_analytical_lcdf(
      1e-6, 2L, c(2.5, 0.4), 2, 0, Inf, 3L, c(0.5, 0.2)
    ),
    "truncated logistic"
  )
  expect_error(
    primarycensored_analytical_lcdf(
      2, 18L, c(3, 2), 2, 0, Inf, 3L, c(0.5, 0)
    ),
    "truncated logistic"
  )
})

test_that("other delays use the ODE path with a truncated logistic primary", {
  for (d in c(0.5, 2, 5)) {
    expect_identical(
      primarycensored_cdf(d, 3L, c(1.5, 2), 2, 0, Inf, 3L, c(0.5, 0.3)),
      primarycensored_numeric_cdf(d, 3L, c(1.5, 2), 2, 3L, c(0.5, 0.3))
    )
  }
})

test_that("normal delays handle negative delays and truncation", {
  pwindow <- 2
  window <- c(0.5, 0.3)
  cdf <- function(x) pnorm(x, 3, 2)
  ref <- function(x) tlogis_reference(x, pwindow, 0.5, 0.3, cdf, FALSE)
  d <- c(-6, -2, -0.5, 0, 1, 3, 6)
  lcdf <- vapply(
    d, primarycensored_lcdf, numeric(1),
    18L, c(3, 2), pwindow, -Inf, Inf, 3L, window
  )
  expect_lt(max_rel_diff(exp(lcdf), ref(d)), 1e-8)
  for (bounds in list(c(-2, 9), c(-Inf, 9), c(-2, Inf), c(0.5, 7))) {
    L <- bounds[1]
    D <- bounds[2]
    lower <- if (is.finite(L)) ref(L) else 0
    upper <- if (is.finite(D)) ref(D) else 1
    x <- c(1, 3, 5)
    expected <- (ref(x) - lower) / (upper - lower)
    actual <- vapply(
      x, primarycensored_cdf, numeric(1),
      18L, c(3, 2), pwindow, L, D, 3L, window
    )
    expect_equal(actual, expected, tolerance = 1e-8)
    actual_l <- vapply(
      x, primarycensored_lcdf, numeric(1),
      18L, c(3, 2), pwindow, L, D, 3L, window
    )
    expect_equal(exp(actual_l), expected, tolerance = 1e-8)
  }
})

per_delay_tlogis_lcdf <- function(delays, dist_id, params, pwindow, window) {
  lower <- if (dist_id == 18L) -Inf else 0
  vapply(
    delays, primarycensored_lcdf, numeric(1), # nolint: object_usage_linter.
    dist_id, params, pwindow, lower, Inf, 3L, window
  )
}

test_that("check_for_tlogis_vectorized needs an integer pwindow", {
  for (dist_id in c(2L, 4L, 18L)) {
    expect_identical(check_for_tlogis_vectorized(dist_id, 3L, 1), 1L)
    expect_identical(check_for_tlogis_vectorized(dist_id, 3L, 7), 1L)
    expect_identical(check_for_tlogis_vectorized(dist_id, 3L, 1.5), 0L)
    expect_identical(check_for_tlogis_vectorized(dist_id, 3L, 0.5), 0L)
    expect_identical(check_for_tlogis_vectorized(dist_id, 2L, 1), 0L)
  }
  for (dist_id in c(1L, 3L, 26L)) {
    expect_identical(check_for_tlogis_vectorized(dist_id, 3L, 1), 0L)
  }
})

test_that("the vectorised CDF matches the per delay CDF", {
  n <- 31L
  for (case in tlogis_stan_cases) {
    for (pwindow in c(1, 2, 7)) {
      # Locations before, inside (integer and not) and after the window
      for (location in c(-1, 0, 0.4, 1, 2.6, pwindow, pwindow + 2.5)) {
        for (scale in c(0.2, 1.5)) {
          if (!tlogis_case_analytic(case, location, scale, pwindow)) {
            next
          }
          for (start in c(1L, 5L)) {
            window <- c(location, scale)
            vectorised <- primarycensored_tlogis_lcdf_vectorized(
              start, n, case$dist_id, case$params, pwindow, location, scale
            )
            expect_length(vectorised, n)
            expect_equal(
              vectorised[start:n],
              per_delay_tlogis_lcdf(
                start:n, case$dist_id, case$params, pwindow, window
              ),
              tolerance = 1e-12,
              info = tlogis_case_label(
                case,
                pwindow = pwindow, location = location, scale = scale,
                start = start
              )
            )
          }
        }
      }
    }
  }
})

test_that("primarycensored_lcdf_vectorized uses the tlogis shared terms", {
  for (case in tlogis_stan_cases[c(4, 7)]) {
    window <- c(4, 0.5)
    expect_identical(
      primarycensored_lcdf_vectorized(
        1L, 20L, case$dist_id, case$params, 3, 3L, window
      ),
      primarycensored_tlogis_lcdf_vectorized(
        1L, 20L, case$dist_id, case$params, 3, 4, 0.5
      )
    )
  }
  # An inadmissible series and a non-integer window use the per delay path
  expect_identical(
    primarycensored_lcdf_vectorized(1L, 10L, 4L, 0.3, 3, 3L, c(0.5, 0.2)),
    per_delay_tlogis_lcdf(1:10, 4L, 0.3, 3, c(0.5, 0.2))
  )
  expect_identical(
    primarycensored_lcdf_vectorized(1L, 10L, 4L, 0.3, 1.5, 3L, c(2, 0.2)),
    per_delay_tlogis_lcdf(1:10, 4L, 0.3, 1.5, c(2, 0.2))
  )
})

test_that("the vectorised PMF matches the per delay PMF with truncation", {
  settings <- list(
    list(max_delay = 10, L = 0, D = 11),
    list(max_delay = 10, L = 0, D = Inf),
    list(max_delay = 20, L = 2, D = 21),
    list(max_delay = 20, L = -Inf, D = Inf)
  )
  for (case in tlogis_stan_cases[c(4, 7, 8)]) {
    for (setting in settings) {
      for (pwindow in c(1, 3)) {
        for (window in list(c(-0.5, 0.4), c(1, 0.3), c(0.5, 1), c(5, 0.6))) {
          if (!tlogis_case_analytic(case, window[1], window[2], pwindow)) {
            next
          }
          max_delay <- setting$max_delay
          vectorised <- primarycensored_sone_lpmf_vectorized(
            max_delay, setting$L, setting$D, case$dist_id, case$params,
            pwindow, 3L, window
          )
          per_delay <- vapply(
            0:max_delay, function(d) {
              primarycensored_lpmf(
                d, case$dist_id, case$params, pwindow, d + 1,
                setting$L, setting$D, 3L, window
              )
            },
            numeric(1)
          )
          expect_equal(
            vectorised, per_delay,
            tolerance = 1e-9,
            info = tlogis_case_label(
              case,
              pwindow = pwindow, location = window[1], scale = window[2],
              L = setting$L, D = setting$D
            )
          )
        }
      }
    }
  }
})

# The finite differences of CmdStan are noisy to about 1e-6, which sets the
# tolerance. The gamma shape gradient has a looser one.
expect_tlogis_gradient_close <- function(res, case, label) {
  tolerance <- rep(5e-4, 4)
  if (case$dist_id == 2L) {
    tolerance[1] <- 2e-3
  }
  allowed <- tolerance * pmax(abs(res$finite_diff), 1e-2)
  testthat::expect_true(
    all(abs(res$gradient - res$finite_diff) <= allowed),
    info = paste0(
      label, ": gradient ", toString(signif(res$gradient, 5)),
      ", finite difference ", toString(signif(res$finite_diff, 5))
    )
  )
}

tlogis_gradient_cases <- tlogis_stan_cases[c(1, 3, 4, 6, 7, 8)]

test_that("truncated logistic log CDFs have finite gradients matching finite
  differences", {
  model <- tlogis_gradient_model()
  # Locations before, inside and after the window, at several delays
  points <- list(
    list(d = 0.3, pwindow = 2, location = 0.5, scale = 0.4),
    list(d = 1, pwindow = 1, location = -0.5, scale = 0.6),
    list(d = 2.5, pwindow = 2, location = 1.2, scale = 0.3),
    list(d = 6, pwindow = 3, location = 4, scale = 0.7),
    list(d = 20, pwindow = 3, location = 1, scale = 1.5),
    list(d = 2.5, pwindow = 2, location = 3, scale = 0.2),
    list(d = 0.0001, pwindow = 2, location = 0.3, scale = 2),
    list(d = 4, pwindow = 2, location = -2, scale = 0.5)
  )
  for (case in tlogis_gradient_cases) {
    for (point in points) {
      if (!tlogis_case_analytic(
        case, point$location, point$scale, point$pwindow
      )) {
        next
      }
      d <- point$d * if (case$params[length(case$params)] > 100) 1e-3 else 1
      label <- tlogis_case_label(
        case,
        d = d, pwindow = point$pwindow, location = point$location,
        scale = point$scale
      )
      res <- tlogis_gradient_at(
        model, case, d, point$pwindow, point$location, point$scale
      )
      expect_false(res$gradient_not_finite, info = label)
      expect_false(res$rejected, info = label)
      expect_length(res$gradient, 4)
      expect_true(all(is.finite(res$gradient)), info = label)
      expect_tlogis_gradient_close(res, case, label)
    }
  }
})

# Gamma shape gradients at points that use the analytic path, with `compare`
# given the result, the point and a label.
expect_shape_gradients <- function(model, points, compare,
                                   vectorised = FALSE) {
  for (point in points) {
    point <- utils::modifyList(list(pwindow = 1), point)
    case <- list(dist_id = 2L, params = c(point$shape, point$rate))
    label <- tlogis_case_label(
      case,
      d = point$d, pwindow = point$pwindow, location = point$location,
      scale = point$scale
    )
    testthat::expect_true(
      tlogis_case_analytic(case, point$location, point$scale, point$pwindow),
      info = label
    )
    res <- tlogis_gradient_at(
      model, case, point$d, point$pwindow, point$location, point$scale,
      vectorised = vectorised
    )
    testthat::expect_false(res$gradient_not_finite, info = label)
    testthat::expect_false(res$rejected, info = label)
    testthat::expect_true(all(is.finite(res$gradient)), info = label)
    compare(res, point, label)
  }
}

expect_close_to_finite_diff <- function(res, point, label) {
  allowed <- 1e-3 * pmax(abs(res$finite_diff), 1e-2)
  testthat::expect_true(
    all(abs(res$gradient - res$finite_diff) <= allowed),
    info = paste0(
      label, ": gradient ", toString(signif(res$gradient, 5)),
      ", finite difference ", toString(signif(res$finite_diff, 5))
    )
  )
}

# The reference is the gradient of the integral, as CmdStan's finite
# differences are not accurate enough at large shapes. The tolerance is 1e-5
# relative, with a floor for gradients close to zero.
expect_close_to_reference <- function(vectorised = FALSE) {
  function(res, point, label) {
    reference <- tlogis_shape_grad_ref(
      point$shape, point$rate, point$d, point$pwindow, point$location,
      point$scale, if (vectorised) point$d
    )
    testthat::expect_true(
      abs(res$gradient[1] - reference) <= 1e-5 * max(abs(reference), 1e-3),
      info = paste0(
        label, ": gradient ", signif(res$gradient[1], 7),
        ", reference ", signif(reference, 7)
      )
    )
  }
}

test_that("the gamma shape gradient matches finite differences where the
  series cancel and in the lower tail", {
  model <- tlogis_gradient_model()
  expect_shape_gradients(
    model,
    list(
      list(shape = 2.1, rate = 5, d = 1.7, location = 0.5, scale = 50),
      list(shape = 6, rate = 2, d = 4, location = -1, scale = 50),
      list(shape = 6, rate = 2, d = 4, location = 3, scale = 2),
      list(shape = 2, rate = 20, d = 0.3, location = 0.5, scale = 2),
      list(shape = 6, rate = 5, d = 1.5, location = 0.3, scale = 20),
      list(shape = 6, rate = 2, d = 0.05, location = 3, scale = 0.5),
      list(shape = 6, rate = 2, d = 0.02, location = 3, scale = 5),
      list(shape = 20, rate = 2, d = 0.4, location = 3, scale = 0.5),
      list(shape = 20, rate = 2, d = 0.8, location = 3, scale = 5)
    ),
    expect_close_to_finite_diff
  )
})

test_that("the gamma shape gradient matches a reference integral at large
  shapes", {
  model <- tlogis_gradient_model()
  expect_shape_gradients(
    model,
    list(
      list(shape = 150, rate = 50, d = 1.0492872, location = -1, scale = 50),
      list(shape = 150, rate = 50, d = 1.0492872, location = 0.4, scale = 5),
      list(shape = 20, rate = 6.6667, d = 0.7187, location = -1, scale = 50),
      list(shape = 80, rate = 26.667, d = 0.95, location = 0.4, scale = 5),
      list(shape = 80, rate = 26.667, d = 3.9, location = 3, scale = 0.5),
      list(shape = 80, rate = 26.667, d = 4.1, location = -1, scale = 50),
      list(shape = 150, rate = 50, d = 4.25, location = 3, scale = 0.5),
      list(shape = 150, rate = 50, d = 4.5, location = 2, scale = 0.05),
      list(shape = 250, rate = 83.333, d = 4.6, location = 3, scale = 0.5),
      list(shape = 250, rate = 83.333, d = 4.8, location = -1, scale = 50),
      list(shape = 400, rate = 133.3, d = 4.6, location = 3, scale = 0.5)
    ),
    expect_close_to_reference()
  )
})

test_that("the vectorised gamma log PMF has a finite, accurate gradient in the
  shape at large shapes", {
  model <- tlogis_gradient_model()
  expect_shape_gradients(
    model,
    list(
      list(shape = 60, rate = 12, d = 9, location = 1, scale = 30),
      list(shape = 80, rate = 16, d = 9, location = 3, scale = 0.5),
      list(shape = 150, rate = 30, d = 9, location = -1, scale = 50),
      list(shape = 250, rate = 27, d = 9, location = 3, scale = 0.5)
    ),
    expect_close_to_reference(vectorised = TRUE),
    vectorised = TRUE
  )
})

test_that("the vectorised truncated logistic log PMF has finite gradients
  matching finite differences", {
  model <- tlogis_gradient_model()
  points <- list(
    list(d = 12, pwindow = 3, location = 4, scale = 0.7),
    list(d = 12, pwindow = 3, location = -0.5, scale = 0.7),
    list(d = 10, pwindow = 2, location = 1.3, scale = 0.4),
    list(d = 10, pwindow = 2, location = 1, scale = 1)
  )
  for (case in tlogis_gradient_cases[c(3, 5, 6)]) {
    for (point in points) {
      if (!tlogis_case_analytic(
        case, point$location, point$scale, point$pwindow
      )) {
        next
      }
      label <- tlogis_case_label(
        case,
        d = point$d, pwindow = point$pwindow, location = point$location,
        scale = point$scale
      )
      res <- tlogis_gradient_at(
        model, case, point$d, point$pwindow, point$location, point$scale,
        vectorised = TRUE
      )
      expect_false(res$gradient_not_finite, info = label)
      expect_false(res$rejected, info = label)
      expect_true(all(is.finite(res$gradient)), info = label)
      expect_tlogis_gradient_close(res, case, label)
    }
  }
})
