skip_on_cran()

# Stan solutions for exponential, gamma and normal delays with a truncated
# logistic primary (primary_id 3 with primary_params = c(location, scale)).
# These tests check the Stan series against the R implementation and a
# reference integral, the dispatch and fallback to the ODE path, the shared
# endpoint vectorised form, and gradients.

tlogis_stan_cases <- list(
  list(dist_id = 4L, params = 2, pdist = pexp, args = list(rate = 2)),
  list(dist_id = 4L, params = 0.3, pdist = pexp, args = list(rate = 0.3)),
  list(
    dist_id = 2L, params = c(0.6, 1.3), pdist = pgamma,
    args = list(shape = 0.6, rate = 1.3)
  ),
  list(
    dist_id = 2L, params = c(2.5, 0.4), pdist = pgamma,
    args = list(shape = 2.5, rate = 0.4)
  ),
  list(
    dist_id = 2L, params = c(20, 4), pdist = pgamma,
    args = list(shape = 20, rate = 4)
  ),
  list(
    dist_id = 2L, params = c(2.5, 1000), pdist = pgamma,
    args = list(shape = 2.5, rate = 1000)
  ),
  list(
    dist_id = 18L, params = c(3, 2), pdist = pnorm,
    args = list(mean = 3, sd = 2)
  ),
  list(
    dist_id = 18L, params = c(-1, 3), pdist = pnorm,
    args = list(mean = -1, sd = 3)
  )
)

tlogis_case_label <- function(case, ...) {
  paste0(
    "dist ", case$dist_id, " params ", toString(case$params), ", ",
    paste(names(list(...)), unlist(list(...)), sep = " = ", collapse = ", ")
  )
}

tlogis_case_cdf <- function(case) {
  function(x) do.call(case$pdist, c(list(x), case$args))
}

# Internal lower bound used by primarycensored_lcdf for each support
tlogis_case_lower <- function(case) {
  if (case$dist_id == 18L) -Inf else 0
}

tlogis_case_obj <- function(case, location, scale) {
  tlogis_object(list(pdist = case$pdist, args = case$args), location, scale)
}

# The window, location and scale sets of the acceptance criteria, and
# delays below and above the window. Stan returns -Inf for probabilities
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
  # The other solutions are unchanged
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
  # check_for_analytical_window adds it to the structural checks
  expect_identical(
    check_for_analytical_window(2L, c(2.5, 0.4), 3L, c(3, 0.5), 2), 1L
  )
  expect_identical(
    check_for_analytical_window(2L, c(2.5, 0.4), 3L, c(0.5, 0.2), 2), 0L
  )
  expect_identical(
    check_for_analytical_window(3L, c(2.5, 0.4), 3L, c(3, 0.5), 2), 0L
  )
  # Unchanged for the other primaries
  expect_identical(
    check_for_analytical_window(2L, c(2.5, 0.4), 1L, numeric(0), 2), 1L
  )
  expect_identical(
    check_for_analytical_window(2L, c(2.5, 0.4), 2L, -0.4, 2), 0L
  )
  expect_identical(
    check_for_analytical_window(2L, c(2.5, 0.4), 2L, 0.5, 2), 1L
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
          check_for_analytical_window(
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
      check_for_analytical_window(
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

# Narrow windows. The integrand of the ODE has a spike of width about
# `scale` where the delay is the window location, and a solver that takes
# large steps across the flat region steps over it. The ODE is solved to a
# relative tolerance of 1e-9 and an absolute tolerance of 1e-10 for this
# primary. The error is an absolute one of about 1e-9 to 1e-7, so these are
# compared in absolute terms at 1e-7 and small CDFs and PMFs are not compared
# in relative terms.
tlogis_narrow_cases <- list(
  list(dist_id = 4L, params = 1.5, pdist = pexp, args = list(rate = 1.5)),
  list(
    dist_id = 2L, params = c(2, 1), pdist = pgamma,
    args = list(shape = 2, rate = 1)
  ),
  list(
    dist_id = 1L, params = c(1, 0.5), pdist = plnorm,
    args = list(meanlog = 1, sdlog = 0.5)
  ),
  list(
    dist_id = 1L, params = c(0, 0.5), pdist = plnorm,
    args = list(meanlog = 0, sdlog = 0.5)
  ),
  list(
    dist_id = 3L, params = c(2, 2), pdist = pweibull,
    args = list(shape = 2, scale = 2)
  ),
  list(
    dist_id = 18L, params = c(3, 2), pdist = pnorm,
    args = list(mean = 3, sd = 2)
  )
)

test_that("the ODE path resolves a narrow truncated logistic primary", {
  pwindow <- 2
  d <- c(0.5, 1, 3, 10)
  for (case in tlogis_narrow_cases) {
    cdf <- tlogis_case_cdf(case)
    for (location in c(-0.5, 0.7, 2.5)) {
      for (scale in c(0.02, 0.005, 0.001)) {
        window <- c(location, scale)
        expected <- tlogis_reference(
          d, pwindow, location, scale, cdf, case$dist_id != 18L
        )
        ode <- vapply(
          d, primarycensored_numeric_cdf, numeric(1),
          case$dist_id, case$params, pwindow, 3L, window
        )
        expect_lt(
          max(abs(ode - expected)), 1e-7,
          label = tlogis_case_label(
            case,
            location = location, scale = scale
          )
        )
      }
    }
  }
})

test_that("the log CDF and log PMF are correct for a narrow primary", {
  pwindow <- 2
  for (case in tlogis_narrow_cases) {
    cdf <- tlogis_case_cdf(case)
    lower <- tlogis_case_lower(case)
    for (location in c(0.7, 1)) {
      for (scale in c(0.02, 0.005)) {
        window <- c(location, scale)
        label <- tlogis_case_label(case, location = location, scale = scale)
        d <- c(0.5, 1, 3, 10)
        expected <- tlogis_reference(
          d, pwindow, location, scale, cdf, case$dist_id != 18L
        )
        lcdf <- vapply(
          d, primarycensored_lcdf, numeric(1),
          case$dist_id, case$params, pwindow, lower, Inf, 3L, window
        )
        expect_lt(max(abs(exp(lcdf) - expected)), 1e-7, label = label)
        # The PMF over integer delays
        ref <- tlogis_reference(
          0:6, pwindow, location, scale, cdf, case$dist_id != 18L
        )
        pmf <- exp(primarycensored_sone_lpmf_vectorized(
          5, lower, Inf, case$dist_id, case$params, pwindow, 3L, window
        ))
        expect_lt(
          max(abs(pmf - diff(ref)[1:6])), 1e-7,
          label = paste("pmf:", label)
        )
      }
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
})

test_that("other delays use the ODE path with a truncated logistic primary", {
  for (d in c(0.5, 2, 5)) {
    expect_identical(
      primarycensored_cdf(d, 3L, c(1.5, 2), 2, 0, Inf, 3L, c(0.5, 0.3)),
      primarycensored_numeric_cdf(d, 3L, c(1.5, 2), 2, 3L, c(0.5, 0.3))
    )
  }
})

test_that("non-parametric delays are analytic with a truncated logistic
  primary", {
  boundaries <- c(0, 1, 3, 6, 10)
  pmf <- c(0.2, 0.3, 0.35, 0.15)
  params <- c(boundaries, pmf)
  window <- c(0.5, 0.3)
  obj <- new_pcens(
    pdist = pdiscretestep, dprimary = dtlogis,
    primary_args = list(location = 0.5, scale = 0.3),
    boundaries = boundaries, pmf = pmf
  )
  d <- c(0.5, 1.5, 3, 5, 8, 12)
  expect_identical(check_for_analytical(26L, 3L), 1L)
  analytic <- vapply(
    d, primarycensored_cdf, numeric(1), 26L, params, 2, -Inf, Inf, 3L, window
  )
  expect_equal(analytic, pcens_cdf(obj, d, 2), tolerance = 1e-9)
  # The step CDF has kinks, so the ODE path is accurate to about 1e-4
  ode <- vapply(
    d, primarycensored_numeric_cdf, numeric(1), 26L, params, 2, 3L, window
  )
  expect_lt(max(abs(analytic - ode)), 1e-4)
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

# Gradients are only observable from a compiled model, so this builds a
# minimal one whose target is the log CDF or the vectorised log PMF, and runs
# `stan_gradient_at()` from helper-stan-gradient.R.
tlogis_gradient_model <- function() {
  testthat::skip_if_not_installed("cmdstanr")
  testthat::skip_if(
    is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))
  )
  functions <- pcd_load_stan_functions(
    wrap_in_block = TRUE, write_to_file = FALSE
  )
  code <- paste0(
    functions, "\n",
    "data {\n",
    "  int dist_id;\n",
    "  int n_params;\n",
    "  int vectorised;\n",
    "  real d;\n",
    "  real pwindow;\n",
    "  real L;\n",
    "}\n",
    "parameters {\n",
    "  real p1;\n",
    "  real<lower=0> p2;\n",
    "  real location;\n",
    "  real<lower=0> scale;\n",
    "}\n",
    "model {\n",
    "  array[2] real all_params = {p1, p2};\n",
    "  array[n_params] real params = all_params[1:n_params];\n",
    "  if (vectorised) {\n",
    "    target += sum(primarycensored_sone_lpmf_vectorized(\n",
    "      to_int(d), L, positive_infinity(), dist_id, params, pwindow, 3,\n",
    "      {location, scale}\n",
    "    ));\n",
    "  } else {\n",
    "    target += primarycensored_lcdf(\n",
    "      d | dist_id, params, pwindow, L, positive_infinity(), 3,\n",
    "      {location, scale}\n",
    "    );\n",
    "  }\n",
    "}\n"
  )
  path <- file.path(tempdir(), "pcd_tlogis_cdf_gradient.stan")
  writeLines(code, path)
  suppressMessages(suppressWarnings(cmdstanr::cmdstan_model(path)))
}

tlogis_gradient_at <- function(model, case, d, pwindow, location, scale,
                               vectorised = FALSE) {
  init <- list(
    p1 = case$params[1],
    p2 = if (length(case$params) > 1) case$params[2] else 1,
    location = location, scale = scale
  )
  stan_gradient_at( # nolint: object_usage_linter.
    model,
    data = list(
      dist_id = case$dist_id, n_params = length(case$params),
      vectorised = as.integer(vectorised), d = d, pwindow = pwindow,
      L = tlogis_case_lower(case)
    ),
    init = init
  )
}

# The finite differences of CmdStan use a step of 1e-6, so the 1e-12 relative
# error of the truncated series shows as noise of about 1e-6 in them, which
# sets the tolerance. The gradient in the shape of a gamma delay comes from
# Stan's gamma_lcdf(), which is accurate to 1e-7 in the body but not far into
# the lower tail, so it has a looser tolerance. Exponential delays have no
# shape.
expect_tlogis_gradient_close <- function(res, case, label) {
  tolerance <- rep(5e-4, 4)
  if (case$dist_id == 2L) {
    tolerance[1] <- 2e-3
  }
  # A parameter the delay does not have has a zero gradient
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

test_that("the gamma shape gradient is accurate where the series cancel", {
  model <- tlogis_gradient_model()
  # A scale that is large relative to the window makes the series cancel
  # heavily, which amplifies any error in the gradient of the shape of the
  # terms. These points had errors of 2% to 7%, and 2.6 times at the last.
  points <- list(
    list(shape = 2.1, rate = 5, d = 1.7, pwindow = 1, location = 0.5, scale = 50),
    list(shape = 6, rate = 2, d = 4, pwindow = 1, location = -1, scale = 50),
    list(shape = 6, rate = 2, d = 4, pwindow = 1, location = 3, scale = 2),
    list(shape = 2, rate = 20, d = 0.3, pwindow = 1, location = 0.5, scale = 2),
    list(shape = 6, rate = 5, d = 1.5, pwindow = 1, location = 0.3, scale = 20)
  )
  for (point in points) {
    case <- list(dist_id = 2L, params = c(point$shape, point$rate))
    label <- tlogis_case_label(
      case,
      d = point$d, pwindow = point$pwindow, location = point$location,
      scale = point$scale
    )
    expect_true(
      tlogis_case_analytic(
        case, point$location, point$scale, point$pwindow
      ),
      info = label
    )
    res <- tlogis_gradient_at(
      model, case, point$d, point$pwindow, point$location, point$scale
    )
    expect_false(res$gradient_not_finite, info = label)
    expect_false(res$rejected, info = label)
    expect_true(all(is.finite(res$gradient)), info = label)
    allowed <- 1e-3 * pmax(abs(res$finite_diff), 1e-2)
    expect_true(
      all(abs(res$gradient - res$finite_diff) <= allowed),
      info = paste0(
        label, ": gradient ", toString(signif(res$gradient, 5)),
        ", finite difference ", toString(signif(res$finite_diff, 5))
      )
    )
  }
})

test_that("the gamma shape gradient is accurate far into the lower tail", {
  model <- tlogis_gradient_model()
  # Stan's gradient of gamma_lcdf() in the shape truncates its series at an
  # absolute tolerance, which is a relative error of up to 80% when the
  # probability is below 1e-10. The log CDFs here are -25 to -115. A
  # location after the window keeps the solution analytic.
  points <- list(
    list(shape = 6, rate = 2, d = 0.05, scale = 0.5),
    list(shape = 6, rate = 2, d = 0.02, scale = 5),
    list(shape = 6, rate = 2, d = 0.05, scale = 5),
    list(shape = 20, rate = 2, d = 0.4, scale = 0.5),
    list(shape = 20, rate = 2, d = 0.4, scale = 5),
    list(shape = 20, rate = 2, d = 0.8, scale = 5)
  )
  for (point in points) {
    case <- list(dist_id = 2L, params = c(point$shape, point$rate))
    label <- tlogis_case_label(
      case,
      d = point$d, pwindow = 1, location = 3, scale = point$scale
    )
    expect_true(
      tlogis_case_analytic(case, 3, point$scale, 1),
      info = label
    )
    res <- tlogis_gradient_at(model, case, point$d, 1, 3, point$scale)
    expect_false(res$gradient_not_finite, info = label)
    expect_true(all(is.finite(res$gradient)), info = label)
    allowed <- 1e-3 * pmax(abs(res$finite_diff), 1e-2)
    expect_true(
      all(abs(res$gradient - res$finite_diff) <= allowed),
      info = paste0(
        label, ": gradient ", toString(signif(res$gradient, 5)),
        ", finite difference ", toString(signif(res$finite_diff, 5))
      )
    )
  }
})

# The gradient in the shape of a gamma delay is compared to the gradient of
# the reference integral, not to CmdStan's finite differences, which have
# errors of 2e-4 where the log CDF is large. The tolerance is 1e-5 relative
# to the gradient, with a floor for gradients that are close to zero.
expect_shape_grad_close <- function(model, point, label,
                                    tolerance = 1e-5,
                                    vectorised = FALSE) {
  case <- list(dist_id = 2L, params = c(point$shape, point$rate))
  max_delay <- if (vectorised) point$d else NULL
  res <- tlogis_gradient_at(
    model, case, point$d, point$pwindow, point$location, point$scale,
    vectorised = vectorised
  )
  testthat::expect_false(res$gradient_not_finite, info = label)
  testthat::expect_false(res$rejected, info = label)
  testthat::expect_true(all(is.finite(res$gradient)), info = label)
  reference <- tlogis_shape_grad_ref(
    point$shape, point$rate, point$d, point$pwindow, point$location,
    point$scale, max_delay
  )
  testthat::expect_true(
    abs(res$gradient[1] - reference) <=
      tolerance * max(abs(reference), 1e-3),
    info = paste0(
      label, ": gradient ", signif(res$gradient[1], 7),
      ", reference ", signif(reference, 7)
    )
  )
}

test_that("the gamma shape gradient is accurate where the lower tail series
  of the gamma CDF meets the cancellation of the series", {
  model <- tlogis_gradient_model()
  # The log CDFs are -19 to -68. The first term of the series of the gamma
  # CDF is 1e-5 to 1 here, where Stan's gamma_lcdf() has a gradient in the
  # shape with a relative error of up to 6e-6 that the cancellation amplified
  # to 1e-3 at shape 20 and 4.6% at shape 150. A location before the window
  # and a scale of 50 cancels the most, and a scale of 5 less.
  points <- list(
    list(shape = 150, rate = 50, d = 1.0492872, location = -1, scale = 50),
    list(shape = 150, rate = 50, d = 1.0492872, location = 0.4, scale = 5),
    list(shape = 20, rate = 6.6667, d = 0.7187, location = -1, scale = 50),
    list(shape = 20, rate = 6.6667, d = 0.7187, location = 0.4, scale = 5),
    list(shape = 80, rate = 26.667, d = 1, location = -1, scale = 50),
    list(shape = 80, rate = 26.667, d = 0.95, location = 0.4, scale = 5)
  )
  for (point in points) {
    point$pwindow <- 1
    case <- list(dist_id = 2L, params = c(point$shape, point$rate))
    label <- tlogis_case_label(
      case,
      d = point$d, pwindow = 1, location = point$location,
      scale = point$scale
    )
    expect_true(
      tlogis_case_analytic(case, point$location, point$scale, 1),
      info = label
    )
    expect_shape_grad_close(model, point, label)
  }
})

test_that("the gamma shape gradient is finite and accurate in the upper tail
  at large shapes", {
  model <- tlogis_gradient_model()
  # The upper tail of the delay is 1e-7 to 1e-16 here, where Stan's
  # gamma_lccdf() has a gradient in the shape with a relative error of
  # 5e-4 at shape 80 and 2e-3 at shape 150, and NaN from shape 200.
  points <- list(
    list(shape = 80, rate = 26.667, d = 3.9, location = 3, scale = 0.5),
    list(shape = 80, rate = 26.667, d = 4.1, location = -1, scale = 50),
    list(shape = 150, rate = 50, d = 4.25, location = 3, scale = 0.5),
    list(shape = 150, rate = 50, d = 4.5, location = 2, scale = 0.05),
    list(shape = 250, rate = 83.333, d = 4.6, location = 3, scale = 0.5),
    list(shape = 250, rate = 83.333, d = 4.8, location = -1, scale = 50),
    list(shape = 400, rate = 133.3, d = 4.6, location = 3, scale = 0.5)
  )
  for (point in points) {
    point$pwindow <- 1
    case <- list(dist_id = 2L, params = c(point$shape, point$rate))
    label <- tlogis_case_label(
      case,
      d = point$d, pwindow = 1, location = point$location,
      scale = point$scale
    )
    expect_true(
      tlogis_case_analytic(case, point$location, point$scale, 1),
      info = label
    )
    expect_shape_grad_close(model, point, label)
  }
})

test_that("the vectorised gamma log PMF has a finite, accurate gradient in the
  shape at large shapes", {
  model <- tlogis_gradient_model()
  points <- list(
    list(shape = 60, rate = 12, d = 9, location = 1, scale = 30),
    list(shape = 80, rate = 16, d = 9, location = 3, scale = 0.5),
    list(shape = 150, rate = 30, d = 9, location = -1, scale = 50),
    list(shape = 250, rate = 27, d = 9, location = 3, scale = 0.5)
  )
  for (point in points) {
    point$pwindow <- 1
    case <- list(dist_id = 2L, params = c(point$shape, point$rate))
    label <- tlogis_case_label(
      case,
      d = point$d, pwindow = 1, location = point$location,
      scale = point$scale
    )
    expect_true(
      tlogis_case_analytic(case, point$location, point$scale, 1),
      info = label
    )
    expect_shape_grad_close(
      model, point, label,
      vectorised = TRUE
    )
  }
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

test_that("the ODE path with a narrow truncated logistic primary has finite
  gradients matching finite differences", {
  model <- tlogis_gradient_model()
  # The gradients pass through the times at which the integral is split,
  # which must cancel. The tolerance is that of the gradient of the ODE.
  cases <- list(tlogis_stan_cases[[4]], tlogis_stan_cases[[2]])
  points <- list(
    list(d = 3, pwindow = 2, location = 0.7, scale = 0.05),
    list(d = 1, pwindow = 2, location = 0.7, scale = 0.02),
    list(d = 3, pwindow = 2, location = -0.5, scale = 0.05),
    list(d = 4, pwindow = 2, location = 2.5, scale = 0.05)
  )
  for (case in cases) {
    for (point in points) {
      # Only the points that use the ODE
      if (tlogis_case_analytic(
        case, point$location, point$scale, point$pwindow
      )) {
        next
      }
      for (vectorised in c(FALSE, TRUE)) {
        label <- tlogis_case_label(
          case,
          d = point$d, pwindow = point$pwindow, location = point$location,
          scale = point$scale, vectorised = vectorised
        )
        res <- tlogis_gradient_at(
          model, case, point$d, point$pwindow, point$location, point$scale,
          vectorised = vectorised
        )
        expect_false(res$gradient_not_finite, info = label)
        expect_false(res$rejected, info = label)
        expect_true(all(is.finite(res$gradient)), info = label)
        allowed <- 1e-2 * pmax(abs(res$finite_diff), 1e-2)
        expect_true(
          all(abs(res$gradient - res$finite_diff) <= allowed),
          info = paste0(
            label, ": gradient ", toString(signif(res$gradient, 5)),
            ", finite difference ", toString(signif(res$finite_diff, 5))
          )
        )
      }
    }
  }
})
