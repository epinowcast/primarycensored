skip_on_cran()

# Stan solution for the lognormal delay (dist_id 1) with an exponentially
# tilted primary (primary_id 2 with primary_params = r). The transform is
# evaluated by Gauss-Legendre panels for a tilt xi = -r < 0 and by a series of
# partial moments for xi > 0. These tests check the transforms against the R
# implementation and a reference integral, the CDF against the reference
# integral and the ODE path, the shared endpoint vectorised form, and
# gradients.

lnorm_stan_cases <- exptilt_lnorm_cases()

lnorm_params <- function(case) c(case$meanlog, case$sdlog)

lnorm_object_stan <- function(case, rho) {
  exptilt_object(exptilt_lnorm_family(case), rho)
}

lnorm_stan_lcdf <- function(d, case, pwindow, rho) {
  vapply(
    d, primarycensored_exptilt_lcdf, numeric(1), # nolint: object_usage_linter.
    1L, lnorm_params(case), pwindow, rho
  )
}

test_that("the lognormal has a tilt transform for every tilt that does not
  overflow", {
  expect_identical(check_for_tilt_transform(1L, -1, c(1.6, 0.5)), 1L)
  expect_identical(check_for_tilt_transform(1L, -50, c(1.6, 0.5)), 1L)
  expect_identical(check_for_tilt_transform(1L, 0, c(1.6, 0.5)), 1L)
  expect_identical(check_for_tilt_transform(1L, 50, c(1.6, 0.5)), 1L)
  expect_identical(check_for_tilt_transform(1L, -1e300, c(650, 1)), 0L)
  expect_identical(check_for_tilt_transform(1L, -1, c(1.6, 0)), 0L)
})

test_that("check_for_analytical includes the lognormal with a tilted
  primary", {
  expect_identical(check_for_analytical(1L, 2L), 1L)
  expect_identical(check_for_exptilt(1L, 2L), 1L)
  expect_identical(check_for_exptilt(1L, 1L), 0L)
  expect_identical(check_for_exptilt_vectorized(1L, 2L, 1), 1L)
  expect_identical(check_for_exptilt_vectorized(1L, 2L, 1.5), 0L)
  # The uniform primary solution is unchanged
  expect_identical(check_for_analytical(1L, 1L), 1L)
  expect_identical(check_for_uniform_terms(1L, 1L), 1L)
  expect_identical(check_for_uniform_terms(1L, 2L), 0L)
  for (rho in c(-2, -0.3, 0.3, 2)) {
    expect_identical(
      check_for_analytical_params(1L, c(1.6, 0.5), 2L, rho), 1L
    )
  }
  expect_identical(
    check_for_analytical_params(1L, c(650, 1), 2L, 1e300), 0L
  )
})

test_that("the shared setup of the tilt transform gives the same pair", {
  params <- c(1.6, 0.5)
  t <- c(-1, 0, 0.05, 0.8, 2, 7, 40)
  for (xi in c(-3, -0.2, 0, 0.3)) {
    context <- log_tilt_transform_context(1L, xi, params)
    expect_length(context, 4)
    for (point in t) {
      expect_identical(
        log_tilt_transform_pair_shared(point, 1L, xi, params, context),
        log_tilt_transform_pair(point, 1L, xi, params)
      )
    }
    # The setup is the mode, the limits and the log total for xi < 0
    if (xi < 0) {
      expect_lt(context[2], context[1])
      expect_lt(context[1], context[3])
      expected <- lnorm_tilt_integral(-45, 45, 1.6, 0.5, xi)
      expect_equal(context[4], expected, tolerance = 1e-10)
    } else {
      expect_identical(context, rep(0, 4))
    }
  }
  # Delays without a setup have an empty one
  expect_identical(log_tilt_transform_context(2L, -0.2, c(2, 1)), rep(0, 4))
  expect_identical(log_tilt_transform_context(18L, -0.2, c(2, 1)), rep(0, 4))
  expect_identical(
    log_tilt_transform_pair_shared(1.5, 2L, -0.2, c(2, 1), rep(0, 4)),
    log_tilt_transform_pair(1.5, 2L, -0.2, c(2, 1))
  )
})

test_that("Stan lognormal tilt transforms match a reference integral", {
  z <- c(-8, -3, -1, 0, 1, 3)
  for (case in lnorm_stan_cases[c(1, 2, 4)]) {
    t <- exp(case$meanlog + case$sdlog * z)
    params <- lnorm_params(case)
    for (xi in c(-4, -0.3, -1e-4, 0.05, 0.6)) {
      lower <- vapply(t, log_tilt_transform, numeric(1), 1L, xi, params)
      expected <- lnorm_tilt_reference(t, case$meanlog, case$sdlog, xi)
      keep <- expected > -700
      expect_lt(
        max(abs(lower[keep] - expected[keep])), 1e-9,
        label = sprintf("xi %g", xi)
      )
      if (xi < 0) {
        upper <- vapply(t, log_tilt_transform_upper, numeric(1), 1L, xi, params)
        expected <- lnorm_tilt_reference(
          t, case$meanlog, case$sdlog, xi,
          upper = TRUE
        )
        keep <- expected > -700
        expect_lt(
          max(abs(upper[keep] - expected[keep])), 1e-9,
          label = sprintf("xi %g", xi)
        )
      }
    }
  }
})

test_that("Stan lognormal tilt transforms are 0 or total below the
  support", {
  params <- c(1.6, 0.5)
  for (t in c(-3, -1e-9, 0)) {
    expect_identical(log_tilt_transform(t, 1L, -0.3, params), -Inf)
    expect_identical(log_tilt_transform(t, 1L, 0.3, params), -Inf)
    expect_identical(log_tilt_transform(t, 1L, 0, params), -Inf)
    expect_identical(log_tilt_transform_upper(t, 1L, 0.3, params), Inf)
    expect_identical(log_tilt_transform_upper(t, 1L, 0, params), 0)
  }
  expected <- lnorm_tilt_integral(-45, 45, 1.6, 0.5, -0.3)
  expect_equal(
    log_tilt_transform_upper(0, 1L, -0.3, params), expected,
    tolerance = 1e-10
  )
  expect_identical(
    log_tilt_transform_upper(-2, 1L, -0.3, params),
    log_tilt_transform_upper(0, 1L, -0.3, params)
  )
})

test_that("Stan lognormal moments match the R moments", {
  ts <- c(1e-3, 0.4, 1, 3.5, 12, 40)
  for (case in lnorm_stan_cases) {
    obj <- lnorm_object_stan(case, 0.1)
    params <- lnorm_params(case)
    expected <- .pcens_tilt_moments(obj, ts)
    actual <- t(vapply(
      ts, primarycensored_tilt_moments, numeric(2), 1L, params
    ))
    keep <- is.finite(expected[, 1]) & expected[, 1] > -700
    expect_equal(
      actual[keep, ], unname(expected[keep, ]),
      tolerance = 1e-9
    )
  }
  params <- c(1.6, 0.5)
  expect_identical(primarycensored_tilt_moments(0, 1L, params), c(-Inf, -Inf))
  expect_identical(primarycensored_tilt_moments(-2, 1L, params), c(-Inf, -Inf))
})

test_that("the lognormal tilted CDF matches a reference integral", {
  for (case in lnorm_stan_cases) {
    cdf <- exptilt_lnorm_cdf(case)
    for (pwindow in c(0.5, 1, 2, 7)) {
      d <- sort(c(
        1e-6, 1e-3, 0.3 * pwindow, pwindow - 1e-3, pwindow, pwindow + 1e-3,
        2, 3, 6, 12, 25
      ))
      for (rho in c(
        -1, -0.5, -0.05, -1e-4, -1e-5, -1e-8, 1e-8, 1e-5, 1e-4,
        0.05, 0.5, 1
      )) {
        expected <- exptilt_reference(d, pwindow, rho, cdf)
        actual <- exp(lnorm_stan_lcdf(d, case, pwindow, rho))
        expect_lt(
          max_rel_diff(actual, expected), 1e-7,
          label = sprintf(
            "meanlog %g, sdlog %g, pwindow %g, r %g",
            case$meanlog, case$sdlog, pwindow, rho
          )
        )
      }
    }
  }
})

test_that("the lognormal tilted CDF matches the R implementation", {
  d <- c(1e-4, 0.3, 1, 2.5, 6, 15, 30)
  for (case in lnorm_stan_cases) {
    for (pwindow in c(0.5, 2, 7)) {
      for (rho in c(-0.2, -1e-6, 1e-6, 0.3)) {
        obj <- lnorm_object_stan(case, rho)
        expect_equal(
          exp(lnorm_stan_lcdf(d, case, pwindow, rho)),
          .pcens_cdf_exptilt(obj, d, pwindow),
          tolerance = 1e-9,
          info = sprintf(
            "meanlog %g, sdlog %g, pwindow %g, r %g",
            case$meanlog, case$sdlog, pwindow, rho
          )
        )
      }
    }
  }
})

test_that("the lognormal tilted CDF is continuous across the small tilt
  forms", {
  d <- c(1e-4, 0.3, 1, 2.5, 6, 15, 30)
  for (case in lnorm_stan_cases) {
    for (pwindow in c(0.5, 2, 7)) {
      for (sign in c(-1, 1)) {
        expect_lt(
          max_rel_diff(
            exp(lnorm_stan_lcdf(d, case, pwindow, sign * 0.9999e-4 / pwindow)),
            exp(lnorm_stan_lcdf(d, case, pwindow, sign * 1.0001e-4 / pwindow))
          ),
          1e-7,
          label = sprintf(
            "meanlog %g, sdlog %g, pwindow %g, sign %g",
            case$meanlog, case$sdlog, pwindow, sign
          )
        )
      }
    }
  }
})

test_that("primarycensored_lcdf and primarycensored_cdf use the lognormal
  solution and agree with the ODE path", {
  d <- c(0.2, 1, 2.5, 6, 15)
  for (case in lnorm_stan_cases) {
    cdf <- exptilt_lnorm_cdf(case)
    params <- lnorm_params(case)
    for (pwindow in c(1, 3)) {
      for (rho in c(-1, -0.5, -1e-8, 1e-8, 0.5, 1)) {
        info <- sprintf(
          "meanlog %g, sdlog %g, pwindow %g, r %g",
          case$meanlog, case$sdlog, pwindow, rho
        )
        expected <- exptilt_reference(d, pwindow, rho, cdf)
        lcdf <- vapply(
          d, primarycensored_lcdf, numeric(1),
          1L, params, pwindow, 0, Inf, 2L, rho
        )
        expect_lt(max_rel_diff(exp(lcdf), expected), 1e-7, label = info)
        plain <- vapply(
          d, primarycensored_cdf, numeric(1),
          1L, params, pwindow, 0, Inf, 2L, rho
        )
        expect_lt(max_rel_diff(plain, expected), 1e-7, label = info)
        # The ODE path has absolute and relative tolerances of 1e-6
        ode <- vapply(
          d, primarycensored_numeric_cdf, numeric(1),
          1L, params, pwindow, 2L, rho
        )
        expect_lt(max(abs(plain - ode)), 1e-4, label = info)
      }
    }
  }
})

lnorm_per_delay_lcdf <- function(delays, params, pwindow, rho) {
  vapply(
    delays, primarycensored_lcdf, numeric(1), # nolint: object_usage_linter.
    1L, params, pwindow, 0, Inf, 2L, rho
  )
}

test_that("the lognormal tilted CDF is accurate far from the origin", {
  pwindow <- 1
  for (m in c(1e5, 1e6, 1e7, 1e8)) {
    d <- m * c(0.5, 1, 2)
    case <- list(meanlog = log(m), sdlog = 0.5)
    for (rho in c(3e-5, -3e-5, 1e-6)) {
      expected <- exptilt_reference(d, pwindow, rho, exptilt_lnorm_cdf(case))
      expect_lt(
        max_rel_diff(
          exp(lnorm_stan_lcdf(d, case, pwindow, rho)), expected
        ),
        1e-7,
        label = sprintf("m = %g, r = %g", m, rho)
      )
    }
  }
})

test_that("the lognormal series limit depends on the tilt and the delay", {
  params <- c(1.6, 0.5)
  # The series is used up to xi t of 60, where it costs as much as an ODE at
  # a tolerance of 1e-10, see the NEWS
  expect_identical(
    vapply(
      c(1, 59, 61, 1e3, 1e6),
      function(t) check_for_tilt_transform_at(1L, 1, params, t, 1),
      integer(1)
    ),
    c(1L, 1L, 0L, 0L, 0L)
  )
  # A window that is wide in tilt terms (xi w above 2) keeps the series,
  # where the ODE is less accurate, up to xi t + 9 sqrt(xi t) + 30 terms, at
  # most 20000
  expect_identical(
    vapply(
      c(61, 1e3, 1.8e4, 1.9e4, 1e6),
      function(t) check_for_tilt_transform_at(1L, 1, params, t, 3),
      integer(1)
    ),
    c(1L, 1L, 1L, 0L, 0L)
  )
  expect_identical(check_for_tilt_transform_at(1L, 1, params, 61, 2), 0L)
  expect_identical(check_for_tilt_transform_at(1L, 1, params, 61, 2.01), 1L)
  # A negative tilt, a small delay and other delays have no limit
  expect_identical(check_for_tilt_transform_at(1L, 1, params, 0, 1), 1L)
  expect_identical(check_for_tilt_transform_at(1L, 1, params, -5, 1), 1L)
  expect_identical(check_for_tilt_transform_at(1L, -1, params, 1e9, 1), 1L)
  expect_identical(
    check_for_tilt_transform_at(18L, 1, c(3, 2), 1e9, 1), 1L
  )
  expect_identical(
    check_for_tilt_transform_at(2L, -3, c(2, 0.4), 1e9, 1), 1L
  )
  # The transform must still exist
  expect_identical(
    check_for_tilt_transform_at(1L, -1e300, c(650, 1), 1, 1), 0L
  )
  expect_identical(check_for_tilt_transform_at(4L, 0.6, 0.3, 1, 1), 0L)
  # The series itself stops where the terms do not fit
  expect_error(
    primarycensored_exptilt_lcdf(1e6, 1L, params, 1, -1),
    "needs more than 20000 terms"
  )
  # and converges just inside the limit
  for (d in c(1.7e4, 1.8e4)) {
    expect_true(is.finite(
      primarycensored_exptilt_lcdf(d, 1L, c(0, 1), 1, -1)
    ))
  }
})

test_that("the lognormal uses the ODE path past the series cut-off", {
  params <- c(log(300), 0.5)
  case <- list(meanlog = log(300), sdlog = 0.5)
  rho <- -1
  d <- c(30, 59, 61, 100, 300)
  expected <- exptilt_reference(d, 1, rho, exptilt_cdf(
    exptilt_lnorm_family(case)
  ))
  lcdf <- vapply(
    d, primarycensored_lcdf, numeric(1),
    1L, params, 1, 0, Inf, 2L, rho
  )
  # The analytical solution where the series is used
  expect_identical(lcdf[1:2], lnorm_stan_lcdf(d[1:2], case, 1, rho))
  # and the ODE, with its tolerance of 1e-6, past it
  ode <- vapply(
    d[3:5], function(x) {
      log(primarycensored_numeric_cdf(x, 1L, params, 1, 2L, rho))
    },
    numeric(1)
  )
  expect_identical(lcdf[3:5], ode)
  expect_lt(max(abs(exp(lcdf) - expected)), 1e-6)
  plain <- vapply(
    d, primarycensored_cdf, numeric(1),
    1L, params, 1, 0, Inf, 2L, rho
  )
  expect_lt(max(abs(plain - expected)), 1e-6)
  # A wide window keeps the series
  lcdf <- vapply(
    d, primarycensored_lcdf, numeric(1),
    1L, params, 3, 0, Inf, 2L, rho
  )
  expect_identical(lcdf, lnorm_stan_lcdf(d, case, 3, rho))
})

test_that("the lognormal uses the ODE path where the series is too long", {
  # r = -20 needs about 20 d terms. The ODE path has tolerances of 1e-6
  params <- c(6, 0.5)
  rho <- -20
  d <- c(200, 900, 950, 1100)
  expected <- exptilt_reference(
    d, 1, rho, function(x) plnorm(x, 6, 0.5)
  )
  lcdf <- vapply(
    d, primarycensored_lcdf, numeric(1),
    1L, params, 1, 0, Inf, 2L, rho
  )
  expect_true(all(is.finite(lcdf)))
  expect_lt(max(abs(exp(lcdf) - expected)), 1e-4)
  # Where the series fits the analytical solution is used
  expect_identical(
    lcdf[1:2], lnorm_stan_lcdf(d[1:2], list(meanlog = 6, sdlog = 0.5), 1, rho)
  )
  plain <- vapply(
    d, primarycensored_cdf, numeric(1),
    1L, params, 1, 0, Inf, 2L, rho
  )
  expect_lt(max(abs(plain - expected)), 1e-4)
  # The vectorised PMF path does the same for the delays it is given
  vectorised <- primarycensored_lcdf_vectorized(
    990L, 1000L, 1L, params, 1, 2L, rho
  )
  expect_identical(
    vectorised[990:1000], lnorm_per_delay_lcdf(990:1000, params, 1, rho)
  )
})

test_that("the vectorised lognormal tilted CDF matches the per delay CDF far
  from the origin", {
  # With |r| w below 1e-4 the small window form applies up to |r| (d + w) of
  # 0.1 and the direct form beyond, so ranges that cross it mix both forms.
  # Each endpoint has its terms computed once for the delays that need them
  pwindow <- 1
  params <- c(log(1e4), 0.5)
  for (rho in c(5e-5, -5e-5)) {
    for (range in list(c(1L, 40L), c(1990L, 2010L), c(19980L, 20000L))) {
      expect_identical(
        primarycensored_exptilt_lcdf_vectorized(
          range[1], range[2], 1L, params, pwindow, rho
        )[range[1]:range[2]],
        lnorm_per_delay_lcdf(range[1]:range[2], params, pwindow, rho),
        info = sprintf("start %g, end %g, r %g", range[1], range[2], rho)
      )
    }
  }
  # Mixed small delay, small window and direct forms
  params <- c(1.6, 0.5)
  for (rho in c(3e-5, -3e-5)) {
    expect_identical(
      primarycensored_exptilt_lcdf_vectorized(
        1L, 5000L, 1L, params, 2, rho
      )[c(1:30, 2990:3010, 4990:5000)],
      lnorm_per_delay_lcdf(
        c(1:30, 2990:3010, 4990:5000), params, 2, rho
      )
    )
  }
})

test_that("the lognormal uses the ODE path where the tilt overflows", {
  params <- c(650, 1)
  expect_identical(check_for_analytical_params(1L, params, 2L, 1e300), 0L)
  expect_error(
    primarycensored_analytical_lcdf(2, 1L, params, 2, 0, Inf, 2L, 1e300),
    "tilted delay distribution"
  )
})

test_that("the lognormal uniform primary solution is unchanged", {
  d <- c(0.5, 2, 6)
  for (case in lnorm_stan_cases[1:3]) {
    params <- lnorm_params(case)
    lcdf <- vapply(
      d, primarycensored_lcdf, numeric(1),
      1L, params, 2, 0, Inf, 1L, numeric(0)
    )
    ode <- vapply(
      d, primarycensored_numeric_cdf, numeric(1),
      1L, params, 2, 1L, numeric(0)
    )
    expect_lt(max(abs(exp(lcdf) - ode)), 1e-5)
  }
})

test_that("the vectorised lognormal tilted CDF matches the per delay CDF", {
  n <- 31L
  for (case in lnorm_stan_cases) {
    params <- lnorm_params(case)
    for (pwindow in c(1, 2, 7)) {
      # Includes the small window form, the small delay form for the first
      # delays, the direct form and both signs of the tilt
      for (rho in c(-0.3, -2e-5, -1e-9, 0, 1e-9, 2e-5, 1e-5, 0.4)) {
        for (start in c(1L, 5L)) {
          vectorised <- primarycensored_exptilt_lcdf_vectorized(
            start, n, 1L, params, pwindow, rho
          )
          expect_length(vectorised, n)
          expect_identical(
            vectorised[start:n],
            lnorm_per_delay_lcdf(start:n, params, pwindow, rho),
            info = sprintf(
              "meanlog %g, sdlog %g, pwindow %g, r %g, start %g",
              case$meanlog, case$sdlog, pwindow, rho, start
            )
          )
        }
      }
    }
  }
})

test_that("primarycensored_lcdf_vectorized uses the lognormal shared terms", {
  params <- c(1.6, 0.5)
  expect_identical(
    primarycensored_lcdf_vectorized(1L, 20L, 1L, params, 3, 2L, 0.25),
    primarycensored_exptilt_lcdf_vectorized(1L, 20L, 1L, params, 3, 0.25)
  )
  expect_identical(
    primarycensored_lcdf_vectorized(1L, 10L, 1L, params, 1.5, 2L, 0.2),
    lnorm_per_delay_lcdf(1:10, params, 1.5, 0.2)
  )
})

test_that("the vectorised lognormal PMF matches the per delay PMF with
  truncation", {
  settings <- list(
    list(max_delay = 10, L = 0, D = 11),
    list(max_delay = 10, L = 0, D = Inf),
    list(max_delay = 20, L = 2, D = 21),
    list(max_delay = 20, L = 2, D = 30)
  )
  for (case in lnorm_stan_cases[1:4]) {
    params <- lnorm_params(case)
    for (setting in settings) {
      for (pwindow in c(1, 3)) {
        for (rho in c(-0.2, 1e-9, 0.3)) {
          max_delay <- setting$max_delay
          vectorised <- primarycensored_sone_lpmf_vectorized(
            max_delay, setting$L, setting$D, 1L, params, pwindow, 2L, rho
          )
          per_delay <- vapply(
            0:max_delay, function(d) {
              primarycensored_lpmf(
                d, 1L, params, pwindow, d + 1, setting$L, setting$D, 2L, rho
              )
            },
            numeric(1)
          )
          expect_equal(
            vectorised, per_delay,
            tolerance = 1e-10,
            info = sprintf(
              "meanlog %g, sdlog %g, pwindow %g, r %g, L %g, D %g",
              case$meanlog, case$sdlog, pwindow, rho, setting$L, setting$D
            )
          )
        }
      }
    }
  }
})

test_that("the vectorised lognormal PMF matches differences of the
  reference CDF", {
  for (case in lnorm_stan_cases[1:3]) {
    params <- lnorm_params(case)
    cdf <- exptilt_lnorm_cdf(case)
    for (rho in c(-0.4, 0.4)) {
      pmf <- primarycensored_sone_pmf_vectorized(
        12, 0, Inf, 1L, params, 3, 2L, rho
      )
      expected <- diff(exptilt_reference(0:13, 3, rho, cdf))
      expect_equal(pmf, expected, tolerance = 1e-8)
    }
  }
})

# Gradients are only observable from a compiled model, so this builds a
# minimal one whose target is the log CDF or the vectorised log PMF of the
# lognormal, and runs `stan_gradient_at()` from helper-stan-gradient.R.
lnorm_gradient_model <- function() {
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
    "  int vectorised;\n",
    "  real d;\n",
    "  real pwindow;\n",
    "}\n",
    "parameters {\n",
    "  real meanlog;\n",
    "  real<lower=0> sdlog;\n",
    "  real rho;\n",
    "}\n",
    "model {\n",
    "  array[2] real params = {meanlog, sdlog};\n",
    "  if (vectorised) {\n",
    "    target += sum(primarycensored_sone_lpmf_vectorized(\n",
    "      to_int(d), 0, positive_infinity(), 1, params, pwindow, 2, {rho}\n",
    "    ));\n",
    "  } else {\n",
    "    target += primarycensored_lcdf(\n",
    "      d | 1, params, pwindow, 0, positive_infinity(), 2, {rho}\n",
    "    );\n",
    "  }\n",
    "}\n"
  )
  path <- file.path(tempdir(), "pcd_lnorm_exptilt_gradient.stan")
  writeLines(code, path)
  suppressMessages(suppressWarnings(cmdstanr::cmdstan_model(path)))
}

lnorm_gradient_at <- function(model, case, d, pwindow, rho,
                              vectorised = FALSE) {
  stan_gradient_at( # nolint: object_usage_linter.
    model,
    data = list(
      vectorised = as.integer(vectorised), d = d, pwindow = pwindow
    ),
    init = list(meanlog = case$meanlog, sdlog = case$sdlog, rho = rho)
  )
}

# Compares the gradient with the finite difference gradient one component at a
# time, relative to the size of the component with a floor for tiny ones.
expect_lnorm_gradient_close <- function(res, label, scale = 1) {
  allowed <- 1e-4 * scale * pmax(abs(res$finite_diff), 1e-2)
  testthat::expect_true(
    all(abs(res$gradient - res$finite_diff) <= allowed),
    info = paste0(
      label, ": gradient ", toString(signif(res$gradient, 5)),
      ", finite difference ", toString(signif(res$finite_diff, 5))
    )
  )
}

test_that("lognormal tilted log CDFs have finite gradients matching finite
  differences", {
  model <- lnorm_gradient_model()
  # Direct, small window and small delay forms; both tails; both signs of the
  # tilt; d below and above pwindow. The points are not at a switch between
  # forms, as finite differences would step across it.
  points <- list(
    list(d = 0.3, pwindow = 2, rho = 0.4),
    list(d = 1, pwindow = 1, rho = -0.2),
    list(d = 2.5, pwindow = 2, rho = 0.4),
    list(d = 6, pwindow = 3, rho = -0.15),
    list(d = 20, pwindow = 3, rho = 0.5),
    list(d = 20, pwindow = 3, rho = -0.5),
    list(d = 2.5, pwindow = 2, rho = 1e-6),
    list(d = 2.5, pwindow = 2, rho = -1e-6),
    list(d = 12, pwindow = 7, rho = 1e-5),
    list(d = 0.0001, pwindow = 2, rho = 0.4),
    list(d = 0.0001, pwindow = 2, rho = -0.4)
  )
  for (case in lnorm_stan_cases[c(1, 2, 4)]) {
    for (point in points) {
      label <- sprintf(
        "meanlog %g, sdlog %g, d %g, pwindow %g, r %g",
        case$meanlog, case$sdlog, point$d, point$pwindow, point$rho
      )
      res <- lnorm_gradient_at(model, case, point$d, point$pwindow, point$rho)
      expect_false(res$gradient_not_finite, info = label)
      expect_false(res$rejected, info = label)
      expect_length(res$gradient, 3)
      expect_true(all(is.finite(res$gradient)), info = label)
      expect_lnorm_gradient_close(res, label)
    }
  }
})

test_that("the vectorised lognormal tilted log PMF has finite gradients
  matching finite differences", {
  model <- lnorm_gradient_model()
  points <- list(
    list(d = 12, pwindow = 3, rho = 0.4),
    list(d = 12, pwindow = 3, rho = -0.15),
    list(d = 6, pwindow = 2, rho = 1e-6)
  )
  for (case in lnorm_stan_cases[c(1, 2, 4)]) {
    for (point in points) {
      label <- sprintf(
        "meanlog %g, sdlog %g, d %g, pwindow %g, r %g",
        case$meanlog, case$sdlog, point$d, point$pwindow, point$rho
      )
      res <- lnorm_gradient_at(
        model, case, point$d, point$pwindow, point$rho,
        vectorised = TRUE
      )
      expect_false(res$gradient_not_finite, info = label)
      expect_false(res$rejected, info = label)
      expect_true(all(is.finite(res$gradient)), info = label)
      expect_lnorm_gradient_close(res, label)
    }
  }
})

test_that("the lognormal tilted CDF is accurate for a large sdlog", {
  # One panel per side lost accuracy for sdlog above about 2. Panels are
  # split in proportion to sdlog
  for (sdlog in c(4, 6, 10, 15)) {
    for (meanlog in c(0, 4)) {
      case <- list(meanlog = meanlog, sdlog = sdlog)
      d <- exp(meanlog + sdlog * seq(-2, 0.5, by = 0.25))
      d <- d[d > 1e-9]
      # The series for a negative tilt needs xi d of at most about 2e4
      for (rho in c(5, 200, -0.005)) {
        expected <- exptilt_reference(
          d, 1, rho, exptilt_lnorm_cdf(case)
        )
        keep <- expected > 1e-100
        actual <- exp(lnorm_stan_lcdf(d, case, 1, rho))
        expect_lt(
          max_rel_diff(actual[keep], expected[keep]), 1e-8,
          label = sprintf(
            "sdlog %g, meanlog %g, r %g", sdlog, meanlog, rho
          )
        )
      }
    }
  }
})

test_that("lognormal tilted log CDFs have finite gradients for a large
  sdlog", {
  model <- lnorm_gradient_model()
  case <- list(meanlog = 2, sdlog = 6)
  points <- list(
    list(d = 0.05, pwindow = 1, rho = 5),
    list(d = 3, pwindow = 2, rho = 0.4),
    list(d = 20, pwindow = 3, rho = -0.3)
  )
  for (point in points) {
    label <- sprintf(
      "meanlog %g, sdlog %g, d %g, pwindow %g, r %g",
      case$meanlog, case$sdlog, point$d, point$pwindow, point$rho
    )
    res <- lnorm_gradient_at(model, case, point$d, point$pwindow, point$rho)
    expect_false(res$gradient_not_finite, info = label)
    expect_false(res$rejected, info = label)
    expect_true(all(is.finite(res$gradient)), info = label)
    expect_lnorm_gradient_close(res, label)
  }
})
