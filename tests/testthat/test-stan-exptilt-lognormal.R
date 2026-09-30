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
})

test_that("a single lognormal delay uses the ODE, as it is faster than the
  quadrature", {
  # The quadrature is faster than the ODE only with the shared terms over
  # integer delays, so the scalar functions keep the ODE.
  for (rho in c(-2, -0.3, 0.3, 2)) {
    expect_identical(
      check_for_analytical_params(1L, c(1.6, 0.5), 2L, rho), 0L
    )
    for (d in c(0.5, 2, 5, 10)) {
      ode <- primarycensored_numeric_cdf(d, 1L, c(1.6, 0.5), 2, 2L, rho)
      expect_identical(
        primarycensored_cdf(d, 1L, c(1.6, 0.5), 2, 0, Inf, 2L, rho), ode
      )
      expect_identical(
        primarycensored_lcdf(d, 1L, c(1.6, 0.5), 2, 0, Inf, 2L, rho), log(ode)
      )
    }
  }
})

test_that("Stan lognormal tilt transforms match the R transforms", {
  z <- c(-12, -6, -3, -1, 0, 0.5, 1, 2, 3)
  for (case in lnorm_stan_cases) {
    obj <- lnorm_object_stan(case, 0.1)
    t <- exp(case$meanlog + case$sdlog * z)
    params <- lnorm_params(case)
    for (xi in c(-5, -1, -0.25, -0.01, 0, 1e-3, 0.1, 0.5, 1)) {
      info <- sprintf(
        "meanlog %g, sdlog %g, xi %g", case$meanlog, case$sdlog, xi
      )
      lower <- vapply(t, log_tilt_transform, numeric(1), 1L, xi, params)
      upper <- vapply(t, log_tilt_transform_upper, numeric(1), 1L, xi, params)
      expected <- .pcens_tilt_pair(obj, t, xi)
      keep <- expected[, 1] > -700
      expect_lt(
        max(abs(lower[keep] - expected[keep, 1])), 1e-9,
        label = info
      )
      if (xi > 0) {
        expect_identical(upper, rep(Inf, length(t)), info = info)
      } else {
        keep <- expected[, 2] > -700
        expect_lt(
          max(abs(upper[keep] - expected[keep, 2])), 1e-9,
          label = info
        )
      }
    }
  }
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

test_that("the lognormal analytical function matches the reference and the
  ODE path", {
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
          d, primarycensored_analytical_lcdf, numeric(1),
          1L, params, pwindow, 0, Inf, 2L, rho
        )
        expect_lt(max_rel_diff(exp(lcdf), expected), 1e-7, label = info)
        # The ODE path has absolute and relative tolerances of 1e-6
        ode <- vapply(
          d, primarycensored_numeric_cdf, numeric(1),
          1L, params, pwindow, 2L, rho
        )
        expect_lt(max(abs(exp(lcdf) - ode)), 1e-4, label = info)
      }
    }
  }
})

test_that("the lognormal analytical function rejects a tilt that overflows", {
  params <- c(650, 1)
  expect_identical(check_for_tilt_transform(1L, -1e300, params), 0L)
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

lnorm_per_delay_lcdf <- function(delays, params, pwindow, rho) {
  # nolint start: object_usage_linter.
  vapply(
    delays, primarycensored_exptilt_lcdf, numeric(1),
    1L, params, pwindow, rho
  )
  # nolint end
}

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
  # A non-integer window uses the scalar function, which is the ODE
  expect_identical(
    primarycensored_lcdf_vectorized(1L, 10L, 1L, params, 1.5, 2L, 0.2),
    vapply(
      1:10, primarycensored_lcdf, numeric(1), # nolint: object_usage_linter.
      1L, params, 1.5, 0, Inf, 2L, 0.2
    )
  )
  # A tilt that overflows uses the scalar function too
  expect_identical(
    primarycensored_lcdf_vectorized(1L, 5L, 1L, c(650, 1), 3, 2L, 1e300),
    vapply(
      1:5, primarycensored_lcdf, numeric(1), # nolint: object_usage_linter.
      1L, c(650, 1), 3, 0, Inf, 2L, 1e300
    )
  )
})

test_that("the vectorised lognormal PMF matches the reference with
  truncation", {
  settings <- list(
    list(max_delay = 10, L = 0, D = 11),
    list(max_delay = 10, L = 0, D = Inf),
    list(max_delay = 20, L = 2, D = 21),
    list(max_delay = 20, L = 2, D = 30)
  )
  for (case in lnorm_stan_cases[1:4]) {
    params <- lnorm_params(case)
    cdf <- exptilt_lnorm_cdf(case)
    for (setting in settings) {
      for (pwindow in c(1, 3)) {
        for (rho in c(-0.2, 1e-9, 0.3)) {
          max_delay <- setting$max_delay
          vectorised <- primarycensored_sone_lpmf_vectorized(
            max_delay, setting$L, setting$D, 1L, params, pwindow, 2L, rho
          )
          ref <- function(x) exptilt_reference(x, pwindow, rho, cdf)
          cdf_lower <- if (setting$L > 0) ref(setting$L) else 0
          cdf_upper <- if (is.finite(setting$D)) ref(setting$D) else 1
          delays <- 0:max_delay
          pmf <- (ref(delays + 1) - ref(delays)) / (cdf_upper - cdf_lower)
          pmf[delays < setting$L] <- 0
          # Differences of the reference CDF are accurate to about 1e-16
          keep <- pmf > 1e-8
          # A finite D beyond max_delay + 1 is normalised with the scalar
          # CDF, which is the ODE with a tolerance of 1e-6
          beyond <- is.finite(setting$D) && setting$D > max_delay + 1
          expect_equal(
            exp(vectorised)[keep], pmf[keep],
            tolerance = if (beyond) 1e-5 else 1e-7,
            info = sprintf(
              "meanlog %g, sdlog %g, pwindow %g, r %g, L %g, D %g",
              case$meanlog, case$sdlog, pwindow, rho, setting$L, setting$D
            )
          )
          expect_true(all(is.infinite(vectorised[delays < setting$L])))
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
# lognormal, and runs `stan_gradient_at()` from helper-stan-gradient.R. The
# scalar function is `primarycensored_exptilt_lcdf()`, as
# `primarycensored_lcdf()` uses the ODE for the lognormal.
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
    "    target += primarycensored_exptilt_lcdf(\n",
    "      d | 1, params, pwindow, rho\n",
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
