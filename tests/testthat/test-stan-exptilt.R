skip_on_cran()

# Stan solutions for exponential, gamma and normal delays with an
# exponentially tilted primary (primary_id 2 with primary_params = r). These
# tests check the Stan transforms and CDFs against the R implementation and a
# reference integral, the dispatch and fallback to the ODE path, the shared
# endpoint vectorised form, and gradients.

exptilt_stan_cases <- list(
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
    dist_id = 18L, params = c(3, 2), pdist = pnorm,
    args = list(mean = 3, sd = 2)
  ),
  list(
    dist_id = 18L, params = c(-1, 3), pdist = pnorm,
    args = list(mean = -1, sd = 3)
  )
)

exptilt_case_rate <- function(case) {
  switch(as.character(case$dist_id),
    "4" = case$params[1],
    "2" = case$params[2],
    Inf
  )
}

exptilt_case_label <- function(case, ...) {
  paste0(
    "dist ", case$dist_id, " params ", toString(case$params), ", ",
    paste(names(list(...)), unlist(list(...)), sep = " = ", collapse = ", ")
  )
}

exptilt_case_cdf <- function(case) {
  function(x) do.call(case$pdist, c(list(x), case$args))
}

# Internal lower bound used by primarycensored_lcdf for each support
exptilt_case_lower <- function(case) {
  if (case$dist_id == 18L) -Inf else 0
}

test_that("check_for_tilt_transform needs the tilted delay to exist", {
  expect_identical(check_for_tilt_transform(4L, -0.5, 0.3), 1L)
  expect_identical(check_for_tilt_transform(4L, 0.29, 0.3), 1L)
  expect_identical(check_for_tilt_transform(4L, 0.3, 0.3), 0L)
  expect_identical(check_for_tilt_transform(4L, 1, 0.3), 0L)
  expect_identical(check_for_tilt_transform(2L, -0.3, c(2, 0.4)), 1L)
  expect_identical(check_for_tilt_transform(2L, 0.39, c(2, 0.4)), 1L)
  expect_identical(check_for_tilt_transform(2L, 0.4, c(2, 0.4)), 0L)
  expect_identical(check_for_tilt_transform(2L, 0, c(2, 0.4)), 1L)
  expect_identical(check_for_tilt_transform(18L, -50, c(3, 2)), 1L)
  expect_identical(check_for_tilt_transform(18L, 50, c(3, 2)), 1L)
  for (dist_id in c(1L, 3L, 5L, 12L, 26L, 27L)) {
    expect_identical(check_for_tilt_transform(dist_id, 0, c(1, 1)), 0L)
  }
})

test_that("check_for_analytical includes the exponentially tilted primary", {
  for (dist_id in c(2L, 4L, 18L)) {
    expect_identical(check_for_analytical(dist_id, 2L), 1L)
    expect_identical(check_for_exptilt(dist_id, 2L), 1L)
    expect_identical(check_for_exptilt(dist_id, 1L), 0L)
  }
  # Other delays stay numerical with an exponentially tilted primary
  for (dist_id in c(1L, 3L, 5L, 9L)) {
    expect_identical(check_for_analytical(dist_id, 2L), 0L)
    expect_identical(check_for_exptilt(dist_id, 2L), 0L)
  }
  # The uniform terms are unchanged
  expect_identical(check_for_analytical(4L, 1L), 0L)
  expect_identical(check_for_analytical(18L, 1L), 0L)
  expect_identical(check_for_analytical(2L, 1L), 1L)
  expect_identical(check_for_uniform_terms(2L, 2L), 0L)
})

test_that("check_for_analytical_params adds the admissibility of the tilt", {
  # tilt r = 0.5 needs rate + r > 0
  expect_identical(
    check_for_analytical_params(2L, c(2, 0.4), 2L, 0.5), 1L
  )
  expect_identical(
    check_for_analytical_params(2L, c(2, 0.4), 2L, -0.3), 1L
  )
  expect_identical(
    check_for_analytical_params(2L, c(2, 0.4), 2L, -0.4), 0L
  )
  expect_identical(
    check_for_analytical_params(4L, 0.3, 2L, -0.5), 0L
  )
  expect_identical(
    check_for_analytical_params(4L, 0.3, 2L, -0.25), 1L
  )
  expect_identical(
    check_for_analytical_params(18L, c(3, 2), 2L, -10), 1L
  )
  # Unchanged for the other solutions
  expect_identical(
    check_for_analytical_params(2L, c(2, 0.4), 1L, numeric(0)), 1L
  )
  expect_identical(
    check_for_analytical_params(3L, c(2, 1), 2L, 0.3), 0L
  )
  expect_identical(
    check_for_analytical_params(26L, c(0, 1, 2, 0.5, 0.5), 2L, 0.3), 1L
  )
})

test_that("Stan tilt transforms match the R transforms", {
  ts <- c(1e-3, 0.4, 1, 3.5, 12, 40)
  for (case in exptilt_stan_cases) {
    obj <- exptilt_object(
      list(pdist = case$pdist, args = case$args), 0.1
    )
    for (xi in c(0, -1, -0.25, 0.1, 0.25)) {
      if (!check_for_tilt_transform(case$dist_id, xi, case$params)) {
        next
      }
      lower <- vapply(
        ts, log_tilt_transform, numeric(1), case$dist_id, xi, case$params
      )
      upper <- vapply(
        ts, log_tilt_transform_upper, numeric(1), case$dist_id, xi,
        case$params
      )
      info <- exptilt_case_label(case, xi = xi)
      expect_equal(lower, .pcens_tilt_transform(obj, ts, xi), info = info)
      expect_equal(
        upper, .pcens_tilt_transform(obj, ts, xi, upper = TRUE),
        info = info
      )
    }
  }
})

test_that("Stan gamma tilt transforms are accurate in both tails", {
  # One tail is evaluated and the other follows from it, so this covers the
  # switch to the upper tail at 1e-8 for small and large shapes. Stan returns
  # -Inf where a term is below the smallest double, where R gives the log.
  ts <- 10^seq(-8, 4, by = 0.5)
  for (shape in c(0.05, 0.3, 1, 7, 100, 1000)) {
    obj <- new_pcens(
      pgamma, dexpgrowth, list(r = 0.1),
      shape = shape, rate = 1.7
    )
    params <- c(shape, 1.7)
    for (xi in c(0, -0.5, 0.9)) {
      for (upper in c(FALSE, TRUE)) {
        stan_fun <- if (upper) log_tilt_transform_upper else log_tilt_transform
        actual <- vapply(
          ts, stan_fun, numeric(1), 2L, xi, params
        )
        expected <- .pcens_tilt_transform(obj, ts, xi, upper = upper)
        # The total is always defined so compare the values that Stan keeps
        keep <- is.finite(actual)
        expect_equal(
          actual[keep], expected[keep],
          tolerance = 1e-9,
          info = sprintf("shape %g, xi %g, upper %s", shape, xi, upper)
        )
        # The transforms of a gamma with a large total can be representable
        # a little below the smallest probability
        expect_true(
          all(expected[!keep] < -600),
          info = sprintf("shape %g, xi %g, upper %s", shape, xi, upper)
        )
      }
    }
  }
})

test_that("Stan tilt transforms are 0 or the total below the support", {
  for (case in exptilt_stan_cases[c(2, 4)]) {
    for (t in c(-3, -1e-9, 0)) {
      expect_identical(
        log_tilt_transform(t, case$dist_id, -0.1, case$params), -Inf
      )
    }
    total <- exp(log_tilt_transform_upper(0, case$dist_id, -0.1, case$params))
    expect_identical(
      exp(log_tilt_transform_upper(-2, case$dist_id, -0.1, case$params)),
      total
    )
  }
  expect_error(
    log_tilt_transform(1, 3L, 0, c(1, 1)),
    "Invalid distribution identifier"
  )
  expect_error(
    log_tilt_transform_upper(1, 3L, 0, c(1, 1)),
    "Invalid distribution identifier"
  )
  expect_error(
    primarycensored_tilt_moments(1, 3L, c(1, 1)),
    "Invalid distribution identifier"
  )
})

test_that("primarycensored_exptilt_lcdf matches a reference integral", {
  for (case in exptilt_stan_cases) {
    cdf <- exptilt_case_cdf(case)
    for (pwindow in c(0.5, 1, 2, 7)) {
      d <- sort(c(
        1e-6, 1e-3, 0.3 * pwindow, pwindow - 1e-3, pwindow, pwindow + 1e-3,
        2, 3, 6, 12, 25,
        if (case$dist_id == 18L) c(-10, -3, -0.5)
      ))
      for (rho in c(
        -1, -0.5, -0.05, -1e-4, -1e-5, -1e-8, 1e-8, 1e-5, 1e-4,
        0.05, 0.5, 1
      )) {
        if (case$dist_id != 18L &&
          exptilt_case_rate(case) + rho <= 0) {
          next
        }
        expected <- exptilt_reference(d, pwindow, rho, cdf)
        actual <- exp(vapply(
          d, primarycensored_exptilt_lcdf, numeric(1),
          case$dist_id, case$params, pwindow, rho
        ))
        expect_lt(
          max_rel_diff(actual, expected), 1e-7,
          label = exptilt_case_label(case, pwindow = pwindow, r = rho)
        )
      }
    }
  }
})

test_that("primarycensored_exptilt_lcdf matches the R implementation", {
  d <- c(1e-4, 0.3, 1, 2.5, 6, 15, 30)
  for (case in exptilt_stan_cases) {
    for (pwindow in c(0.5, 2, 7)) {
      for (rho in c(-0.2, -1e-6, 1e-6, 0.3)) {
        if (case$dist_id != 18L &&
          exptilt_case_rate(case) + rho <= 0) {
          next
        }
        obj <- exptilt_object(
          list(pdist = case$pdist, args = case$args), rho
        )
        expect_equal(
          exp(vapply(
            d, primarycensored_exptilt_lcdf, numeric(1),
            case$dist_id, case$params, pwindow, rho
          )),
          pcens_cdf(obj, d, pwindow),
          tolerance = 1e-9,
          info = exptilt_case_label(case, pwindow = pwindow, r = rho)
        )
      }
    }
  }
})

test_that("the tilted CDF is continuous across the small tilt forms", {
  d <- c(1e-4, 0.3, 1, 2.5, 6, 15, 30)
  for (case in exptilt_stan_cases) {
    for (pwindow in c(0.5, 2, 7)) {
      for (sign in c(-1, 1)) {
        expect_lt(
          max_rel_diff(
            exp(vapply(
              d, primarycensored_exptilt_lcdf, numeric(1),
              case$dist_id, case$params, pwindow, sign * 0.9999e-4 / pwindow
            )),
            exp(vapply(
              d, primarycensored_exptilt_lcdf, numeric(1),
              case$dist_id, case$params, pwindow, sign * 1.0001e-4 / pwindow
            ))
          ),
          1e-7,
          label = exptilt_case_label(case, pwindow = pwindow, sign = sign)
        )
      }
    }
  }
})

test_that("primarycensored_lcdf and primarycensored_cdf use the analytical
  solution and agree with the ODE path", {
  d <- c(0.2, 1, 2.5, 6, 15)
  for (case in exptilt_stan_cases) {
    lower <- exptilt_case_lower(case)
    cdf <- exptilt_case_cdf(case)
    for (pwindow in c(1, 3)) {
      for (rho in c(-1, -0.5, -1e-8, 1e-8, 0.5, 1)) {
        if (case$dist_id != 18L &&
          exptilt_case_rate(case) + rho <= 0) {
          next
        }
        info <- exptilt_case_label(case, pwindow = pwindow, r = rho)
        expect_identical(
          check_for_analytical_params(case$dist_id, case$params, 2L, rho), 1L
        )
        expected <- exptilt_reference(d, pwindow, rho, cdf)
        lcdf <- vapply(
          d, primarycensored_lcdf, numeric(1),
          case$dist_id, case$params, pwindow, lower, Inf, 2L, rho
        )
        expect_lt(max_rel_diff(exp(lcdf), expected), 1e-7, label = info)
        plain <- vapply(
          d, primarycensored_cdf, numeric(1),
          case$dist_id, case$params, pwindow, lower, Inf, 2L, rho
        )
        expect_lt(max_rel_diff(plain, expected), 1e-7, label = info)
        # The ODE path has absolute and relative tolerances of 1e-6, and
        # about 1e-4 for a shape below 1 where the density is singular at 0
        ode <- vapply(
          d, primarycensored_numeric_cdf, numeric(1),
          case$dist_id, case$params, pwindow, 2L, rho
        )
        expect_lt(max(abs(plain - ode)), 1e-4, label = info)
      }
    }
  }
})

test_that("inadmissible tilts use the ODE path", {
  cases <- list(
    list(dist_id = 4L, params = 0.3, rho = -0.5),
    list(dist_id = 4L, params = 0.3, rho = -0.3),
    list(dist_id = 2L, params = c(2.5, 0.4), rho = -0.5),
    list(dist_id = 2L, params = c(2.5, 0.4), rho = -0.4)
  )
  for (case in cases) {
    expect_identical(
      check_for_analytical_params(
        case$dist_id, case$params, 2L, case$rho
      ), 0L
    )
    for (d in c(0.5, 2, 5, 10)) {
      ode <- primarycensored_numeric_cdf(
        d, case$dist_id, case$params, 2, 2L, case$rho
      )
      expect_identical(
        primarycensored_cdf(
          d, case$dist_id, case$params, 2, 0, Inf, 2L, case$rho
        ),
        ode
      )
      expect_identical(
        primarycensored_lcdf(
          d, case$dist_id, case$params, 2, 0, Inf, 2L, case$rho
        ),
        log(ode)
      )
    }
  }
})

test_that("the analytical function rejects an inadmissible tilt", {
  expect_error(
    primarycensored_analytical_lcdf(
      2, 2L, c(2.5, 0.4), 2, 0, Inf, 2L, -0.5
    ),
    "tilted delay distribution"
  )
})

test_that("primarycensored_numeric_cdf is the CDF for other delays", {
  # Not analytical with an exponentially tilted primary, so the public CDF
  # is the ODE result
  for (d in c(0.5, 2, 5)) {
    expect_identical(
      primarycensored_cdf(d, 3L, c(1.5, 2), 2, 0, Inf, 2L, 0.3),
      primarycensored_numeric_cdf(d, 3L, c(1.5, 2), 2, 2L, 0.3)
    )
  }
})

test_that("normal delays handle negative delays and truncation", {
  pwindow <- 2
  rho <- 0.3
  cdf <- function(x) pnorm(x, 3, 2)
  ref <- function(x) exptilt_reference(x, pwindow, rho, cdf)
  d <- c(-6, -2, -0.5, 0, 1, 3, 6)
  lcdf <- vapply(
    d, primarycensored_lcdf, numeric(1),
    18L, c(3, 2), pwindow, -Inf, Inf, 2L, rho
  )
  expect_lt(max_rel_diff(exp(lcdf), ref(d)), 1e-7)
  for (bounds in list(c(-2, 9), c(-Inf, 9), c(-2, Inf), c(0.5, 7))) {
    L <- bounds[1]
    D <- bounds[2]
    lower <- if (is.finite(L)) ref(L) else 0
    upper <- if (is.finite(D)) ref(D) else 1
    x <- c(1, 3, 5)
    expected <- (ref(x) - lower) / (upper - lower)
    actual <- vapply(
      x, primarycensored_cdf, numeric(1),
      18L, c(3, 2), pwindow, L, D, 2L, rho
    )
    expect_equal(actual, expected, tolerance = 1e-7)
    actual_l <- vapply(
      x, primarycensored_lcdf, numeric(1),
      18L, c(3, 2), pwindow, L, D, 2L, rho
    )
    expect_equal(exp(actual_l), expected, tolerance = 1e-7)
  }
})

per_delay_exptilt_lcdf <- function(delays, dist_id, params, pwindow, rho) {
  lower <- if (dist_id == 18L) -Inf else 0
  vapply(
    delays, primarycensored_lcdf, numeric(1), # nolint: object_usage_linter.
    dist_id, params, pwindow, lower, Inf, 2L, rho
  )
}

test_that("check_for_exptilt_vectorized needs an integer pwindow", {
  for (dist_id in c(2L, 4L, 18L)) {
    expect_identical(check_for_exptilt_vectorized(dist_id, 2L, 1), 1L)
    expect_identical(check_for_exptilt_vectorized(dist_id, 2L, 7), 1L)
    expect_identical(check_for_exptilt_vectorized(dist_id, 2L, 1.5), 0L)
    expect_identical(check_for_exptilt_vectorized(dist_id, 2L, 0.5), 0L)
    expect_identical(check_for_exptilt_vectorized(dist_id, 1L, 1), 0L)
  }
  for (dist_id in c(1L, 3L, 26L)) {
    expect_identical(check_for_exptilt_vectorized(dist_id, 2L, 1), 0L)
  }
})

test_that("the vectorised tilted CDF matches the per delay CDF", {
  n <- 31L
  for (case in exptilt_stan_cases) {
    for (pwindow in c(1, 2, 7)) {
      # Includes the small window form (r * pwindow < 1e-4), the small delay
      # form for the first delays, and the direct form
      for (rho in c(-0.3, -2e-5, -1e-9, 0, 1e-9, 2e-5, 1e-5, 0.4)) {
        if (case$dist_id != 18L &&
          exptilt_case_rate(case) + rho <= 0) {
          next
        }
        for (start in c(1L, 5L)) {
          vectorised <- primarycensored_exptilt_lcdf_vectorized(
            start, n, case$dist_id, case$params, pwindow, rho
          )
          expect_length(vectorised, n)
          expect_identical(
            vectorised[start:n],
            per_delay_exptilt_lcdf(
              start:n, case$dist_id, case$params, pwindow, rho
            ),
            info = exptilt_case_label(
              case,
              pwindow = pwindow, r = rho, start = start
            )
          )
        }
      }
    }
  }
})

test_that("the vectorised tilted CDF mixes the small delay and direct forms", {
  # r * pwindow = 2e-4 is above the small window threshold but r * d is below
  # the small delay threshold for d < 5
  pwindow <- 10
  rho <- 2e-5
  for (case in exptilt_stan_cases[c(2, 4, 5)]) {
    expect_identical(
      primarycensored_exptilt_lcdf_vectorized(
        1L, 25L, case$dist_id, case$params, pwindow, rho
      )[1:25],
      per_delay_exptilt_lcdf(1:25, case$dist_id, case$params, pwindow, rho)
    )
  }
})

test_that("primarycensored_lcdf_vectorized uses the tilted shared terms", {
  for (case in exptilt_stan_cases) {
    rho <- 0.25
    expect_identical(
      primarycensored_lcdf_vectorized(
        1L, 20L, case$dist_id, case$params, 3, 2L, rho
      ),
      primarycensored_exptilt_lcdf_vectorized(
        1L, 20L, case$dist_id, case$params, 3, rho
      )
    )
  }
  # An inadmissible tilt and a non-integer window use the per delay path
  expect_identical(
    primarycensored_lcdf_vectorized(1L, 10L, 4L, 0.3, 3, 2L, -0.5),
    per_delay_exptilt_lcdf(1:10, 4L, 0.3, 3, -0.5)
  )
  expect_identical(
    primarycensored_lcdf_vectorized(1L, 10L, 4L, 0.3, 1.5, 2L, 0.2),
    per_delay_exptilt_lcdf(1:10, 4L, 0.3, 1.5, 0.2)
  )
})

test_that("the vectorised PMF matches the per delay PMF with truncation", {
  settings <- list(
    list(max_delay = 10, L = 0, D = 11),
    list(max_delay = 10, L = 0, D = Inf),
    list(max_delay = 20, L = 2, D = 21),
    list(max_delay = 20, L = 2, D = 30),
    list(max_delay = 20, L = -Inf, D = Inf)
  )
  for (case in exptilt_stan_cases) {
    lower_support <- case$dist_id == 18L
    for (setting in settings) {
      for (pwindow in c(1, 3)) {
        for (rho in c(-0.2, 1e-9, 0.3)) {
          if (!lower_support && exptilt_case_rate(case) + rho <= 0) {
            next
          }
          max_delay <- setting$max_delay
          vectorised <- primarycensored_sone_lpmf_vectorized(
            max_delay, setting$L, setting$D, case$dist_id, case$params,
            pwindow, 2L, rho
          )
          per_delay <- vapply(
            0:max_delay, function(d) {
              primarycensored_lpmf(
                d, case$dist_id, case$params, pwindow, d + 1,
                setting$L, setting$D, 2L, rho
              )
            },
            numeric(1)
          )
          expect_equal(
            vectorised, per_delay,
            tolerance = 1e-10,
            info = exptilt_case_label(
              case,
              pwindow = pwindow, r = rho, L = setting$L, D = setting$D
            )
          )
        }
      }
    }
  }
})

test_that("the vectorised PMF of a normal delay sums to the CDF", {
  # The intervals start at 0 and so miss the mass of negative delays
  pmf <- primarycensored_sone_pmf_vectorized(
    60, -Inf, Inf, 18L, c(3, 2), 2, 2L, 0.2
  )
  cdf <- function(x) pnorm(x, 3, 2)
  expected <- diff(exptilt_reference(c(0, 61), 2, 0.2, cdf))
  expect_equal(sum(pmf), expected, tolerance = 1e-9)
})

# Gradients are only observable from a compiled model, so this builds a
# minimal one whose target is the log CDF or the vectorised log PMF, and runs
# `stan_gradient_at()` from helper-stan-gradient.R.
exptilt_gradient_model <- function() {
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
    "  real rho;\n",
    "}\n",
    "model {\n",
    "  array[2] real all_params = {p1, p2};\n",
    "  array[n_params] real params = all_params[1:n_params];\n",
    "  if (vectorised) {\n",
    "    target += sum(primarycensored_sone_lpmf_vectorized(\n",
    "      to_int(d), L, positive_infinity(), dist_id, params, pwindow, 2,\n",
    "      {rho}\n",
    "    ));\n",
    "  } else {\n",
    "    target += primarycensored_lcdf(\n",
    "      d | dist_id, params, pwindow, L, positive_infinity(), 2, {rho}\n",
    "    );\n",
    "  }\n",
    "}\n"
  )
  path <- file.path(tempdir(), "pcd_exptilt_gradient.stan")
  writeLines(code, path)
  suppressMessages(suppressWarnings(cmdstanr::cmdstan_model(path)))
}

exptilt_gradient_at <- function(model, case, d, pwindow, rho,
                                vectorised = FALSE) {
  init <- list(
    p1 = case$params[1],
    p2 = if (length(case$params) > 1) case$params[2] else 1,
    rho = rho
  )
  stan_gradient_at( # nolint: object_usage_linter.
    model,
    data = list(
      dist_id = case$dist_id, n_params = length(case$params),
      vectorised = as.integer(vectorised), d = d, pwindow = pwindow,
      L = exptilt_case_lower(case)
    ),
    init = init
  )
}

# The shape gradient of Stan's gamma_lcdf has a relative error of 1.7e-2 for
# shape 20 at 2 and of 0.5 for shape 100 at 30, well below the shape, so the
# tilt transforms take the lower tail from a series there, see
# primarycensored_log_gamma_p(). The shape gradient of gamma_lccdf is
# inaccurate in the bulk (1e-3 at shape 2.5 and 7, 5e-3 at shape 20 and 24),
# so the tilt transforms take the upper tail from gamma_lcdf and use
# gamma_lccdf only beyond the point where the upper tail is below 1e-8. The
# points are also not at a switch between forms, as finite differences would
# step across it. The last case is a shape of 100.
exptilt_gradient_cases <- c(
  exptilt_stan_cases,
  list(list(
    dist_id = 2L, params = c(100, 10), pdist = pgamma,
    args = list(shape = 100, rate = 10)
  ))
)

# Compares the gradient with the finite difference gradient one component at a
# time, relative to the size of the component with a floor for tiny ones.
expect_gradient_close <- function(res, case, label, scale = 1) {
  tolerance <- rep(1e-4, 3) * scale
  allowed <- tolerance * pmax(abs(res$finite_diff), 1e-2)
  testthat::expect_true(
    all(abs(res$gradient - res$finite_diff) <= allowed),
    info = paste0(
      label, ": gradient ", toString(signif(res$gradient, 5)),
      ", finite difference ", toString(signif(res$finite_diff, 5))
    )
  )
}

test_that("tilted log CDFs have finite gradients matching finite
  differences", {
  model <- exptilt_gradient_model()
  # Direct, small window and small delay forms; both tails; d below and
  # above pwindow
  points <- list(
    list(d = 0.3, pwindow = 2, rho = 0.4),
    list(d = 1, pwindow = 1, rho = -0.2),
    list(d = 2.5, pwindow = 2, rho = 0.4),
    list(d = 6, pwindow = 3, rho = -0.15),
    list(d = 20, pwindow = 3, rho = 0.5),
    list(d = 7, pwindow = 1, rho = 0.3),
    list(d = 7, pwindow = 1, rho = -0.1),
    list(d = 2.5, pwindow = 2, rho = 1e-6),
    list(d = 2.5, pwindow = 2, rho = -1e-6),
    list(d = 12, pwindow = 7, rho = 1e-5),
    list(d = 0.0001, pwindow = 2, rho = 0.4),
    list(d = 0.0001, pwindow = 2, rho = -0.4),
    list(d = 4, pwindow = 2, rho = 0),
    # Lower tail of the larger shapes, where the shape gradient of
    # gamma_lcdf is inaccurate
    list(d = 1, pwindow = 1, rho = 0.3),
    list(d = 3, pwindow = 3, rho = -0.15),
    list(d = 2.5, pwindow = 1, rho = 0.3),
    list(d = 4, pwindow = 2, rho = 0.3)
  )
  for (case in exptilt_gradient_cases) {
    for (point in points) {
      if (case$dist_id != 18L &&
        exptilt_case_rate(case) + point$rho <= 0) {
        next
      }
      # The CDF underflows to zero, so there is no log CDF to differentiate
      if (case$dist_id == 2L && case$params[1] >= 100 && point$d < 0.5) {
        next
      }
      label <- exptilt_case_label(
        case,
        d = point$d, pwindow = point$pwindow, r = point$rho
      )
      res <- exptilt_gradient_at(
        model, case, point$d, point$pwindow, point$rho
      )
      expect_false(res$gradient_not_finite, info = label)
      expect_false(res$rejected, info = label)
      expect_length(res$gradient, 3)
      expect_true(all(is.finite(res$gradient)), info = label)
      expect_gradient_close(
        res, case, label,
        scale = if (is.null(point$scale)) 1 else point$scale
      )
    }
  }
})

test_that("the vectorised tilted log PMF has finite gradients matching
  finite differences", {
  model <- exptilt_gradient_model()
  # The small tilt form has an absolute error of about 1e-14 in the upper
  # tail, which finite differences of tiny PMF values amplify, so those
  # points stop at a delay where the PMF is not tiny. The small delay form
  # truncates at (r d)^2, so its gradient in r has a relative error of about
  # 1e-4, which the point that uses it allows for.
  points <- list(
    list(d = 12, pwindow = 3, rho = 0.4),
    list(d = 12, pwindow = 3, rho = -0.15),
    list(d = 6, pwindow = 2, rho = 1e-6),
    list(d = 6, pwindow = 10, rho = 1.5e-5, scale = 5),
    # Lower tail of the larger shapes
    list(d = 1, pwindow = 1, rho = 0.3),
    list(d = 3, pwindow = 3, rho = -0.15)
  )
  for (case in exptilt_gradient_cases[c(2, 4, 5, 8)]) {
    for (point in points) {
      if (case$dist_id != 18L &&
        exptilt_case_rate(case) + point$rho <= 0) {
        next
      }
      label <- exptilt_case_label(
        case,
        d = point$d, pwindow = point$pwindow, r = point$rho
      )
      res <- exptilt_gradient_at(
        model, case, point$d, point$pwindow, point$rho,
        vectorised = TRUE
      )
      expect_false(res$gradient_not_finite, info = label)
      expect_false(res$rejected, info = label)
      expect_true(all(is.finite(res$gradient)), info = label)
      expect_gradient_close(
        res, case, label,
        scale = if (is.null(point$scale)) 1 else point$scale
      )
    }
  }
})

test_that("the gamma tilt transform gradient is accurate in the bulk", {
  # The shape gradient of gamma_lccdf is off by about 1e-2 here, see the
  # note above exptilt_gradient_cases. The likelihood of several delays
  # is the case that showed it.
  model <- exptilt_gradient_model()
  case <- exptilt_stan_cases[[4]]
  for (pwindow in c(1, 3)) {
    for (rho in c(0.3, -0.1)) {
      for (d in c(5, 8, 12)) {
        res <- exptilt_gradient_at(
          model, case, d, pwindow, rho,
          vectorised = TRUE
        )
        label <- exptilt_case_label(case, d = d, pwindow = pwindow, r = rho)
        expect_true(all(is.finite(res$gradient)), info = label)
        expect_gradient_close(res, case, label)
      }
    }
  }
})

test_that("the gamma tilt transform is finite far in the upper tail", {
  # Beyond an upper tail of 1e-8 the lower tail is 1 to rounding and the
  # transform uses gamma_lccdf, which stays finite where log1m_exp of the
  # log CDF would not
  model <- exptilt_gradient_model()
  for (case in exptilt_stan_cases[c(3, 4)]) {
    for (d in c(40, 80, 200)) {
      res <- exptilt_gradient_at(model, case, d, 2, 0.3)
      label <- exptilt_case_label(case, d = d)
      expect_false(res$gradient_not_finite, info = label)
      expect_false(res$rejected, info = label)
      expect_true(all(is.finite(res$gradient)), info = label)
    }
  }
  upper <- vapply(
    c(30, 60, 120), log_tilt_transform_upper, numeric(1), 2L, -0.3,
    c(2.5, 0.4)
  )
  obj <- new_pcens(
    pgamma, dexpgrowth, list(r = 0.3),
    shape = 2.5, rate = 0.4
  )
  expect_equal(
    upper, .pcens_tilt_transform(obj, c(30, 60, 120), -0.3, upper = TRUE),
    tolerance = 1e-8
  )
})

test_that("a zero width primary window gives the delay CDF", {
  # The exact limit of the tilted window as the width goes to 0, which R
  # returns. The direct and small tilt forms divide by the width.
  d <- c(0.2, 1, 2.5, 6, 15)
  for (case in exptilt_stan_cases) {
    lower <- exptilt_case_lower(case)
    expected <- exptilt_case_cdf(case)(d)
    for (rho in c(-0.2, 0, 1e-9, 0.3)) {
      if (case$dist_id != 18L &&
        exptilt_case_rate(case) + rho <= 0) {
        next
      }
      info <- exptilt_case_label(case, r = rho)
      direct <- vapply(
        d, primarycensored_exptilt_lcdf, numeric(1),
        case$dist_id, case$params, 0, rho
      )
      expect_equal(exp(direct), expected, tolerance = 1e-12, info = info)
      lcdf <- vapply(
        d, primarycensored_lcdf, numeric(1),
        case$dist_id, case$params, 0, lower, Inf, 2L, rho
      )
      expect_equal(exp(lcdf), expected, tolerance = 1e-12, info = info)
      cdf <- vapply(
        d, primarycensored_cdf, numeric(1),
        case$dist_id, case$params, 0, lower, Inf, 2L, rho
      )
      expect_equal(cdf, expected, tolerance = 1e-12, info = info)
    }
  }
})

test_that("a zero width primary window gives the delay PMF", {
  for (case in exptilt_stan_cases[c(2, 4, 5)]) {
    cdf <- exptilt_case_cdf(case)
    lower <- if (case$dist_id == 18L) -Inf else 0
    pmf <- exp(primarycensored_sone_lpmf_vectorized(
      10, lower, Inf, case$dist_id, case$params, 0, 2L, 0.3
    ))
    expect_equal(
      pmf, diff(cdf(0:11)),
      tolerance = 1e-10,
      info = exptilt_case_label(case)
    )
  }
})

test_that("the gamma tilt transform gradients are accurate in the lower
  tail", {
  # The summed shape gradient of the vectorised PMF for delays 0 to 12,
  # where the lower tail is far below the shape for the early delays. The
  # PMF sums to a log probability of about -70 to -400, so finite
  # differences of it in Stan are noisy. The reference is central
  # differences of the log PMF from the reference integral.
  model <- exptilt_gradient_model()
  log_pmf_sum <- function(params, pwindow, rho) {
    cdf <- function(x) stats::pgamma(x, params[1], params[2])
    sum(log(diff(c(0, exptilt_reference(1:13, pwindow, rho, cdf)))))
  }
  for (case in exptilt_gradient_cases[c(5, 8)]) {
    for (pwindow in c(1, 3)) {
      rho <- 0.3
      theta <- c(case$params, rho)
      expected <- vapply(seq_along(theta), function(i) {
        h <- 1e-5 * theta[i]
        up <- down <- theta
        up[i] <- theta[i] + h
        down[i] <- theta[i] - h
        (log_pmf_sum(up[1:2], pwindow, up[3]) -
          log_pmf_sum(down[1:2], pwindow, down[3])) / (2 * h)
      }, numeric(1))
      # The rate is on the log scale in the model, with a Jacobian term
      expected[2] <- expected[2] * theta[2] + 1
      res <- exptilt_gradient_at(
        model, case, 12, pwindow, rho,
        vectorised = TRUE
      )
      label <- exptilt_case_label(case, pwindow = pwindow)
      expect_true(all(is.finite(res$gradient)), info = label)
      expect_true(
        all(abs(res$gradient - expected) <= 1e-5 * pmax(abs(expected), 1e-2)),
        info = paste0(
          label, ": gradient ", toString(signif(res$gradient, 6)),
          ", reference ", toString(signif(expected, 6))
        )
      )
    }
  }
})

# Gradient of a one argument function of the gamma lower tail in the shape,
# from a compiled model.
exptilt_log_gamma_p_model <- function() {
  testthat::skip_if_not_installed("cmdstanr")
  testthat::skip_if(
    is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))
  )
  functions <- pcd_load_stan_functions(
    wrap_in_block = TRUE, write_to_file = FALSE
  )
  code <- paste0(
    functions, "\n",
    "data {\n  real x;\n}\n",
    "parameters {\n  real a;\n}\n",
    "model {\n  target += primarycensored_log_gamma_p(x, a);\n}\n"
  )
  path <- file.path(tempdir(), "pcd_log_gamma_p.stan")
  writeLines(code, path)
  suppressMessages(suppressWarnings(cmdstanr::cmdstan_model(path)))
}

test_that("primarycensored_log_gamma_p is accurate in value and shape
  gradient well below the shape", {
  shapes <- c(0.3, 2.5, 20, 100, 1000)
  fractions <- c(0.001, 0.02, 0.1, 0.2, 0.3, 0.45, 0.55, 0.8, 1, 1.5)
  for (shape in shapes) {
    x <- shape * fractions
    expected <- stats::pgamma(x, shape, log.p = TRUE)
    actual <- vapply(x, primarycensored_log_gamma_p, numeric(1), shape)
    # Skip values the reference cannot represent
    keep <- is.finite(expected) & expected > -700
    expect_equal(
      actual[keep], expected[keep],
      tolerance = 1e-12, info = paste("shape", shape)
    )
  }
  model <- exptilt_log_gamma_p_model()
  for (shape in c(2.5, 20, 100)) {
    for (x in shape * c(0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 1)) {
      res <- stan_gradient_at( # nolint: object_usage_linter.
        model,
        data = list(x = x), init = list(a = shape)
      )
      h <- 1e-5 * shape
      expected <- (
        stats::pgamma(x, shape + h, log.p = TRUE) -
          stats::pgamma(x, shape - h, log.p = TRUE)
      ) / (2 * h)
      info <- paste("shape", shape, "x", x)
      expect_false(res$gradient_not_finite, info = info)
      # CmdStan prints the gradient to 6 significant digits
      expect_equal(res$gradient, expected, tolerance = 1e-5, info = info)
    }
  }
})

test_that("the normal tilted CDF is accurate in the moderate lower tail
  for a small tilt", {
  # Phi() loses relative precision for negative arguments down to -5, and
  # the direct form amplifies it by 1 / (|rho| w) for a small tilt
  case <- exptilt_stan_cases[[6]]
  cdf <- exptilt_case_cdf(case)
  for (pwindow in c(0.1, 1, 2.83)) {
    for (scaled in c(1.1e-4, 1e-3, 5e-3)) {
      for (sign in c(-1, 1)) {
        rho <- sign * scaled / pwindow
        d <- case$params[1] + case$params[2] * c(-4.95, -4.9, -4.6, -4.2)
        expected <- exptilt_reference(d, pwindow, rho, cdf)
        actual <- exp(vapply(
          d, primarycensored_exptilt_lcdf, numeric(1),
          case$dist_id, case$params, pwindow, rho
        ))
        expect_lt(
          max_rel_diff(actual, expected), 1e-9,
          label = exptilt_case_label(case, pwindow = pwindow, r = rho)
        )
      }
    }
  }
})
