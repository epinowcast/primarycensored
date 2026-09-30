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
  expect_identical(check_for_tilt_transform(2L, -0.4, c(2, 0.4)), 0L)
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
      for (rho in c(-0.2, 1e-8, 0.4)) {
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
        # The ODE path has absolute and relative tolerances of 1e-6
        ode <- vapply(
          d, primarycensored_numeric_cdf, numeric(1),
          case$dist_id, case$params, pwindow, 2L, rho
        )
        expect_equal(plain, ode, tolerance = 1e-5, info = info)
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
    list(d = 2.5, pwindow = 2, rho = 1e-6),
    list(d = 2.5, pwindow = 2, rho = -1e-6),
    list(d = 12, pwindow = 7, rho = 1e-5),
    list(d = 0.0001, pwindow = 2, rho = 0.4),
    list(d = 0.0001, pwindow = 2, rho = -0.4),
    list(d = 4, pwindow = 2, rho = 0)
  )
  for (case in exptilt_stan_cases) {
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
        model, case, point$d, point$pwindow, point$rho
      )
      expect_false(res$gradient_not_finite, info = label)
      expect_false(res$rejected, info = label)
      expect_length(res$gradient, 3)
      expect_true(all(is.finite(res$gradient)), info = label)
      expect_equal(
        res$gradient, res$finite_diff,
        tolerance = 1e-4, info = label
      )
    }
  }
})

test_that("the vectorised tilted log PMF has finite gradients matching
  finite differences", {
  model <- exptilt_gradient_model()
  points <- list(
    list(d = 12, pwindow = 3, rho = 0.4),
    list(d = 12, pwindow = 3, rho = -0.15),
    list(d = 10, pwindow = 2, rho = 1e-6),
    list(d = 10, pwindow = 10, rho = 2e-5)
  )
  for (case in exptilt_stan_cases[c(2, 4, 6)]) {
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
      expect_equal(
        res$gradient, res$finite_diff,
        tolerance = 1e-4, info = label
      )
    }
  }
})
