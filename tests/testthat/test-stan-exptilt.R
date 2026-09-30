skip_on_cran()

# Stan solutions for exponential, gamma and normal delays with an
# exponentially tilted primary (primary_id 2 with primary_params = r)

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

# Lower and upper transform from the pair
tilt_lower <- function(t, dist_id, xi, params) {
  log_tilt_transform_pair(t, dist_id, xi, params)[1]
}

tilt_upper <- function(t, dist_id, xi, params) {
  log_tilt_transform_pair(t, dist_id, xi, params)[2]
}

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
  for (dist_id in c(1L, 3L, 5L, 9L)) {
    expect_identical(check_for_analytical(dist_id, 2L), 0L)
    expect_identical(check_for_exptilt(dist_id, 2L), 0L)
  }
  expect_identical(check_for_analytical(4L, 1L), 0L)
  expect_identical(check_for_analytical(18L, 1L), 0L)
  expect_identical(check_for_analytical(2L, 1L), 1L)
  expect_identical(check_for_uniform_terms(2L, 2L), 0L)
})

test_that("check_for_analytical_params adds the admissibility of the tilt", {
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
        ts, tilt_lower, numeric(1), case$dist_id, xi, case$params
      )
      upper <- vapply(
        ts, tilt_upper, numeric(1), case$dist_id, xi,
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
  # Stan returns -Inf where a term is below the smallest double
  ts <- 10^seq(-8, 4, by = 0.5)
  for (shape in c(0.05, 0.3, 1, 7, 100, 1000)) {
    obj <- new_pcens(
      pgamma, dexpgrowth, list(r = 0.1),
      shape = shape, rate = 1.7
    )
    params <- c(shape, 1.7)
    for (xi in c(0, -0.5, 0.9)) {
      for (upper in c(FALSE, TRUE)) {
        stan_fun <- if (upper) tilt_upper else tilt_lower
        actual <- vapply(
          ts, stan_fun, numeric(1), 2L, xi, params
        )
        expected <- .pcens_tilt_transform(obj, ts, xi, upper = upper)
        keep <- is.finite(actual)
        expect_equal(
          actual[keep], expected[keep],
          tolerance = 1e-9,
          info = sprintf("shape %g, xi %g, upper %s", shape, xi, upper)
        )
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
        tilt_lower(t, case$dist_id, -0.1, case$params), -Inf
      )
    }
    total <- exp(tilt_upper(0, case$dist_id, -0.1, case$params))
    expect_identical(
      exp(tilt_upper(-2, case$dist_id, -0.1, case$params)),
      total
    )
  }
  expect_error(
    tilt_lower(1, 3L, 0, c(1, 1)),
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

test_that("primarycensored_exptilt_lcdf matches Monte Carlo samples", {
  set.seed(202)
  n <- 5000
  pwindow <- 2
  for (case in exptilt_stan_cases[c(2, 4, 5, 6, 7)]) {
    delay <- switch(as.character(case$dist_id),
      "4" = stats::rexp(n, case$params[1]),
      "2" = stats::rgamma(n, case$params[1], case$params[2]),
      "18" = stats::rnorm(n, case$params[1], case$params[2])
    )
    for (rho in c(-0.2, 0.5)) {
      primary <- vapply(
        seq_len(n), function(i) expgrowth_rng(0, pwindow, rho), numeric(1)
      )
      ks <- stats::ks.test(delay + primary, function(x) {
        exp(vapply(
          x, primarycensored_exptilt_lcdf, numeric(1),
          case$dist_id, case$params, pwindow, rho
        ))
      })
      expect_gt(
        ks$p.value, 1e-3,
        label = exptilt_case_label(case, pwindow = pwindow, r = rho)
      )
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
        # ODE tolerances are 1e-6, and about 1e-4 for a shape below 1
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
      # Covers the small window, small delay and direct forms
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
  # The intervals start at 0 and miss the mass of negative delays
  pmf <- primarycensored_sone_pmf_vectorized(
    60, -Inf, Inf, 18L, c(3, 2), 2, 2L, 0.2
  )
  cdf <- function(x) pnorm(x, 3, 2)
  expected <- diff(exptilt_reference(c(0, 61), 2, 0.2, cdf))
  expect_equal(sum(pmf), expected, tolerance = 1e-9)
})

# A minimal model whose target is the log CDF or the vectorised log PMF, for
# `stan_gradient_at()` from helper-stan-gradient.R
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

# Gradient cases, with a shape of 100 added
exptilt_gradient_cases <- c(
  exptilt_stan_cases,
  list(list(
    dist_id = 2L, params = c(100, 10), pdist = pgamma,
    args = list(shape = 100, rate = 10)
  ))
)

# Relative to each component, with a floor for tiny ones
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
    # Lower tail of the larger shapes
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
      # The CDF underflows
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
  # Finite differences amplify the small tilt error for tiny PMF values, and
  # the small delay form has a gradient in r with a relative error of about
  # 1e-4, which the scaled point allows for
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

test_that("the gamma tilt transform is finite far in the upper tail", {
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
    c(30, 60, 120), tilt_upper, numeric(1), 2L, -0.3,
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
  # Finite differences in Stan are noisy, so the reference is central
  # differences of the log PMF from the reference integral
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
      # Log scale rate with a Jacobian term
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

# Gradient of the gamma lower tail in the shape
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

# Shapes for which Stan's gamma_lcdf and gamma_lccdf shape gradients are NaN
exptilt_large_shape_cases <- list(
  list(
    dist_id = 2L, params = c(200, 20), pdist = pgamma,
    args = list(shape = 200, rate = 20)
  ),
  list(
    dist_id = 2L, params = c(250, 25), pdist = pgamma,
    args = list(shape = 250, rate = 25)
  ),
  list(
    dist_id = 2L, params = c(1000, 200), pdist = pgamma,
    args = list(shape = 1000, rate = 200)
  )
)

test_that("the gamma tilt transform gradients are finite and accurate for
  large shapes", {
  model <- exptilt_gradient_model()
  # Delays across the bulk and far into the upper tail
  points <- list(
    list(d = 14, pwindow = 1, rho = 0.2),
    list(d = 20, pwindow = 1, rho = 0.2),
    list(d = 20, pwindow = 3, rho = -0.1),
    list(d = 8, pwindow = 2, rho = 0.3),
    list(d = 11, pwindow = 1, rho = -0.1),
    list(d = 30, pwindow = 2, rho = 0.2)
  )
  for (case in exptilt_large_shape_cases) {
    for (point in points) {
      label <- exptilt_case_label(
        case,
        d = point$d, pwindow = point$pwindow, r = point$rho
      )
      res <- exptilt_gradient_at(
        model, case, point$d, point$pwindow, point$rho
      )
      expect_false(res$gradient_not_finite, info = label)
      expect_false(res$rejected, info = label)
      expect_true(all(is.finite(res$gradient)), info = label)
      expect_gradient_close(res, case, label)
    }
  }
})

test_that("the vectorised tilted log PMF has finite gradients for large
  shapes", {
  model <- exptilt_gradient_model()
  # A shape of 1000 has no delay range with every PMF above the smallest
  # double
  for (case in exptilt_large_shape_cases[1:2]) {
    for (point in list(
      list(d = 20, pwindow = 1, rho = 0.2),
      list(d = 12, pwindow = 1, rho = 0.2),
      list(d = 14, pwindow = 2, rho = -0.1)
    )) {
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
      expect_gradient_close(res, case, label)
    }
  }
})

# A model whose target is one statement in the data `x` and parameter `a`
exptilt_unary_model <- function(statement, name) {
  testthat::skip_if_not_installed("cmdstanr")
  testthat::skip_if(
    is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))
  )
  functions <- pcd_load_stan_functions(
    wrap_in_block = TRUE, write_to_file = FALSE
  )
  code <- paste0(
    functions, "\n",
    "data {\n  real x;\n  int which;\n}\n",
    "parameters {\n  real a;\n}\n",
    "model {\n  target += ", statement, ";\n}\n"
  )
  path <- file.path(tempdir(), paste0(name, ".stan"))
  writeLines(code, path)
  suppressMessages(suppressWarnings(cmdstanr::cmdstan_model(path)))
}

test_that("primarycensored_log_gamma_pq gives both tails in value and shape
  gradient for any shape", {
  fractions <- c(0.001, 0.3, 0.7, 0.9, 1, 1.05, 1.2, 1.6, 3, 8)
  shapes <- c(0.3, 2.5, 20, 200, 1000, 1e4)
  for (shape in shapes) {
    x <- shape * fractions
    lower <- stats::pgamma(x, shape, log.p = TRUE)
    upper <- stats::pgamma(x, shape, lower.tail = FALSE, log.p = TRUE)
    actual <- vapply(x, primarycensored_log_gamma_pq, numeric(2), shape)
    for (k in 1:2) {
      expected <- if (k == 1) lower else upper
      keep <- is.finite(expected) & expected > -1e5
      expect_equal(
        actual[k, keep], expected[keep],
        tolerance = 1e-10, info = paste("tail", k, "shape", shape)
      )
    }
  }
  model <- exptilt_unary_model(
    "primarycensored_log_gamma_pq(x, a)[which]", "pcd_log_gamma_pq"
  )
  for (shape in c(2.5, 20, 200, 1000)) {
    for (x in shape * c(0.05, 0.3, 0.7, 1, 1.3, 2, 4)) {
      for (k in 1:2) {
        res <- stan_gradient_at( # nolint: object_usage_linter.
          model,
          data = list(x = x, which = k), init = list(a = shape)
        )
        h <- 1e-4 * shape
        reference <- function(s) {
          stats::pgamma(x, s, lower.tail = k == 1, log.p = TRUE)
        }
        expected <- (-reference(shape + 2 * h) + 8 * reference(shape + h) -
          8 * reference(shape - h) + reference(shape - 2 * h)) / (12 * h)
        info <- paste("tail", k, "shape", shape, "x", x)
        expect_false(res$gradient_not_finite, info = info)
        expect_false(res$rejected, info = info)
        # CmdStan prints 6 significant digits
        expect_lte(
          abs(res$gradient - expected),
          1e-5 * max(abs(expected), 1e-6), label = info
        )
      }
    }
  }
})

test_that("the normal log CDF has an exact gradient in the deep lower tail", {
  model <- exptilt_unary_model(
    "primarycensored_log_std_normal_cdf(a)", "pcd_log_std_normal_cdf"
  )
  z <- c(-5, -20, -36.9, -37.1, -38, -45, -60, -150)
  actual <- vapply(z, primarycensored_log_std_normal_cdf, numeric(1))
  expect_equal(actual, stats::pnorm(z, log.p = TRUE), tolerance = 1e-13)
  for (zz in z) {
    res <- stan_gradient_at( # nolint: object_usage_linter.
      model,
      data = list(x = 0, which = 1L), init = list(a = zz)
    )
    # d/dz log Phi(z) = phi(z) / Phi(z)
    expected <- exp(stats::dnorm(zz, log = TRUE) -
      stats::pnorm(zz, log.p = TRUE))
    expect_false(res$gradient_not_finite, info = as.character(zz))
    # CmdStan prints the gradient to 6 significant digits
    expect_equal(
      res$gradient, expected,
      tolerance = 1e-5, info = as.character(zz)
    )
  }
})

# Reference log CDF for a normal delay, scaled to keep a CDF far below the
# smallest double, and its five point central difference in the parameters
exptilt_normal_log_reference <- function(d, pwindow, rho, mu, sigma) {
  log_integrand <- function(z) {
    window_density <- exptilt_window_density( # nolint: object_usage_linter.
      z, pwindow, rho
    )
    stats::pnorm((d - z - mu) / sigma, log.p = TRUE) + log(window_density)
  }
  shift <- max(log_integrand(0), log_integrand(pwindow))
  integral <- stats::integrate(
    function(z) exp(vapply(z, log_integrand, numeric(1)) - shift),
    0, pwindow,
    rel.tol = 1e-13, abs.tol = 0, subdivisions = 2000L
  )$value
  shift + log(integral)
}

exptilt_normal_log_gradient <- function(d, pwindow, rho, mu, sigma) {
  theta <- c(mu, sigma, rho)
  f <- function(theta) {
    exptilt_normal_log_reference(d, pwindow, theta[3], theta[1], theta[2])
  }
  grad <- vapply(1:3, function(i) {
    h <- 1e-3 * abs(theta[i])
    at <- function(step) {
      theta[i] <- theta[i] + step * h
      f(theta)
    }
    (-at(2) + 8 * at(1) - 8 * at(-1) + at(-2)) / (12 * h)
  }, numeric(1))
  # Log scale standard deviation with a Jacobian term
  grad[2] <- grad[2] * sigma + 1
  grad
}

test_that("the normal tilted log CDF gradients are accurate in the deep
  lower tail", {
  # Arguments below -37, where std_normal_lcdf() is inaccurate
  model <- exptilt_gradient_model()
  case <- exptilt_stan_cases[[6]]
  case$params <- c(8, 3)
  for (pwindow in c(1, 3, 1e-3)) {
    rho <- -3
    d <- -96.5
    expected <- exptilt_normal_log_gradient(d, pwindow, rho, 8, 3)
    res <- exptilt_gradient_at(model, case, d, pwindow, rho)
    label <- paste("pwindow", pwindow)
    expect_false(res$gradient_not_finite, info = label)
    expect_true(all(is.finite(res$gradient)), info = label)
    expect_true(
      all(abs(res$gradient - expected) <= 1e-5 * pmax(abs(expected), 1e-2)),
      info = paste0(
        label, ": gradient ", toString(signif(res$gradient, 6)),
        ", reference ", toString(signif(expected, 6))
      )
    )
  }
})

# Reference log of the upper tail 1 - F_rho(d) of a gamma delay
exptilt_gamma_log_upper <- function(d, pwindow, rho, shape, rate) {
  log_integrand <- function(z) {
    window_density <- exptilt_window_density( # nolint: object_usage_linter.
      z, pwindow, rho
    )
    stats::pgamma(d - z, shape, rate, lower.tail = FALSE, log.p = TRUE) +
      log(window_density)
  }
  shift <- log_integrand(0)
  integral <- stats::integrate(
    function(z) exp(vapply(z, log_integrand, numeric(1)) - shift),
    0, pwindow,
    rel.tol = 1e-13, abs.tol = 0
  )$value
  shift + log(integral)
}

test_that("the rate gradient of the tilted log CDF is accurate where the CDF
  is close to 1", {
  # The log CDF is about -U for an upper tail U far below 1
  model <- exptilt_unary_model(
    paste0(
      "primarycensored_lcdf(x | 2, {200.0, a}, 1.0, 0.0, ",
      "positive_infinity(), 2, {0.2})"
    ),
    "pcd_exptilt_lcdf_rate"
  )
  for (d in c(12, 14, 16, 18, 20)) {
    res <- stan_gradient_at( # nolint: object_usage_linter.
      model,
      data = list(x = d, which = 1L), init = list(a = 20)
    )
    log_upper <- function(rate) {
      exptilt_gamma_log_upper(d, 1, 0.2, 200, rate)
    }
    h <- 1e-4
    slope <- (
      -log_upper(20 * exp(2 * h)) + 8 * log_upper(20 * exp(h)) -
        8 * log_upper(20 * exp(-h)) + log_upper(20 * exp(-2 * h))
    ) / (12 * h * 20)
    upper <- exp(log_upper(20))
    expected <- -upper * slope / (1 - upper)
    expect_equal(
      res$gradient, expected,
      tolerance = 1e-4, info = as.character(d)
    )
  }
})

test_that("the small tilt forms have a tilt gradient within the documented
  bound", {
  # The relative error is about |rho| w / 6, up to 2e-5 at the threshold
  model <- exptilt_gradient_model()
  points <- list(
    list(params = c(100, 10), d = 3.744, pwindow = 3, rho = 1e-5),
    list(params = c(100, 10), d = 12, pwindow = 3, rho = 1e-5),
    list(params = c(2.5, 0.4), d = 6, pwindow = 10, rho = 9.9e-6),
    list(params = c(2.5, 0.4), d = 6, pwindow = 2, rho = -4.9e-5),
    list(params = c(20, 4), d = 5, pwindow = 1, rho = 9e-5),
    list(params = c(2.5, 0.4), d = 0.0009, pwindow = 2, rho = 0.04)
  )
  for (point in points) {
    case <- list(dist_id = 2L, params = point$params)
    cdf <- function(x) stats::pgamma(x, point$params[1], point$params[2])
    log_cdf <- function(rho) {
      log(exptilt_reference(point$d, point$pwindow, rho, cdf))
    }
    h <- max(abs(point$rho) * 0.05, 1e-6)
    expected <- (
      -log_cdf(point$rho + 2 * h) + 8 * log_cdf(point$rho + h) -
        8 * log_cdf(point$rho - h) + log_cdf(point$rho - 2 * h)
    ) / (12 * h)
    res <- exptilt_gradient_at(
      model, case, point$d, point$pwindow, point$rho
    )
    label <- exptilt_case_label(
      case,
      d = point$d, pwindow = point$pwindow, r = point$rho
    )
    expect_lt(abs(res$gradient[3] / expected - 1), 3e-5, label = label)
  }
})

test_that("the small tilt and direct forms have accurate gradients for large
  shapes", {
  # CmdStan finite differences are inaccurate at a shape of 1000, so the
  # reference is a central difference of the reference integral
  model <- exptilt_gradient_model()
  points <- list(
    list(case = 3, d = 5, pwindow = 2, rho = -1e-5),
    list(case = 3, d = 5.3, pwindow = 3, rho = -1e-3),
    list(case = 3, d = 4.6, pwindow = 1, rho = 1e-5),
    list(case = 1, d = 10, pwindow = 2, rho = -1e-5),
    list(case = 1, d = 11, pwindow = 1, rho = 1e-5),
    list(case = 1, d = 9, pwindow = 3, rho = -1e-3)
  )
  for (point in points) {
    case <- exptilt_large_shape_cases[[point$case]]
    theta <- c(case$params, point$rho)
    log_cdf <- function(theta) {
      cdf <- function(x) stats::pgamma(x, theta[1], theta[2])
      log(exptilt_reference(point$d, point$pwindow, theta[3], cdf))
    }
    steps <- c(
      1e-5 * theta[1:2], max(abs(theta[3]) * 0.05, 1e-7)
    )
    expected <- vapply(1:3, function(i) {
      at <- function(step) {
        theta[i] <- theta[i] + step * steps[i]
        log_cdf(theta)
      }
      (-at(2) + 8 * at(1) - 8 * at(-1) + at(-2)) / (12 * steps[i])
    }, numeric(1))
    # Log scale rate with a Jacobian term
    expected[2] <- expected[2] * theta[2] + 1
    res <- exptilt_gradient_at(
      model, case, point$d, point$pwindow, point$rho
    )
    label <- exptilt_case_label(
      case,
      d = point$d, pwindow = point$pwindow, r = point$rho
    )
    expect_true(all(is.finite(res$gradient)), info = label)
    expect_true(
      all(abs(res$gradient - expected) <= 3e-5 * pmax(abs(expected), 1e-2)),
      info = paste0(
        label, ": gradient ", toString(signif(res$gradient, 6)),
        ", reference ", toString(signif(expected, 6))
      )
    )
  }
})
