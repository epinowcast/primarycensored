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

# Shapes for which Stan's gamma_lcdf and gamma_lccdf shape gradients are NaN
exptilt_large_shape_cases <- list(
  list(dist_id = 2L, params = c(200, 20)),
  list(dist_id = 2L, params = c(250, 25)),
  list(dist_id = 2L, params = c(1000, 200))
)

# Gradient cases, with a shape of 100 added
exptilt_gradient_cases <- c(
  exptilt_stan_cases,
  list(list(dist_id = 2L, params = c(100, 10)))
)

# Lower and upper transform from the pair
tilt_lower <- function(t, dist_id, xi, params) {
  log_tilt_transform_pair( # nolint: object_usage_linter.
    t, dist_id, xi, params
  )[1]
}

tilt_upper <- function(t, dist_id, xi, params) {
  log_tilt_transform_pair( # nolint: object_usage_linter.
    t, dist_id, xi, params
  )[2]
}

exptilt_case_rate <- function(case) {
  switch(as.character(case$dist_id),
    "4" = case$params[1],
    "2" = case$params[2],
    Inf
  )
}

# Whether the tilt is admissible for the delay
exptilt_case_ok <- function(case, rho) {
  case$dist_id == 18L || exptilt_case_rate(case) + rho > 0
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

per_delay_exptilt_lcdf <- function(delays, dist_id, params, pwindow, rho) {
  lower <- if (dist_id == 18L) -Inf else 0
  vapply(
    delays, primarycensored_lcdf, numeric(1), # nolint: object_usage_linter.
    dist_id, params, pwindow, lower, Inf, 2L, rho
  )
}

# A model with the Stan functions followed by `code`
exptilt_stan_model <- function(name, code) {
  testthat::skip_if_not_installed("cmdstanr")
  testthat::skip_if(
    is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))
  )
  functions <- pcd_load_stan_functions(
    wrap_in_block = TRUE, write_to_file = FALSE
  )
  path <- file.path(tempdir(), paste0(name, ".stan"))
  writeLines(paste0(functions, "\n", code), path)
  suppressMessages(suppressWarnings(cmdstanr::cmdstan_model(path)))
}

# The ODE branch of primarycensored_cdf() for given delays, as the CDF of a
# fixed parameter run of a model that calls the ODE solver directly
exptilt_ode_model <- function() {
  exptilt_stan_model("pcd_exptilt_ode", paste(
    "data {",
    "  int n;",
    "  array[n] real d;",
    "  int dist_id;",
    "  int n_params;",
    "  array[n_params] real params;",
    "  real pwindow;",
    "  array[1] real primary_params;",
    "}",
    "generated quantities {",
    "  array[n] real cdf;",
    "  for (i in 1:n) {",
    "    real lower_bound = dist_has_positive_support(dist_id)",
    "      ? fmax(d[i] - pwindow, 0) : d[i] - pwindow;",
    "    array[n_params + 1] real theta =",
    "      append_array(params, primary_params);",
    "    array[4] int ids = {dist_id, 2, n_params, 1};",
    "    cdf[i] = ode_rk45(",
    "      primarycensored_ode, rep_vector(0.0, 1), lower_bound, {d[i]},",
    "      theta, {d[i], pwindow}, ids",
    "    )[1, 1];",
    "  }",
    "}",
    sep = "\n"
  ))
}

exptilt_ode_cdf <- function(model, case, d, pwindow, rho) {
  fit <- model$sample(
    data = list(
      n = length(d), d = as.array(d), dist_id = case$dist_id,
      n_params = length(case$params), params = as.array(case$params),
      pwindow = pwindow, primary_params = as.array(rho)
    ),
    fixed_param = TRUE, chains = 1, iter_sampling = 1, refresh = 0,
    show_messages = FALSE, sig_figs = 18
  )
  as.numeric(fit$draws("cdf", format = "matrix"))
}

# A minimal model whose target is the log CDF or the vectorised log PMF, for
# `stan_gradient_at()` from helper-stan-gradient.R
exptilt_gradient_model <- function() {
  exptilt_stan_model("pcd_exptilt_gradient", paste(
    "data {",
    "  int dist_id;",
    "  int n_params;",
    "  int vectorised;",
    "  real d;",
    "  real pwindow;",
    "  real L;",
    "}",
    "parameters {",
    "  real p1;",
    "  real<lower=0> p2;",
    "  real rho;",
    "}",
    "model {",
    "  array[2] real all_params = {p1, p2};",
    "  array[n_params] real params = all_params[1:n_params];",
    "  if (vectorised) {",
    "    target += sum(primarycensored_sone_lpmf_vectorized(",
    "      to_int(d), L, positive_infinity(), dist_id, params, pwindow, 2,",
    "      {rho}",
    "    ));",
    "  } else {",
    "    target += primarycensored_lcdf(",
    "      d | dist_id, params, pwindow, L, positive_infinity(), 2, {rho}",
    "    );",
    "  }",
    "}",
    sep = "\n"
  ))
}

# A model whose target is one statement in the data `x` and parameter `a`
exptilt_unary_model <- function(statement, name) {
  exptilt_stan_model(name, paste0(
    "data {\n  real x;\n  int which;\n}\n",
    "parameters {\n  real a;\n}\n",
    "model {\n  target += ", statement, ";\n}\n"
  ))
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

# Gradient of a unary model at `a`, with the two flags of stan_gradient_at()
exptilt_unary_gradient <- function(model, x, a, which = 1L) {
  stan_gradient_at( # nolint: object_usage_linter.
    model,
    data = list(x = x, which = which), init = list(a = a)
  )
}

test_that("the checks follow the delay, the primary and the tilt", {
  admissible <- function(dist_id, params, xi) {
    vapply(
      xi, function(x) check_for_tilt_transform(dist_id, x, params),
      integer(1)
    )
  }
  expect_identical(
    admissible(4L, 0.3, c(-0.5, 0.29, 0.3, 1)), c(1L, 1L, 0L, 0L)
  )
  expect_identical(
    admissible(2L, c(2, 0.4), c(-0.3, 0.39, 0.4, 0)), c(1L, 1L, 0L, 1L)
  )
  expect_identical(admissible(18L, c(3, 2), c(-50, 50)), c(1L, 1L))
  for (dist_id in c(1L, 3L, 5L, 12L, 26L, 27L)) {
    expect_identical(check_for_tilt_transform(dist_id, 0, c(1, 1)), 0L)
  }
  for (dist_id in c(2L, 4L, 18L)) {
    expect_identical(check_for_exptilt(dist_id, 2L), 1L)
    expect_identical(check_for_exptilt(dist_id, 1L), 0L)
    for (pwindow in c(1, 7)) {
      expect_identical(
        check_for_analytical_vectorized(dist_id, 2L, pwindow), 1L
      )
    }
    for (pwindow in c(0.5, 1.5)) {
      expect_identical(
        check_for_analytical_vectorized(dist_id, 2L, pwindow), 0L
      )
    }
  }
  for (dist_id in c(1L, 3L, 5L, 9L)) {
    expect_identical(check_for_exptilt(dist_id, 2L), 0L)
  }
  expect_identical(check_for_analytical_vectorized(3L, 2L, 1), 0L)
  expect_identical(check_for_analytical(4L, 1L), 0L)
  expect_identical(check_for_analytical(18L, 1L), 0L)
  expect_identical(check_for_analytical(2L, 1L), 1L)
  expect_identical(check_for_uniform_terms(2L, 2L), 0L)
  checks <- list(
    list(2L, c(2, 0.4), 2L, 0.5, 1L),
    list(2L, c(2, 0.4), 2L, -0.3, 1L),
    list(2L, c(2, 0.4), 2L, -0.4, 0L),
    list(4L, 0.3, 2L, -0.5, 0L),
    list(4L, 0.3, 2L, -0.25, 1L),
    list(18L, c(3, 2), 2L, -10, 1L),
    list(2L, c(2, 0.4), 1L, numeric(0), 1L),
    list(3L, c(2, 1), 2L, 0.3, 0L),
    list(26L, c(0, 1, 2, 0.5, 0.5), 2L, 0.3, 1L)
  )
  for (check in checks) {
    expect_identical(
      do.call(check_for_analytical_params, check[1:4]), check[[5]]
    )
  }
})

test_that("Stan tilt transforms and moments match the R versions", {
  ts <- c(1e-3, 0.4, 1, 3.5, 12, 40)
  for (case in exptilt_stan_cases) {
    obj <- exptilt_object(list(pdist = case$pdist, args = case$args), 0.1)
    for (xi in c(0, -1, -0.25, 0.1, 0.25)) {
      if (!check_for_tilt_transform(case$dist_id, xi, case$params)) {
        next
      }
      info <- exptilt_case_label(case, xi = xi)
      expect_equal(
        vapply(ts, tilt_lower, numeric(1), case$dist_id, xi, case$params),
        .pcens_tilt_transform(obj, ts, xi), info = info
      )
      expect_equal(
        vapply(ts, tilt_upper, numeric(1), case$dist_id, xi, case$params),
        .pcens_tilt_transform(obj, ts, xi, upper = TRUE), info = info
      )
    }
    expect_equal(
      t(vapply(
        ts, primarycensored_tilt_moments, numeric(3), case$dist_id,
        case$params
      )),
      unname(.pcens_tilt_moments(obj, ts)),
      tolerance = 1e-9, info = case$dist_id
    )
  }
  # Stan returns -Inf where a gamma term is below the smallest double
  ts <- 10^seq(-8, 4, by = 0.5)
  for (shape in c(0.05, 0.3, 1, 7, 100, 1000)) {
    obj <- new_pcens(
      pgamma, dexpgrowth, list(r = 0.1), shape = shape, rate = 1.7
    )
    for (xi in c(0, -0.5, 0.9)) {
      for (upper in c(FALSE, TRUE)) {
        actual <- vapply(
          ts, if (upper) tilt_upper else tilt_lower, numeric(1), 2L, xi,
          c(shape, 1.7)
        )
        expected <- .pcens_tilt_transform(obj, ts, xi, upper = upper)
        keep <- is.finite(actual)
        info <- sprintf("shape %g, xi %g, upper %s", shape, xi, upper)
        expect_equal(
          actual[keep], expected[keep], tolerance = 1e-9, info = info
        )
        expect_true(all(expected[!keep] < -600), info = info)
      }
    }
  }
  # The gamma upper tail far above the mean
  ts <- c(30, 60, 120)
  obj <- new_pcens(
    pgamma, dexpgrowth, list(r = 0.3), shape = 2.5, rate = 0.4
  )
  expect_equal(
    vapply(ts, tilt_upper, numeric(1), 2L, -0.3, c(2.5, 0.4)),
    .pcens_tilt_transform(obj, ts, -0.3, upper = TRUE),
    tolerance = 1e-8
  )
})

test_that("Stan tilt transforms are 0 or the total below the support", {
  for (case in exptilt_stan_cases[c(2, 4)]) {
    for (t in c(-3, -1e-9, 0)) {
      expect_identical(
        tilt_lower(t, case$dist_id, -0.1, case$params), -Inf
      )
    }
    expect_identical(
      exp(tilt_upper(-2, case$dist_id, -0.1, case$params)),
      exp(tilt_upper(0, case$dist_id, -0.1, case$params))
    )
  }
  expect_error(
    tilt_lower(1, 3L, 0, c(1, 1)), "Invalid distribution identifier"
  )
  expect_error(
    primarycensored_tilt_moments(1, 3L, c(1, 1)),
    "Invalid distribution identifier"
  )
})

test_that("primarycensored_exptilt_lcdf matches a reference integral and R", {
  rhos <- c(
    -1, -0.5, -0.05, -1e-4, -1e-5, -1e-8, 1e-8, 1e-5, 1e-4, 0.05, 0.5, 1
  )
  for (case in exptilt_stan_cases) {
    cdf <- exptilt_case_cdf(case)
    for (pwindow in c(0.5, 1, 2, 7)) {
      d <- sort(c(
        1e-6, 1e-3, 0.3 * pwindow, pwindow - 1e-3, pwindow, pwindow + 1e-3,
        2, 3, 6, 12, 25,
        if (case$dist_id == 18L) c(-10, -3, -0.5)
      ))
      for (rho in c(rhos, -0.2, 0.3, -1e-6, 1e-6)) {
        if (!exptilt_case_ok(case, rho)) {
          next
        }
        actual <- exp(vapply(
          d, primarycensored_exptilt_lcdf, numeric(1),
          case$dist_id, case$params, pwindow, rho
        ))
        label <- exptilt_case_label(case, pwindow = pwindow, r = rho)
        expect_lt(
          max_rel_diff(actual, exptilt_reference(d, pwindow, rho, cdf)),
          1e-7, label = label
        )
        expect_equal(
          actual,
          pcens_cdf(exptilt_object(
            list(pdist = case$pdist, args = case$args), rho
          ), d, pwindow),
          tolerance = 1e-9, info = label
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

test_that("the normal tilted CDF is accurate in the lower tail for small
  tilts", {
  grid <- exptilt_normal_tail_grid()
  expected <- exptilt_normal_tail_reference(grid)
  actual <- vapply(seq_len(nrow(grid)), function(i) {
    primarycensored_lcdf(
      grid$d[i], 18L, c(-4, 0.3), grid$pwindow[i], -Inf, Inf, 2L, grid$rho[i]
    )
  }, numeric(1))
  expect_lt(max(abs(expm1(actual - expected))), 3e-7)
  # Moderate lower tail, with a tighter tolerance
  case <- exptilt_stan_cases[[6]]
  cdf <- exptilt_case_cdf(case)
  for (pwindow in c(0.1, 1, 2.83)) {
    for (rho in c(-1, 1) * rep(c(1.1e-4, 1e-3, 5e-3), each = 2) / pwindow) {
      d <- case$params[1] + case$params[2] * c(-4.95, -4.9, -4.6, -4.2)
      actual <- exp(vapply(
        d, primarycensored_exptilt_lcdf, numeric(1),
        case$dist_id, case$params, pwindow, rho
      ))
      expect_lt(
        max_rel_diff(actual, exptilt_reference(d, pwindow, rho, cdf)), 1e-9,
        label = exptilt_case_label(case, pwindow = pwindow, r = rho)
      )
    }
  }
})

test_that("the tilted log CDF is accurate for gamma delays with large
  shapes", {
  q <- c(6, 8, 9.5, 10, 10.5)
  for (shape in c(500, 1000, 5000)) {
    for (pwindow in c(1, 7)) {
      for (rho in c(-1e-3, -1e-5, 2e-5, 1.5e-4, 1e-3, 1e-2)) {
        expected <- exptilt_gamma_log_reference(
          q, pwindow, rho, shape, shape / 10
        )
        actual <- vapply(
          q, primarycensored_lcdf, numeric(1),
          2L, c(shape, shape / 10), pwindow, 0, Inf, 2L, rho
        )
        expect_lt(
          max(abs(expm1(actual - expected))), 1e-6,
          label = sprintf(
            "shape = %g, pwindow = %g, r = %g", shape, pwindow, rho
          )
        )
      }
    }
  }
})

test_that("primarycensored_lcdf and primarycensored_cdf use the analytical
  solution and agree with the ODE path", {
  d <- c(0.2, 1, 2.5, 6, 15)
  ode_model <- exptilt_ode_model()
  for (case in exptilt_stan_cases) {
    lower <- exptilt_case_lower(case)
    cdf <- exptilt_case_cdf(case)
    for (pwindow in c(1, 3)) {
      for (rho in c(-1, -0.5, -1e-8, 1e-8, 0.5, 1)) {
        if (!exptilt_case_ok(case, rho)) {
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
        ode <- exptilt_ode_cdf(ode_model, case, d, pwindow, rho)
        expect_lt(max(abs(plain - ode)), 1e-4, label = info)
      }
    }
  }
})

test_that("inadmissible tilts use the ODE path", {
  ode_model <- exptilt_ode_model()
  cases <- list(
    list(dist_id = 4L, params = 0.3, rho = -0.5),
    list(dist_id = 4L, params = 0.3, rho = -0.3),
    list(dist_id = 2L, params = c(2.5, 0.4), rho = -0.5),
    list(dist_id = 2L, params = c(2.5, 0.4), rho = -0.4)
  )
  d <- c(0.5, 2, 5, 10)
  for (case in cases) {
    expect_identical(
      check_for_analytical_params(case$dist_id, case$params, 2L, case$rho),
      0L
    )
    ode <- exptilt_ode_cdf(ode_model, case, d, 2, case$rho)
    plain <- vapply(
      d, primarycensored_cdf, numeric(1),
      case$dist_id, case$params, 2, 0, Inf, 2L, case$rho
    )
    expect_equal(plain, ode, tolerance = 1e-12)
    lcdf <- vapply(
      d, primarycensored_lcdf, numeric(1),
      case$dist_id, case$params, 2, 0, Inf, 2L, case$rho
    )
    expect_equal(lcdf, log(ode), tolerance = 1e-12)
  }
  expect_error(
    primarycensored_analytical_lcdf(2, 2L, c(2.5, 0.4), 2, 0, Inf, 2L, -0.5),
    "tilted delay distribution"
  )
})

test_that("the numerical path returns 0 for non-positive delays without a
  lower truncation", {
  cases <- list(
    list(dist_id = 3L, params = c(1.5, 2), rho = 0.2),
    list(dist_id = 4L, params = 1, rho = -2)
  )
  for (case in cases) {
    for (d in c(0, -1)) {
      expect_identical(
        primarycensored_cdf(
          d, case$dist_id, case$params, 1, -Inf, Inf, 2L, case$rho
        ), 0
      )
    }
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

test_that("a zero width primary window gives the delay CDF and PMF", {
  d <- c(0.2, 1, 2.5, 6, 15)
  for (case in exptilt_stan_cases) {
    lower <- exptilt_case_lower(case)
    expected <- exptilt_case_cdf(case)(d)
    for (rho in c(-0.2, 0, 1e-9, 0.3)) {
      if (!exptilt_case_ok(case, rho)) {
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
  for (case in exptilt_stan_cases[c(2, 4, 5)]) {
    pmf <- exp(primarycensored_sone_lpmf_vectorized(
      10, 0, Inf, case$dist_id, case$params, 0, 2L, 0.3
    ))
    expect_equal(
      pmf, diff(exptilt_case_cdf(case)(0:11)),
      tolerance = 1e-10, info = exptilt_case_label(case)
    )
  }
})

test_that("the vectorised tilted CDF matches the per delay CDF", {
  n <- 31L
  for (case in exptilt_stan_cases) {
    for (pwindow in c(1, 2, 7)) {
      # Covers the small window, small delay and direct forms
      rhos <- c(-0.3, -2e-4, -2e-5, -1e-9, 0, 1e-9, 1e-5, 2e-5, 2e-4, 0.4)
      for (rho in rhos) {
        if (!exptilt_case_ok(case, rho)) {
          next
        }
        for (start in c(1L, 5L)) {
          vectorised <- primarycensored_analytical_lcdf_vectorized(
            start, n, case$dist_id, case$params, pwindow, 2L, rho
          )
          expect_length(vectorised, n)
          expect_identical(
            vectorised[start:n],
            per_delay_exptilt_lcdf(
              start:n, case$dist_id, case$params, pwindow, rho
            ),
            info = exptilt_case_label(
              case, pwindow = pwindow, r = rho, start = start
            )
          )
        }
      }
    }
  }
  # Delays with the small delay form and the direct form in one call
  for (case in exptilt_stan_cases[c(2, 4, 5)]) {
    expect_identical(
      primarycensored_analytical_lcdf_vectorized(
        1L, 25L, case$dist_id, case$params, 10, 2L, 2e-5
      )[1:25],
      per_delay_exptilt_lcdf(1:25, case$dist_id, case$params, 10, 2e-5)
    )
  }
})

test_that("primarycensored_lcdf_vectorized uses the tilted shared terms", {
  for (case in exptilt_stan_cases) {
    expect_identical(
      primarycensored_lcdf_vectorized(
        1L, 20L, case$dist_id, case$params, 3, 2L, 0.25
      ),
      primarycensored_analytical_lcdf_vectorized(
        1L, 20L, case$dist_id, case$params, 3, 2L, 0.25
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
    for (setting in settings) {
      for (pwindow in c(1, 3)) {
        for (rho in c(-0.2, 1e-9, 0.3)) {
          if (!exptilt_case_ok(case, rho)) {
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
  # The intervals start at 0 and miss the mass of negative delays
  pmf <- primarycensored_sone_pmf_vectorized(
    60, -Inf, Inf, 18L, c(3, 2), 2, 2L, 0.2
  )
  expected <- diff(exptilt_reference(
    c(0, 61), 2, 0.2, function(x) pnorm(x, 3, 2)
  ))
  expect_equal(sum(pmf), expected, tolerance = 1e-9)
})

# nolint start: object_usage_linter.
# Gradient of the gradient model at a point against the CmdStan finite
# differences, relative to each component with a floor for tiny ones
expect_fd_gradient <- function(model, case, point, vectorised, scale = 1) {
  label <- exptilt_case_label(
    case, d = point$d, pwindow = point$pwindow, r = point$rho
  )
  res <- exptilt_gradient_at(
    model, case, point$d, point$pwindow, point$rho,
    vectorised = vectorised
  )
  expect_false(res$gradient_not_finite, info = label)
  expect_false(res$rejected, info = label)
  expect_length(res$gradient, 3)
  expect_true(all(is.finite(res$gradient)), info = label)
  allowed <- 1e-4 * scale * pmax(abs(res$finite_diff), 1e-2)
  expect_true(
    all(abs(res$gradient - res$finite_diff) <= allowed),
    info = paste0(
      label, ": gradient ", toString(signif(res$gradient, 5)),
      ", finite difference ", toString(signif(res$finite_diff, 5))
    )
  )
}
# nolint end

# Finite differences amplify the small tilt error for tiny PMF values, and
# the small delay form has a gradient in r with a relative error of about
# 1e-4, which the scaled point allows for. A shape of 1000 has no delay
# range with every PMF above the smallest double.
exptilt_fd_groups <- function() {
  list(
    list(
      cases = exptilt_gradient_cases, vectorised = FALSE,
      points = data.frame(
        d = c(
          0.3, 1, 2.5, 6, 20, 7, 7, 2.5, 2.5, 12, 1e-4, 1e-4, 4, 1, 3, 2.5, 4
        ),
        pwindow = c(2, 1, 2, 3, 3, 1, 1, 2, 2, 7, 2, 2, 2, 1, 3, 1, 2),
        rho = c(
          0.4, -0.2, 0.4, -0.15, 0.5, 0.3, -0.1, 1e-6, -1e-6, 1e-5, 0.4, -0.4,
          0, 0.3, -0.15, 0.3, 0.3
        )
      )
    ),
    list(
      cases = exptilt_gradient_cases[c(2, 4, 5, 8)], vectorised = TRUE,
      points = data.frame(
        d = c(12, 12, 6, 6, 1, 3), pwindow = c(3, 3, 2, 10, 1, 3),
        rho = c(0.4, -0.15, 1e-6, 1.5e-5, 0.3, -0.15),
        scale = c(1, 1, 1, 5, 1, 1)
      )
    ),
    list(
      cases = exptilt_large_shape_cases, vectorised = FALSE,
      points = data.frame(
        d = c(14, 20, 20, 8, 11, 30), pwindow = c(1, 1, 3, 2, 1, 2),
        rho = c(0.2, 0.2, -0.1, 0.3, -0.1, 0.2)
      )
    ),
    list(
      cases = exptilt_large_shape_cases[1:2], vectorised = TRUE,
      points = data.frame(
        d = c(20, 12, 14), pwindow = c(1, 1, 2), rho = c(0.2, 0.2, -0.1)
      )
    )
  )
}

expect_fd_group <- function(model, group) {
  for (case in group$cases) {
    for (i in seq_len(nrow(group$points))) {
      point <- group$points[i, ]
      # The CDF underflows for large shapes at small delays
      underflow <- !group$vectorised && case$dist_id == 2L &&
        case$params[1] >= 100 && point$d < 0.5
      if (exptilt_case_ok(case, point$rho) && !underflow) {
        expect_fd_gradient(
          model, case, point, group$vectorised,
          if (is.null(point$scale)) 1 else point$scale
        )
      }
    }
  }
}

test_that("tilted log CDFs have finite gradients matching finite
  differences", {
  model <- exptilt_gradient_model()
  for (group in exptilt_fd_groups()) {
    expect_fd_group(model, group)
  }
  # The gamma upper tail far above the mean
  for (case in exptilt_stan_cases[c(3, 4)]) {
    for (d in c(40, 80, 200)) {
      res <- exptilt_gradient_at(model, case, d, 2, 0.3)
      label <- exptilt_case_label(case, d = d)
      expect_false(res$gradient_not_finite, info = label)
      expect_false(res$rejected, info = label)
      expect_true(all(is.finite(res$gradient)), info = label)
    }
  }
})

exptilt_gamma_case <- function(params) list(dist_id = 2L, params = params)

exptilt_normal_case <- function(params) list(dist_id = 18L, params = params)

# A point with relative tolerances for the delay parameters and the tilt
exptilt_ref_point <- function(case, d, pwindow, rho, vectorised = FALSE,
                              tol = rep(1e-5, 3)) {
  list(
    case = case, d = d, pwindow = pwindow, rho = rho,
    vectorised = vectorised, tol = tol
  )
}

exptilt_reference_points <- function() {
  large <- exptilt_large_shape_cases
  tilts <- expand.grid(
    pwindow = c(1, 7), scaled = c(1.1e-5, 1e-4, 1e-3, 3e-3),
    z = c(-6, -12), sign = c(-1, 1)
  )
  c(
    # Gamma lower tail of the vectorised PMF
    Map(
      function(params, pwindow) {
        exptilt_ref_point(exptilt_gamma_case(params), 12, pwindow, 0.3, TRUE)
      },
      params = list(c(20, 4), c(20, 4), c(100, 10), c(100, 10)),
      pwindow = c(1, 3, 1, 3)
    ),
    # Large shapes in the small tilt and direct forms
    Map(
      function(case, d, pwindow, rho) {
        exptilt_ref_point(large[[case]], d, pwindow, rho, tol = rep(3e-5, 3))
      },
      case = c(3, 3, 3, 1, 1, 1), d = c(5, 5.3, 4.6, 10, 11, 9),
      pwindow = c(2, 3, 1, 2, 1, 3),
      rho = c(-1e-5, -1e-3, 1e-5, -1e-5, 1e-5, -1e-3)
    ),
    # The tilt gradient near the small window and small delay limits
    Map(
      function(params, d, pwindow, rho) {
        exptilt_ref_point(
          exptilt_gamma_case(params), d, pwindow, rho,
          tol = c(Inf, Inf, 3e-5)
        )
      },
      params = list(
        c(100, 10), c(100, 10), c(2.5, 0.4), c(2.5, 0.4), c(20, 4),
        c(2.5, 0.4)
      ),
      d = c(3.744, 12, 6, 6, 5, 0.0009), pwindow = c(3, 3, 10, 2, 1, 2),
      rho = c(1e-5, 1e-5, 9.9e-6, -4.9e-5, 9e-5, 0.04)
    ),
    # The small window form in the vectorised PMF
    lapply(c(-3e-5, 3e-5, -1e-3, 4e-3), function(rho) {
      exptilt_ref_point(
        exptilt_gamma_case(c(3, 1)), 12, 2, rho, TRUE, c(1e-4, 1e-4, 1e-5)
      )
    }),
    # Normal delays far below the mean, across the small window limit
    lapply(c(1, 3, 1e-3), function(pwindow) {
      exptilt_ref_point(exptilt_normal_case(c(8, 3)), -96.5, pwindow, -3)
    }),
    Map(
      function(pwindow, scaled, z, sign) {
        exptilt_ref_point(
          exptilt_normal_case(c(0, 1)), z, pwindow, sign * scaled / pwindow
        )
      },
      tilts$pwindow, tilts$scaled, tilts$z, tilts$sign
    ),
    Map(
      function(rho, pwindow, z) {
        exptilt_ref_point(
          exptilt_normal_case(c(-4, 0.3)), -4 + 0.3 * z, pwindow, rho
        )
      },
      rho = c(2e-5, -1e-4, 3e-4), pwindow = c(7, 7, 1), z = c(-30, -12, -12)
    )
  )
}

# nolint start: object_usage_linter.
expect_reference_gradient <- function(model, point) {
  case <- point$case
  expected <- exptilt_log_gradient(
    case$dist_id, case$params, point$d, point$pwindow, point$rho,
    point$vectorised
  )
  res <- exptilt_gradient_at(
    model, case, point$d, point$pwindow, point$rho, point$vectorised
  )
  label <- exptilt_case_label(
    case, d = point$d, pwindow = point$pwindow, r = point$rho,
    vectorised = point$vectorised
  )
  expect_false(res$gradient_not_finite, info = label)
  expect_true(all(is.finite(res$gradient)), info = label)
  expect_true(
    all(abs(res$gradient - expected) <= point$tol * pmax(abs(expected), 1e-2)),
    info = paste0(
      label, ": gradient ", toString(signif(res$gradient, 6)),
      ", reference ", toString(signif(expected, 6))
    )
  )
}
# nolint end

test_that("tilted log CDF gradients match differences of the reference
  integral", {
  model <- exptilt_gradient_model()
  for (point in exptilt_reference_points()) {
    expect_reference_gradient(model, point)
  }
})

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
  # Log of the upper tail 1 - F_rho(d) of a gamma delay
  log_upper <- function(d, rate) {
    log_integrand <- function(z) {
      stats::pgamma(d - z, 200, rate, lower.tail = FALSE, log.p = TRUE) +
        log(exptilt_window_density(z, 1, 0.2))
    }
    shift <- log_integrand(0)
    shift + log(stats::integrate(
      function(z) exp(vapply(z, log_integrand, numeric(1)) - shift),
      0, 1, rel.tol = 1e-13, abs.tol = 0
    )$value)
  }
  for (d in c(12, 14, 16, 18, 20)) {
    res <- exptilt_unary_gradient(model, d, 20)
    h <- 1e-4
    slope <- (
      -log_upper(d, 20 * exp(2 * h)) + 8 * log_upper(d, 20 * exp(h)) -
        8 * log_upper(d, 20 * exp(-h)) + log_upper(d, 20 * exp(-2 * h))
    ) / (12 * h * 20)
    upper <- exp(log_upper(d, 20))
    expect_equal(
      res$gradient, -upper * slope / (1 - upper),
      tolerance = 1e-4, info = as.character(d)
    )
  }
})

test_that("primarycensored_log_gamma_pq gives both tails in value and shape
  gradient for any shape", {
  fractions <- c(0.001, 0.3, 0.7, 0.9, 1, 1.05, 1.2, 1.6, 3, 8)
  for (shape in c(0.3, 2.5, 20, 200, 1000, 1e4)) {
    x <- shape * fractions
    actual <- vapply(x, primarycensored_log_gamma_pq, numeric(2), shape)
    for (k in 1:2) {
      expected <- stats::pgamma(x, shape, lower.tail = k == 1, log.p = TRUE)
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
        res <- exptilt_unary_gradient(model, x, shape, k)
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
  expect_equal(
    vapply(z, primarycensored_log_std_normal_cdf, numeric(1)),
    stats::pnorm(z, log.p = TRUE),
    tolerance = 1e-13
  )
  for (zz in z) {
    res <- exptilt_unary_gradient(model, 0, zz)
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
