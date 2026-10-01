skip_on_cran()
skip_if_not_installed("flexsurv")

# Stan solutions for Weibull (dist_id 3, params = [shape, scale]) and
# generalised gamma (dist_id 5, params = [shape, scale, k]) delays with an
# exponentially tilted primary (primary_id 2 with primary_params = r).

stacy_cases <- lapply(exptilt_stacy_families(), function(family) {
  p <- family$stacy
  weibull <- startsWith(family$label, "weibull")
  list(
    label = family$label, family = family,
    dist_id = if (weibull) 3L else 5L,
    params = if (weibull) c(p$shape, p$scale) else c(p$shape, p$scale, p$k)
  )
})

stacy_label <- function(case, ...) {
  paste0(
    case$label, ", ",
    paste(names(list(...)), unlist(list(...)), sep = " = ", collapse = ", ")
  )
}

stacy_lcdf <- function(d, case, pwindow, rho) {
  vapply(
    d, primarycensored_exptilt_lcdf, numeric(1),
    case$dist_id, case$params, pwindow, rho
  )
}

# The ODE path is used where the log CDF is exactly that of the ODE
stacy_uses_ode <- function(d, case, pwindow, rho) {
  ode <- vapply(
    d, primarycensored_numeric_cdf, numeric(1),
    case$dist_id, case$params, pwindow, 2L, rho
  )
  stacy_lcdf(d, case, pwindow, rho) == log(ode)
}

test_that("Weibull and generalised gamma use the exponentially tilted
  solution", {
  for (dist_id in c(3L, 5L)) {
    for (xi in c(-50, 0, 0.3, 50)) {
      expect_identical(check_for_tilt_transform(dist_id, xi, c(1.5, 5, 2)), 1L)
    }
    expect_identical(check_for_exptilt(dist_id, 2L), 1L)
    expect_identical(check_for_exptilt(dist_id, 1L), 0L)
    expect_identical(check_for_analytical_vectorized(dist_id, 2L, 3), 1L)
    expect_identical(check_for_analytical_vectorized(dist_id, 2L, 1.5), 0L)
  }
  expect_identical(
    check_for_analytical_params(5L, c(1.3, 4, 2.5), 2L, 7), 1L
  )
})

test_that("Stan Stacy transforms match the R transforms", {
  ts <- c(1e-3, 0.4, 1, 3.5, 12, 40)
  for (case in stacy_cases) {
    obj <- exptilt_object(case$family, 0.1)
    for (xi in c(0, -1, -0.25, -0.05, 0.1, 0.25)) {
      stan <- t(vapply(
        ts, log_tilt_transform_pair, numeric(3), case$dist_id, xi, case$params
      ))
      expected <- .pcens_tilt_transform(obj, ts, xi)
      info <- stacy_label(case, xi = xi)
      keep <- !is.na(stan[, 1]) & !is.na(expected)
      expect_true(any(keep), info = info)
      expect_equal(
        stan[keep, 1], c(expected)[keep],
        tolerance = 1e-9, info = info
      )
      expect_equal(
        stan[keep, 3], attr(expected, "log_loss")[keep],
        tolerance = 1e-9, info = info
      )
      # The upper transform is the survival function for xi = 0, else `nan`
      upper <- c(.pcens_tilt_transform(obj, ts, xi, upper = TRUE))
      expect_identical(is.nan(stan[, 2]), is.nan(upper), info = info)
      if (xi == 0) {
        expect_equal(
          stan[is.finite(stan[, 2]), 2], upper[is.finite(stan[, 2])],
          tolerance = 1e-9, info = info
        )
      }
    }
    # Below the support the lower transform is -Inf
    expect_identical(
      log_tilt_transform_pair(-1, case$dist_id, -0.1, case$params)[1], -Inf
    )
  }
})

test_that("the Stan series is nan where it cancels or does not converge", {
  lower <- function(t, xi, params) {
    log_tilt_transform_pair(t, 3L, xi, params)[1]
  }
  expect_true(is.nan(lower(30, -1, c(0.7, 5))))
  expect_true(is.nan(lower(1000, 0.5, c(1.5, 5))))
  expect_false(is.nan(lower(3, -1, c(0.7, 5))))
})

test_that("Stan Stacy moments match the R moments", {
  ts <- c(-1, 0, 1e-2, 0.3, 4, 20)
  for (case in stacy_cases) {
    obj <- exptilt_object(case$family, 0.2)
    actual <- t(vapply(
      ts, primarycensored_tilt_moments, numeric(3), case$dist_id, case$params
    ))
    expect_equal(
      actual, unname(.pcens_tilt_moments(obj, ts)),
      tolerance = 1e-9, info = case$label
    )
  }
})

test_that("primarycensored_exptilt_lcdf matches a reference integral and
  the R implementation", {
  # The small tilt forms have a truncation error of up to about 3e-8
  for (case in stacy_cases) {
    for (pwindow in c(0.5, 2, 7)) {
      d <- c(1e-6, 1e-3, 0.3 * pwindow, pwindow, pwindow + 1e-3, 3, 6, 12, 25)
      for (rho in c(-1, -0.3, -1e-5, 1e-5, 1e-3, 0.05, 0.3, 1)) {
        info <- stacy_label(case, pwindow = pwindow, r = rho)
        actual <- stacy_lcdf(d, case, pwindow, rho)
        ode <- stacy_uses_ode(d, case, pwindow, rho)
        expected <- exptilt_reference(
          d, pwindow, rho, exptilt_cdf(case$family)
        )
        expect_lt(
          max_rel_diff(exp(actual)[!ode], expected[!ode]), 1e-7,
          label = info
        )
        expect_equal(
          exp(actual)[!ode],
          pcens_cdf(exptilt_object(case$family, rho), d, pwindow)[!ode],
          tolerance = 1e-9, info = info
        )
        # The ODE has tolerances of 1e-6
        expect_lt(max(abs(exp(actual) - expected)[ode], 0), 1e-5, label = info)
      }
    }
  }
})

test_that("the ODE path is used only where the series is not reliable", {
  d <- c(0.3, 1, 2.5, 6, 12)
  for (case in stacy_cases[c(1, 4)]) {
    for (rho in c(-0.5, -0.2, 0.05, 0.2)) {
      expect_false(
        any(stacy_uses_ode(d, case, 3, rho)),
        info = stacy_label(case, r = rho)
      )
    }
  }
  expect_true(stacy_uses_ode(30, stacy_cases[[2]], 2, 1))
})

test_that("primarycensored_cdf matches the empirical CDF of
  rprimarycensored", {
  set.seed(2025)
  n <- 1e5
  d <- c(1, 2.5, 5, 8, 14)
  for (case in stacy_cases[c(1, 4)]) {
    for (rho in c(-0.5, 0.3)) {
      samples <- do.call(
        rprimarycensored,
        c(
          list(
            n, case$family$rdist, pwindow = 2, swindow = 0,
            rprimary = rexpgrowth, rprimary_args = list(r = rho)
          ),
          case$family$args
        )
      )
      expected <- vapply(
        d, primarycensored_cdf, numeric(1),
        case$dist_id, case$params, 2, 0, Inf, 2L, rho
      )
      empirical <- vapply(d, function(dd) mean(samples <= dd), numeric(1))
      body <- expected > 0.01 & expected < 0.99
      z <- (empirical - expected)[body] /
        sqrt(expected * (1 - expected) / n)[body]
      expect_lt(max(abs(z)), 4.5, label = stacy_label(case, r = rho))
    }
  }
})

test_that("truncation is normalised for the series solutions", {
  rho <- 0.3
  for (case in stacy_cases[c(1, 4)]) {
    ref <- function(x) {
      exptilt_reference(x, 2, rho, exptilt_cdf(case$family))
    }
    for (bounds in list(c(0, 9), c(1, Inf), c(0.5, 7))) {
      upper <- if (is.finite(bounds[2])) ref(bounds[2]) else 1
      x <- c(1.5, 3, 5)
      x <- x[x > bounds[1] & x < bounds[2]]
      expected <- (ref(x) - ref(bounds[1])) / (upper - ref(bounds[1]))
      actual <- vapply(
        x, primarycensored_cdf, numeric(1),
        case$dist_id, case$params, 2, bounds[1], bounds[2], 2L, rho
      )
      expect_equal(actual, expected, tolerance = 1e-7, info = case$label)
    }
  }
})

test_that("the vectorised log CDF and PMF match the per delay forms", {
  for (case in stacy_cases[c(1, 2, 4)]) {
    for (pwindow in c(1, 3)) {
      for (rho in c(-0.3, -2e-5, 1e-9, 2e-5, 0.4)) {
        info <- stacy_label(case, pwindow = pwindow, r = rho)
        vectorised <- primarycensored_analytical_lcdf_vectorized(
          1L, 25L, case$dist_id, case$params, pwindow, 2L, rho
        )
        per_delay <- stacy_lcdf(1:25, case, pwindow, rho)
        # ODE values are shifted to the series value before them
        ode <- stacy_uses_ode(1:25, case, pwindow, rho)
        expect_identical(vectorised[!ode], per_delay[!ode], info = info)
        expect_equal(
          exp(vectorised), exp(per_delay), tolerance = 1e-6, info = info
        )
        vectorised <- primarycensored_sone_lpmf_vectorized(
          10, 2, 21, case$dist_id, case$params, pwindow, 2L, rho
        )
        per_delay <- vapply(0:10, function(d) {
          primarycensored_lpmf(
            d, case$dist_id, case$params, pwindow, d + 1, 2, 21, 2L, rho
          )
        }, numeric(1))
        expect_equal(vectorised, per_delay, tolerance = 1e-10, info = info)
      }
    }
  }
})

test_that("the vectorised upper tail PMF matches a survival based
  reference", {
  x <- 0:40
  for (case in stacy_cases) {
    for (pwindow in c(1, 3)) {
      for (rho in c(0.3, 0.05, -0.3)) {
        expected <- exptilt_pmf_reference(case$family, x, pwindow, rho)
        actual <- exp(primarycensored_sone_lpmf_vectorized(
          max(x), 0, Inf, case$dist_id, case$params, pwindow, 2L, rho
        ))
        info <- stacy_label(case, pwindow = pwindow, r = rho)
        # The ODE is used where the terms of the direct form are large
        expect_lt(
          max_rel_diff(actual[expected > 1e-6], expected[expected > 1e-6]),
          1e-4, label = info
        )
        if (case$label %in% c("weibull 3 2", "gengamma 3 1 0.2")) {
          keep <- expected > 1e-8
          expect_lt(
            max_rel_diff(actual[keep], expected[keep]), 1e-6, label = info
          )
        }
      }
    }
  }
})

test_that("the vectorised PMF is accurate across the switch to the ODE", {
  x <- 0:40
  for (case in stacy_cases[c(2, 4, 5)]) {
    expect_true(any(stacy_uses_ode(1:41, case, 3, 0.3)), info = case$label)
    expect_false(all(stacy_uses_ode(1:41, case, 3, 0.3)), info = case$label)
    expected <- exptilt_pmf_reference(case$family, x, 3, 0.3)
    actual <- exp(primarycensored_sone_lpmf_vectorized(
      max(x), 0, Inf, case$dist_id, case$params, 3, 2L, 0.3
    ))
    keep <- expected > 1e-4
    expect_lt(
      max_rel_diff(actual[keep], expected[keep]), 1e-5, label = case$label
    )
  }
})

# Gradients are only observable from a compiled model, so this builds a
# minimal one whose target is the log CDF or the summed vectorised log PMF.
stacy_gradient_model <- function() {
  testthat::skip_if_not_installed("cmdstanr")
  testthat::skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)))
  code <- paste0(
    pcd_load_stan_functions(wrap_in_block = TRUE, write_to_file = FALSE),
    "\ndata {\n",
    "  int dist_id;\n",
    "  int n_params;\n",
    "  int vectorised;\n",
    "  real d;\n",
    "  real pwindow;\n",
    "}\n",
    "parameters {\n",
    "  real<lower=0> p1;\n",
    "  real<lower=0> p2;\n",
    "  real<lower=0> p3;\n",
    "  real rho;\n",
    "}\n",
    "model {\n",
    "  array[3] real all_params = {p1, p2, p3};\n",
    "  array[n_params] real params = all_params[1:n_params];\n",
    "  if (vectorised) {\n",
    "    target += sum(primarycensored_sone_lpmf_vectorized(\n",
    "      to_int(d), 0.0, positive_infinity(), dist_id, params, pwindow, 2,\n",
    "      {rho}\n",
    "    ));\n",
    "  } else {\n",
    "    target += primarycensored_lcdf(\n",
    "      d | dist_id, params, pwindow, 0.0, positive_infinity(), 2, {rho}\n",
    "    );\n",
    "  }\n",
    "}\n"
  )
  path <- file.path(tempdir(), "pcdstacygradient.stan")
  writeLines(code, path)
  suppressMessages(suppressWarnings(cmdstanr::cmdstan_model(path)))
}

expect_stacy_gradient <- function(model, case, d, pwindow, rho,
                                  vectorised = FALSE) {
  res <- stan_gradient_at( # nolint: object_usage_linter.
    model,
    data = list(
      dist_id = case$dist_id, n_params = length(case$params),
      vectorised = as.integer(vectorised), d = d, pwindow = pwindow
    ),
    init = list(
      p1 = case$params[1], p2 = case$params[2],
      p3 = if (length(case$params) > 2) case$params[3] else 1, rho = rho
    )
  )
  info <- stacy_label(case, d = d, pwindow = pwindow, r = rho)
  expect_false(res$gradient_not_finite, info = info)
  expect_false(res$rejected, info = info)
  expect_length(res$gradient, 4)
  expect_true(all(is.finite(res$gradient)), info = info)
  # Relative to the size of each component, with a floor for tiny ones
  expect_true(
    all(abs(res$gradient - res$finite_diff) <=
      1e-4 * pmax(abs(res$finite_diff), 1e-2)),
    info = paste0(
      info, ": gradient ", toString(signif(res$gradient, 5)),
      ", finite difference ", toString(signif(res$finite_diff, 5))
    )
  )
}

test_that("tilted log CDFs have gradients matching finite differences", {
  model <- stacy_gradient_model()
  points <- list(
    list(d = 0.3, pwindow = 2, rho = 0.4),
    list(d = 2.5, pwindow = 2, rho = -0.2),
    list(d = 6, pwindow = 3, rho = 0.3),
    list(d = 14, pwindow = 3, rho = -0.15),
    list(d = 2.5, pwindow = 2, rho = 3e-4),
    list(d = 12, pwindow = 1, rho = -2e-4),
    list(d = 4, pwindow = 2, rho = 1e-6),
    list(d = 40, pwindow = 3, rho = 0.1)
  )
  for (case in stacy_cases[c(1, 2, 4)]) {
    for (point in points) {
      expect_stacy_gradient(
        model, case, point$d, point$pwindow, point$rho
      )
    }
  }
})

# Gradient of the summed log PMF on the log scale of the delay parameters, by
# central differences of the reference PMF. The upper tail PMF is a difference
# of CDFs close to 1, so finite differences of the Stan density are noisy
stacy_pmf_gradient <- function(case, pwindow, rho, delays) {
  theta <- c(log(case$params), rho)
  n <- length(case$params)
  total <- function(theta) {
    delay <- case$family
    delay$args <- as.list(stats::setNames(
      exp(theta[seq_len(n)]), c("shape", "scale", "k")[seq_len(n)]
    ))
    sum(log(exptilt_pmf_reference(delay, delays, pwindow, theta[n + 1])))
  }
  vapply(seq_along(theta), function(i) {
    shift <- replace(numeric(length(theta)), i, 1e-5)
    (total(theta + shift) - total(theta - shift)) / 2e-5
  }, numeric(1))
}

test_that("the summed vectorised log PMF has the gradient of the reference
  PMF", {
  model <- stacy_gradient_model()
  points <- list(
    list(case = stacy_cases[[1]], pwindow = 1, rho = 0.3),
    list(case = stacy_cases[[1]], pwindow = 1, rho = -0.15),
    list(case = stacy_cases[[4]], pwindow = 1, rho = 0.2),
    list(case = stacy_cases[[4]], pwindow = 3, rho = 0.3)
  )
  for (point in points) {
    case <- point$case
    n <- length(case$params)
    res <- stan_gradient_at( # nolint: object_usage_linter.
      model,
      data = list(
        dist_id = case$dist_id, n_params = n, vectorised = 1L, d = 40,
        pwindow = point$pwindow
      ),
      init = list(
        p1 = case$params[1], p2 = case$params[2],
        p3 = if (n > 2) case$params[3] else 1, rho = point$rho
      )
    )
    # The Jacobian of the log transform adds 1 to each delay parameter
    actual <- res$gradient[c(seq_len(n), 4)] - c(rep(1, n), 0)
    expected <- stacy_pmf_gradient(case, point$pwindow, point$rho, 0:40)
    expect_false(res$rejected)
    expect_equal(
      actual, expected, tolerance = 1e-3,
      info = stacy_label(case, pwindow = point$pwindow, r = point$rho)
    )
  }
})
