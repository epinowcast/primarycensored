skip_on_cran()

# Uniform primary analytical solutions in Stan for the Exponential (4),
# Beta (9), Chi-square (13), Inverse gamma (16), Normal (18), Inverse
# chi-square (19), Pareto (21) and Scaled inverse chi-square (22) delays.
#
# The numerical path integrates the delay CDF over the primary window, see
# primarycensored_ode(). Its integrand is the exposed Stan function
# dist_lcdf(), so the reference here integrates exp(dist_lcdf()) with
# stats::integrate at a tolerance near double precision (see
# helper-uniform-reference.R). The Stan ODE solver itself only reaches an
# absolute 1e-6, which is used where the ODE path is run (the fallback
# cases). Analytical solutions are compared at a relative 1e-8, looser than
# the R tests because the incomplete gamma and beta functions in Stan lose
# a few digits in the far lower tail.

# nolint start: object_usage_linter.
stan_reference <- function(dist_id, params, d, pwindow, kinks = numeric(0)) {
  reference_uniform_cdf(
    function(t) exp(dist_lcdf(t, params, dist_id)),
    d, pwindow, kinks
  )
}

support_lower <- function(dist_id) {
  if (dist_has_positive_support(dist_id) == 1L) 0 else -Inf
}

analytical_lcdf <- function(d, dist_id, params, pwindow) {
  vapply(d, function(di) {
    primarycensored_analytical_lcdf(
      di, dist_id, params, pwindow, support_lower(dist_id), Inf, 1L,
      numeric(0)
    )
  }, numeric(1))
}
# nolint end

delays_positive <- c(0.001, 0.01, 0.1, 0.5, 1, 2, 5, 10, 30, 100)

# Each case is a delay with its parameter sets, the delays at which to
# evaluate it, and the kinks of F_T in t
uniform_cases <- list(
  list(
    name = "exponential", dist_id = 4L,
    params = list(0.001, 0.05, 0.5, 3, 20),
    delays = delays_positive, kinks = 0
  ),
  list(
    name = "normal", dist_id = 18L,
    params = list(c(-3, 0.1), c(0, 1), c(2, 3), c(10, 0.3), c(10, 3)),
    delays = c(-40, -20, -10, -5, -2, -1, 0, 0.5, 2, 5, 10, 20, 40),
    kinks = numeric(0)
  ),
  list(
    name = "chi-square", dist_id = 13L,
    params = list(0.5, 1, 3, 10),
    delays = delays_positive, kinks = 0
  ),
  list(
    name = "beta", dist_id = 9L,
    params = list(c(0.5, 0.5), c(1, 1), c(2, 3), c(5, 10), c(0.5, 10)),
    delays = c(0.001, 0.1, 0.3, 0.5, 0.9, 1, 1.5, 2, 3, 10), kinks = c(0, 1)
  ),
  list(
    name = "inverse gamma", dist_id = 16L,
    params = list(c(1.01, 1), c(1.5, 0.1), c(3, 2), c(3, 50), c(10, 5)),
    delays = c(delays_positive, 1000), kinks = 0
  ),
  list(
    name = "inverse chi-square", dist_id = 19L,
    params = list(2.05, 3, 10),
    delays = c(delays_positive, 1000), kinks = 0
  ),
  list(
    name = "scaled inverse chi-square", dist_id = 22L,
    params = list(c(2.05, 0.3), c(3, 1), c(10, 4)),
    delays = c(delays_positive, 1000), kinks = 0
  )
)

test_that("check_for_uniform_terms and check_for_analytical cover the new
  delays for a uniform primary only", {
  for (dist_id in c(1L, 2L, 3L, 4L, 5L, 9L, 13L, 16L, 18L, 19L, 21L, 22L)) {
    expect_identical(check_for_uniform_terms(dist_id, 1L), 1L)
    expect_identical(check_for_analytical(dist_id, 1L), 1L)
    expect_identical(check_for_uniform_terms(dist_id, 2L), 0L)
    expect_identical(check_for_analytical(dist_id, 2L), 0L)
  }
  for (dist_id in c(12L, 15L, 17L, 20L, 23L, 24L, 25L)) {
    expect_identical(check_for_uniform_terms(dist_id, 1L), 0L)
    expect_identical(check_for_analytical(dist_id, 1L), 0L)
  }
})

test_that("check_uniform_terms_params applies the shape rules and
  check_for_analytical_params the dispatch", {
  # Inverse gamma needs shape > 1, inverse and scaled inverse chi-square
  # need nu > 2
  expect_identical(check_uniform_terms_params(16L, c(1.0001, 1)), 1L)
  expect_identical(check_uniform_terms_params(16L, c(1, 1)), 0L)
  expect_identical(check_uniform_terms_params(16L, c(0.5, 1)), 0L)
  expect_identical(check_uniform_terms_params(19L, 2.0001), 1L)
  expect_identical(check_uniform_terms_params(19L, 2), 0L)
  expect_identical(check_uniform_terms_params(22L, c(2.0001, 1)), 1L)
  expect_identical(check_uniform_terms_params(22L, c(1.5, 1)), 0L)
  # Every other delay is valid for all admissible parameters
  expect_identical(check_uniform_terms_params(4L, 0.1), 1L)
  expect_identical(check_uniform_terms_params(21L, c(1, 0.5)), 1L)
  expect_identical(check_uniform_terms_params(18L, c(0, 1)), 1L)

  expect_identical(check_for_analytical_params(16L, 1L, c(2, 1)), 1L)
  expect_identical(check_for_analytical_params(16L, 1L, c(0.8, 1)), 0L)
  expect_identical(check_for_analytical_params(16L, 2L, c(2, 1)), 0L)
  expect_identical(check_for_analytical_params(4L, 1L, 0.3), 1L)
})

test_that("check_for_analytical_vectorized covers the new delays for integer
  windows", {
  for (dist_id in c(4L, 9L, 13L, 16L, 18L, 19L, 21L, 22L)) {
    expect_identical(check_for_analytical_vectorized(dist_id, 1L, 3), 1L)
    expect_identical(check_for_analytical_vectorized(dist_id, 1L, 1.5), 0L)
    expect_identical(check_for_analytical_vectorized(dist_id, 2L, 3), 0L)
  }
})

test_that("uniform primary terms match their closed forms", {
  ts <- c(0.05, 0.4, 1, 3.5, 12)
  # Exponential, the terms are [log G(t), -Inf] where G is the antiderivative
  # of the CDF, t - 1 / rate + exp(-rate * t) / rate
  for (rate in c(0.01, 0.5, 4)) {
    for (t in ts) {
      terms <- primarycensored_uniform_terms(t, 4L, rate)
      expect_equal(
        exp(terms[[1]]), t - 1 / rate + exp(-rate * t) / rate,
        tolerance = 1e-9
      )
      expect_identical(terms[[2]], -Inf)
    }
  }
  # Normal, G = sigma * (z * Phi(z) + phi(z)) = E[(t - T)^+]
  for (t in c(-6, -1, 0, 0.7, 4)) {
    terms <- primarycensored_uniform_terms(t, 18L, c(1, 2))
    z <- (t - 1) / 2
    expect_equal(
      exp(terms[[1]]), 2 * (z * pnorm(z) + dnorm(z)),
      tolerance = 1e-10
    )
    expect_identical(terms[[2]], -Inf)
  }
  # Inverse gamma, [log(t F(t)), log(E F~(t))] with E = beta / (alpha - 1)
  # and F~ the inverse gamma CDF with shape alpha - 1
  alpha <- 3
  beta <- 2
  for (t in ts) {
    terms <- primarycensored_uniform_terms(t, 16L, c(alpha, beta))
    expect_equal(
      exp(terms[[1]]), t * pgamma(beta / t, alpha, lower.tail = FALSE),
      tolerance = 1e-10
    )
    expect_equal(
      exp(terms[[2]]),
      beta / (alpha - 1) * pgamma(beta / t, alpha - 1, lower.tail = FALSE),
      tolerance = 1e-10
    )
  }
  # Beta, F = F~ = 1 for t >= 1
  terms <- primarycensored_uniform_terms(0.3, 9L, c(2, 3))
  expect_equal(exp(terms[[1]]), 0.3 * pbeta(0.3, 2, 3), tolerance = 1e-10)
  expect_equal(exp(terms[[2]]), 2 / 5 * pbeta(0.3, 3, 3), tolerance = 1e-10)
  terms <- primarycensored_uniform_terms(4, 9L, c(2, 3))
  expect_equal(exp(terms), c(4, 2 / 5), tolerance = 1e-12)
  # Pareto, G = int_{y_min}^t F, zero at and below y_min
  y_min <- 0.5
  a <- 2.5
  for (t in c(0.6, 1, 4, 30)) {
    terms <- primarycensored_uniform_terms(t, 21L, c(y_min, a))
    expect_equal(
      exp(terms[[1]]),
      stats::integrate(
        function(u) 1 - (y_min / u)^a, y_min, t,
        rel.tol = 1.2e-14
      )$value,
      tolerance = 1e-9
    )
  }
  expect_identical(
    primarycensored_uniform_terms(y_min, 21L, c(y_min, a)),
    c(-Inf, -Inf)
  )
  # Delays that map to another delay give the same terms
  for (t in ts) {
    expect_identical(
      primarycensored_uniform_terms(t, 13L, 7),
      primarycensored_uniform_terms(t, 2L, c(3.5, 0.5))
    )
    expect_identical(
      primarycensored_uniform_terms(t, 19L, 7),
      primarycensored_uniform_terms(t, 16L, c(3.5, 0.5))
    )
    expect_identical(
      primarycensored_uniform_terms(t, 22L, c(7, 1.5)),
      primarycensored_uniform_terms(t, 16L, c(3.5, 7 * 1.5^2 / 2))
    )
  }
})

test_that("uniform primary terms are -Inf for t <= 0 for non-negative delays", {
  params <- list(
    `4` = 0.5, `9` = c(2, 3), `13` = 3, `16` = c(3, 2), `19` = 5,
    `21` = c(0.5, 2), `22` = c(5, 1)
  )
  for (id in names(params)) {
    for (t in c(0, -0.5, -3)) {
      expect_identical(
        primarycensored_uniform_terms(t, as.integer(id), params[[id]]),
        c(-Inf, -Inf), info = paste("dist", id, "t", t)
      )
    }
  }
})

test_that("primarycensored_uniform_lower_bound clips only for non-negative
  delays", {
  expect_identical(primarycensored_uniform_lower_bound(0.5, 4L, 2), 0)
  expect_identical(primarycensored_uniform_lower_bound(3, 4L, 2), 1)
  expect_identical(primarycensored_uniform_lower_bound(0.5, 18L, 2), -1.5)
  expect_identical(primarycensored_uniform_lower_bound(-1, 18L, 2), -3)
})

test_that("the analytical log CDF matches the numerical path integrand for
  each new delay", {
  for (case in uniform_cases) {
    for (params in case$params) {
      for (pwindow in c(0.5, 1, 3, 10)) {
        info <- sprintf(
          "%s params = %s pwindow = %g", case$name, toString(params), pwindow
        )
        analytic <- exp(
          analytical_lcdf(case$delays, case$dist_id, params, pwindow)
        )
        reference <- stan_reference(
          case$dist_id, params, case$delays, pwindow, case$kinks
        )
        expect_rel_equal(analytic, reference, tolerance = 1e-8, info = info)
      }
    }
  }
})

test_that("primarycensored_lcdf and primarycensored_cdf dispatch to the
  analytical solution for the new delays", {
  for (case in uniform_cases) {
    params <- case$params[[2]]
    pwindow <- 2
    lower <- support_lower(case$dist_id)
    for (d in case$delays[case$delays > 0 & case$delays < 30]) {
      analytic <- primarycensored_analytical_lcdf(
        d, case$dist_id, params, pwindow, lower, Inf, 1L, numeric(0)
      )
      expect_identical(
        primarycensored_lcdf(
          d, case$dist_id, params, pwindow, lower, Inf, 1L, numeric(0)
        ),
        analytic, info = paste(case$name, "d", d)
      )
      expect_equal(
        primarycensored_cdf(
          d, case$dist_id, params, pwindow, lower, Inf, 1L, numeric(0)
        ),
        exp(analytic), tolerance = 1e-12, info = paste(case$name, "d", d)
      )
    }
  }
})

test_that("the Stan analytical solution matches the R solution for the
  exponential, normal, chi-square and beta", {
  cases <- list(
    list(
      dist_id = 4L, pdist = pexp, stan = 0.4, args = list(rate = 0.4),
      delays = c(0.2, 1, 3, 10, 40)
    ),
    list(
      dist_id = 18L, pdist = pnorm, stan = c(2, 1.5),
      args = list(mean = 2, sd = 1.5), delays = c(-8, -2, 0, 1, 4, 12)
    ),
    list(
      dist_id = 13L, pdist = pchisq, stan = 4, args = list(df = 4),
      delays = c(0.2, 1, 3, 10, 40)
    ),
    list(
      dist_id = 9L, pdist = pbeta, stan = c(2, 3),
      args = list(shape1 = 2, shape2 = 3), delays = c(0.1, 0.5, 0.9, 1.5, 3)
    )
  )
  for (case in cases) {
    for (pwindow in c(0.5, 1, 2, 5)) {
      obj <- do.call(new_pcens, c(list(case$pdist, dunif, list()), case$args))
      r_result <- pcens_cdf(obj, case$delays, pwindow)
      stan_result <- exp(
        analytical_lcdf(case$delays, case$dist_id, case$stan, pwindow)
      )
      expect_rel_equal(
        stan_result, r_result, tolerance = 1e-8,
        info = sprintf("dist %d pwindow %g", case$dist_id, pwindow)
      )
    }
  }
})

test_that("the exponential solution is the gamma solution with shape 1", {
  for (rate in c(0.05, 0.5, 4)) {
    for (pwindow in c(0.5, 2, 7)) {
      delays <- c(0.01, 0.5, 2, 10, 40)
      expect_equal(
        analytical_lcdf(delays, 4L, rate, pwindow),
        analytical_lcdf(delays, 2L, c(1, rate), pwindow),
        tolerance = 1e-9
      )
    }
  }
})

test_that("the chi-square and inverse chi-square solutions are the gamma and
  inverse gamma solutions", {
  delays <- c(0.01, 0.5, 2, 10, 40)
  for (pwindow in c(0.5, 2, 7)) {
    for (nu in c(1, 4, 12)) {
      expect_identical(
        analytical_lcdf(delays, 13L, nu, pwindow),
        analytical_lcdf(delays, 2L, c(nu / 2, 0.5), pwindow)
      )
    }
    for (nu in c(2.5, 4, 12)) {
      expect_identical(
        analytical_lcdf(delays, 19L, nu, pwindow),
        analytical_lcdf(delays, 16L, c(nu / 2, 0.5), pwindow)
      )
      expect_identical(
        analytical_lcdf(delays, 22L, c(nu, 1.3), pwindow),
        analytical_lcdf(delays, 16L, c(nu / 2, nu * 1.3^2 / 2), pwindow)
      )
    }
  }
})

test_that("Stan's chi-square and inverse chi-square CDFs are the gamma and
  inverse gamma CDFs the mapping relies on", {
  for (t in c(0.3, 1, 4, 20)) {
    expect_equal(
      dist_lcdf(t, 6, 13L), dist_lcdf(t, c(3, 0.5), 2L),
      tolerance = 1e-12
    )
    expect_equal(
      dist_lcdf(t, 6, 19L), dist_lcdf(t, c(3, 0.5), 16L),
      tolerance = 1e-12
    )
    expect_equal(
      dist_lcdf(t, c(6, 1.3), 22L),
      dist_lcdf(t, c(3, 6 * 1.3^2 / 2), 16L),
      tolerance = 1e-12
    )
  }
})

test_that("the normal solution handles a window that starts below zero", {
  # The primary window is [d - pwindow, d] and is not clipped at 0, so the
  # CDF is positive at d <= 0 and equals the Phi-based closed form
  params <- c(1, 2)
  pwindow <- 3
  delays <- c(-4, -1, 0, 0.5, 2)
  analytic <- exp(analytical_lcdf(delays, 18L, params, pwindow))
  expect_true(all(analytic > 0))
  expect_rel_equal(
    analytic,
    stan_reference(18L, params, delays, pwindow),
    tolerance = 1e-9
  )
  # The numerical and analytical paths give the same truncated CDF, with a
  # lower truncation point below 0
  for (L in c(-5, -1, 0, 2)) {
    for (D in c(6, 12, Inf)) {
      for (d in c(L + 0.5, 1, 3, 5)) {
        if (d <= L || d >= D) next
        truncated <- primarycensored_lcdf(
          d, 18L, params, pwindow, L, D, 1L, numeric(0)
        )
        cdf <- function(x) {
          stan_reference(18L, params, x, pwindow)
        }
        upper <- if (is.infinite(D)) 1 else cdf(D)
        expected <- (cdf(d) - cdf(L)) / (upper - cdf(L))
        expect_equal(
          exp(truncated), expected, tolerance = 1e-8,
          info = sprintf("L = %g, D = %g, d = %g", L, D, d)
        )
      }
    }
  }
})

test_that("the normal solution is accurate in the lower tail and across the
  switch to the asymptotic series", {
  # z = (t - mu) / sigma runs from well below the switch at -10 to above 0.
  # The lower tail CDFs are tiny, so the values are checked on the log scale
  # against the R solution, which is itself checked against quadrature in
  # test-pcens_cdf-uniform-solutions.R.
  mu <- 10
  sigma <- 0.5
  pwindow <- 1
  delays <- mu + sigma * c(-30, -20, -12, -10.5, -10, -9.5, -5, 0, 4)
  obj <- new_pcens(pnorm, dunif, list(), mean = mu, sd = sigma)
  r_result <- pcens_cdf(obj, delays, pwindow)
  stan_result <- exp(analytical_lcdf(delays, 18L, c(mu, sigma), pwindow))
  expect_rel_equal(stan_result, r_result, tolerance = 1e-9)
  # Deep in the tail the log CDF stays finite and decreases
  deep <- analytical_lcdf(mu + sigma * c(-60, -45, -40, -30), 18L,
    c(mu, sigma), pwindow
  )
  expect_true(all(is.finite(deep)))
  expect_true(all(diff(deep) > 0))
})

test_that("the solutions are -Inf below the support and 0 far above it", {
  for (id in c(4L, 13L)) {
    expect_identical(analytical_lcdf(0, id, 2, 1), -Inf)
    expect_identical(
      primarycensored_lcdf(-1, id, 2, 1, 0, Inf, 1L, numeric(0)),
      -Inf
    )
    # log(1 - tiny), limited by the rounding of the terms at d = 1e4
    expect_lt(abs(analytical_lcdf(1e4, id, 2, 1)), 1e-9)
  }
  # The Pareto has no mass below y_min
  expect_identical(analytical_lcdf(0.4, 21L, c(0.5, 2), 1), -Inf)
  expect_identical(analytical_lcdf(0.5, 21L, c(0.5, 2), 1), -Inf)
  expect_true(is.finite(analytical_lcdf(0.6, 21L, c(0.5, 2), 1)))
  # The Beta is 1 once the window is above 1
  expect_lt(abs(analytical_lcdf(3, 9L, c(2, 3), 1)), 1e-12)
})

test_that("the Pareto solution holds for shapes at and below 1 with no mean", {
  for (a in c(0.2, 0.9, 1, 1 + 1e-9, 1.1)) {
    for (y_min in c(0.05, 1)) {
      delays <- c(y_min * c(1.001, 1.5, 3), 10, 100)
      for (pwindow in c(0.5, 3)) {
        info <- sprintf(
          "y_min = %g, alpha = %g, pwindow = %g", y_min, a, pwindow
        )
        expect_rel_equal(
          exp(analytical_lcdf(delays, 21L, c(y_min, a), pwindow)),
          stan_reference(21L, c(y_min, a), delays, pwindow, y_min),
          tolerance = 1e-8, info = info
        )
      }
    }
  }
})

test_that("inverse gamma delays without a finite mean use the numerical
  path", {
  pwindow <- 2
  for (params in list(c(0.6, 1), c(1, 2), c(1, 0.5))) {
    expect_identical(check_for_analytical_params(16L, 1L, params), 0L)
    expect_error(
      primarycensored_uniform_terms(3, 16L, params),
      "shape > 1"
    )
    for (d in c(0.5, 2, 6, 20)) {
      lcdf <- primarycensored_lcdf(
        d, 16L, params, pwindow, 0, Inf, 1L, numeric(0)
      )
      # The ODE solver's absolute tolerance is 1e-6
      expect_lt(
        abs(exp(lcdf) - stan_reference(16L, params, d, pwindow, 0)), 1e-5
      )
    }
  }
  # The same holds for the chi-square variants with nu <= 2
  expect_identical(check_for_analytical_params(19L, 1L, 2), 0L)
  expect_identical(check_for_analytical_params(22L, 1L, c(1.5, 1)), 0L)
  expect_lt(
    abs(
      exp(primarycensored_lcdf(3, 19L, 1.5, pwindow, 0, Inf, 1L, numeric(0))) -
        stan_reference(19L, 1.5, 3, pwindow, 0)
    ),
    1e-5
  )
})

test_that("the inverse gamma solution stays finite and exact just above shape 1
  and deep in the lower tail", {
  # alpha - 1 = 1e-3 makes E = beta / (alpha - 1) large
  expect_rel_equal(
    exp(analytical_lcdf(c(1, 10, 100), 16L, c(1.001, 1), 2)),
    stan_reference(16L, c(1.001, 1), c(1, 10, 100), 2, 0),
    tolerance = 1e-8
  )
  # beta / t well above the point where the CDF underflows
  for (d in c(0.01, 0.05, 0.2)) {
    value <- analytical_lcdf(d, 16L, c(3, 50), 0.5)
    expect_true(value == -Inf || is.finite(value))
    expect_false(is.nan(value))
  }
})

test_that("the beta solution handles q below the support and d above it", {
  params <- c(2, 3)
  # q < 0 and d > 1 in the same window
  expect_rel_equal(
    exp(analytical_lcdf(c(0.2, 0.9, 1.2, 1.9), 9L, params, 2)),
    stan_reference(9L, params, c(0.2, 0.9, 1.2, 1.9), 2, c(0, 1)),
    tolerance = 1e-9
  )
  expect_equal(
    primarycensored_cdf(0.4, 9L, params, 0.5, 0, Inf, 1L, numeric(0)),
    stan_reference(9L, params, 0.4, 0.5, c(0, 1)),
    tolerance = 1e-9
  )
})

test_that("primarycensored_analytical_lcdf_vectorized matches the per-delay
  log CDF for the new delays", {
  per_delay <- function(delays, dist_id, params, pwindow) {
    vapply(delays, function(d) {
      primarycensored_lcdf(
        d, dist_id, params, pwindow, support_lower(dist_id), Inf, 1L,
        numeric(0)
      )
    }, numeric(1))
  }
  dists <- list(
    list(dist_id = 4L, params = list(0.05, 0.5, 3)),
    list(dist_id = 18L, params = list(c(10, 3), c(0, 1), c(-2, 0.4))),
    list(dist_id = 13L, params = list(3, 10)),
    list(dist_id = 9L, params = list(c(2, 3), c(0.5, 0.5))),
    list(dist_id = 16L, params = list(c(3, 2), c(1.5, 0.3))),
    list(dist_id = 19L, params = list(5)),
    list(dist_id = 21L, params = list(c(0.5, 2), c(2, 0.8))),
    list(dist_id = 22L, params = list(c(5, 1.5)))
  )
  n <- 25L
  for (dist in dists) {
    for (params in dist$params) {
      for (pwindow in c(1, 2, 3, 7)) {
        for (start in c(1L, 4L, 10L)) {
          info <- paste(
            "dist", dist$dist_id, "params", toString(params),
            "pwindow", pwindow, "start", start
          )
          vectorised <- primarycensored_analytical_lcdf_vectorized(
            start, n, dist$dist_id, params, pwindow
          )
          expect_length(vectorised, n)
          expect_identical(
            vectorised[start:n],
            per_delay(start:n, dist$dist_id, params, pwindow),
            info = info
          )
        }
      }
    }
  }
})

test_that("primarycensored_sone_lpmf_vectorized matches primarycensored_lpmf
  for the new delays, including the normal below 0", {
  dists <- list(
    list(dist_id = 4L, params = 0.3),
    list(dist_id = 18L, params = c(6, 2)),
    list(dist_id = 18L, params = c(0.5, 3)),
    list(dist_id = 13L, params = 6),
    list(dist_id = 9L, params = c(2, 3)),
    list(dist_id = 16L, params = c(3, 5)),
    list(dist_id = 21L, params = c(0.5, 1.5)),
    list(dist_id = 22L, params = c(5, 2))
  )
  settings <- list(
    list(max_delay = 15, L = 0, D = 16),
    list(max_delay = 15, L = -Inf, D = 16),
    list(max_delay = 15, L = -3, D = Inf),
    list(max_delay = 15, L = 3, D = 30)
  )
  for (dist in dists) {
    lower <- support_lower(dist$dist_id)
    for (s in settings) {
      for (pwindow in c(1, 3)) {
        info <- paste(
          "dist", dist$dist_id, "params", toString(dist$params),
          "L", s$L, "D", s$D, "pwindow", pwindow
        )
        vectorised <- primarycensored_sone_lpmf_vectorized(
          s$max_delay, s$L, s$D, dist$dist_id, dist$params, pwindow, 1L,
          numeric(0)
        )
        delays <- 0:s$max_delay
        full <- delays >= s$L & delays + 1 > s$L
        per_delay <- vapply(delays[full], function(d) {
          primarycensored_lpmf(
            d, dist$dist_id, dist$params, pwindow, d + 1, s$L, s$D, 1L,
            numeric(0)
          )
        }, numeric(1))
        # The CDF difference cannot resolve a PMF below the rounding error
        # of the CDF near 1, where it can be NaN for any delay, so compare
        # where the PMF is above 1e-9
        resolved <- is.finite(per_delay) & per_delay > log(1e-9)
        if (!(dist$dist_id == 9L && s$L >= 3)) {
          expect_true(any(resolved), info = info)
        }
        expect_equal(
          vectorised[full][resolved], per_delay[resolved],
          tolerance = 1e-9, info = info
        )
      }
    }
  }
  # Without truncation, for non-negative delays the PMF sums to the mass on
  # [0, max_delay + 1]
  pmf <- exp(primarycensored_sone_lpmf_vectorized(
    20, 0, 21, 4L, 0.3, 3, 1L, numeric(0)
  ))
  expect_equal(sum(pmf), 1, tolerance = 1e-10)
})

test_that("primarycensored_sone_pmf_vectorized matches R dprimarycensored for
  the new delays", {
  cases <- list(
    list(dist_id = 4L, params = 0.3, pdist = pexp, args = list(rate = 0.3)),
    list(
      dist_id = 18L, params = c(6, 2), pdist = pnorm,
      args = list(mean = 6, sd = 2)
    ),
    list(dist_id = 13L, params = 6, pdist = pchisq, args = list(df = 6)),
    list(
      dist_id = 9L, params = c(2, 3), pdist = pbeta,
      args = list(shape1 = 2, shape2 = 3)
    )
  )
  for (case in cases) {
    # The Beta has no mass above 1 + pwindow, where the CDF difference is
    # at the rounding error of 1
    max_delay <- if (case$dist_id == 9L) 2 else 15
    for (pwindow in c(1, 3)) {
      for (D in c(max_delay + 1, Inf)) {
        info <- paste("dist", case$dist_id, "pwindow", pwindow, "D", D)
        stan_pmf <- primarycensored_sone_pmf_vectorized(
          max_delay, 0, D, case$dist_id, case$params, pwindow, 1L,
          numeric(0)
        )
        r_pmf <- do.call(dprimarycensored, c(
          list(
            0:max_delay, case$pdist,
            pwindow = pwindow, swindow = 1, L = 0, D = D
          ),
          case$args
        ))
        expect_equal(stan_pmf, r_pmf, tolerance = 1e-6, info = info)
      }
    }
  }
})

# Gradients are only observable from a compiled model. As in
# test-stan-lognormal-tail-gradient.R, this drives `diagnose test=gradient`
# and compares the autodiff gradient with finite differences.
uniform_gradient_model <- function() {
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
    "  int n_params;\n",
    "  real d;\n",
    "  real pwindow;\n",
    "  int dist_id;\n",
    "}\n",
    "parameters {\n",
    "  array[n_params] real params;\n",
    "}\n",
    "model {\n",
    "  target += primarycensored_lcdf(\n",
    "    d | dist_id, params, pwindow,\n",
    "    dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),\n",
    "    positive_infinity(), 1, {0.0}[1:0]\n",
    "  );\n",
    "}\n"
  )
  path <- file.path(tempdir(), "pcd_uniform_gradient.stan")
  writeLines(code, path)
  suppressMessages(suppressWarnings(cmdstanr::cmdstan_model(path)))
}

uniform_gradient_at <- function(model, dist_id, params, d, pwindow) {
  stan_gradient_at( # nolint: object_usage_linter.
    model,
    data = list(
      n_params = length(params), d = d, pwindow = pwindow, dist_id = dist_id
    ),
    init = list(params = as.array(params))
  )
}

test_that("the new analytical solutions have finite gradients that match
  finite differences", {
  model <- uniform_gradient_model()
  # Stan's gradient of the incomplete gamma and beta functions with respect
  # to a shape has a limited precision, so the inverse gamma, inverse
  # chi-square and beta cases use 1e-3 (5e-3 with the partial expectation
  # shape 0.2) and the others 1e-4
  cases <- list(
    list(name = "exponential", id = 4L, params = 0.5, tol = 1e-4,
      delays = c(0.01, 0.3, 1, 3, 40)),
    list(name = "exponential", id = 4L, params = 0.001, tol = 1e-4,
      delays = c(0.05, 5, 40)),
    list(name = "normal", id = 18L, params = c(2, 1), tol = 1e-4,
      delays = c(-12, -8, -3, 0, 1, 4)),
    list(name = "normal", id = 18L, params = c(-3, 0.5), tol = 1e-4,
      delays = c(-8, -4, -2.5, 0, 3)),
    list(name = "normal", id = 18L, params = c(10, 0.3), tol = 1e-4,
      delays = c(3, 6, 8, 10, 12)),
    list(name = "chi-square", id = 13L, params = 3, tol = 1e-4,
      delays = c(0.5, 3, 10, 20)),
    list(name = "chi-square", id = 13L, params = 10, tol = 1e-4,
      delays = c(3, 10, 40)),
    list(name = "beta", id = 9L, params = c(2, 3), tol = 1e-3,
      delays = c(0.05, 0.5, 0.95, 1.5)),
    list(name = "beta", id = 9L, params = c(0.6, 1.5), tol = 1e-3,
      delays = c(0.05, 0.5, 0.95, 1.5)),
    list(name = "inverse gamma", id = 16L, params = c(3, 2), tol = 1e-3,
      delays = c(0.3, 1, 5, 50)),
    list(name = "inverse gamma", id = 16L, params = c(1.2, 5), tol = 5e-3,
      delays = c(0.5, 5, 50)),
    list(name = "inverse chi-square", id = 19L, params = 5, tol = 1e-3,
      delays = c(0.3, 1, 5, 50)),
    list(name = "scaled inverse chi-square", id = 22L, params = c(5, 1.5),
      tol = 1e-3, delays = c(0.3, 1, 5, 50)),
    list(name = "pareto", id = 21L, params = c(0.5, 2), tol = 1e-4,
      delays = c(0.6, 1, 5, 50)),
    list(name = "pareto", id = 21L, params = c(0.5, 0.8), tol = 1e-4,
      delays = c(0.6, 1, 5, 50)),
    list(name = "pareto", id = 21L, params = c(0.5, 1), tol = 1e-4,
      delays = c(0.6, 1, 5, 50)),
    list(name = "pareto", id = 21L, params = c(0.05, 1.3), tol = 1e-4,
      delays = c(0.06, 0.5, 5))
  )
  for (case in cases) {
    for (d in case$delays) {
      for (pwindow in c(0.5, 1, 3)) {
        info <- sprintf(
          "%s params = %s d = %g pwindow = %g", case$name,
          toString(case$params), d, pwindow
        )
        res <- uniform_gradient_at(model, case$id, case$params, d, pwindow)
        expect_false(res$rejected, info = info)
        expect_false(res$gradient_not_finite, info = info)
        expect_length(res$gradient, length(case$params))
        expect_true(all(is.finite(res$gradient)), info = info)
        expect_equal(
          res$gradient, res$finite_diff,
          tolerance = case$tol, info = info
        )
      }
    }
  }
})

test_that("deep in the inverse gamma lower tail the result is usable or
  log(0), never a non-finite gradient", {
  model <- uniform_gradient_model()
  # beta / t is large enough that the upper incomplete gamma underflows
  for (d in c(0.05, 0.02, 0.01)) {
    res <- uniform_gradient_at(model, 16L, c(3, 50), d, 0.5)
    expect_false(res$gradient_not_finite, info = paste("d", d))
    if (!res$rejected) {
      expect_true(all(is.finite(res$gradient)), info = paste("d", d))
    }
  }
})
