skip_on_cran()

# Log-logistic delay (dist_id 31, params = [scale, shape]) with a uniform
# primary (primary_id 1).


test_that("the log-logistic delay is registered in Stan", {
  expect_identical(dist_has_positive_support(31L), 1L)
  x <- c(1e-6, 0.3, 1, 4, 30, 1e3)
  for (case in ll_cases) {
    expect_equal(
      vapply(x, dist_lcdf, numeric(1), case$params, 31L),
      stats::plogis(
        case$params[2] * (log(x) - log(case$params[1])),
        log.p = TRUE
      ),
      tolerance = 1e-13
    )
    expect_identical(dist_lcdf(0, case$params, 31L), -Inf)
    expect_identical(dist_lcdf(-2, case$params, 31L), -Inf)
  }
  # No underflow to a zero CDF far below the scale
  expect_equal(
    dist_lcdf(1e-200, c(2, 4.5), 31L),
    4.5 * (log(1e-200) - log(2)),
    tolerance = 1e-12
  )
})

test_that("dispatch checks include the log-logistic delay", {
  expect_identical(check_for_uniform_terms(31L, 1L), 1L)
  expect_identical(check_for_uniform_terms(31L, 2L), 0L)
  expect_identical(check_for_analytical(31L, 1L), 1L)
  expect_identical(check_for_analytical_vectorized(31L, 1L, 3), 1L)
  expect_identical(check_for_analytical_vectorized(31L, 1L, 2.5), 0L)
  expect_identical(
    check_for_analytical_params(31L, c(3, 0.5), 1L, numeric(0)), 1L
  )
  # A tiny shape has no solution at all
  expect_identical(
    check_for_analytical_params(31L, c(3, 0.005), 1L, numeric(0)), 0L
  )
  expect_identical(
    check_for_analytical_params(31L, c(3, 0.01), 1L, numeric(0)), 1L
  )
})

test_that("the Stan partial moment ratio matches the R series", {
  # Across the regimes (A up to 1, up to 3, and beyond), at the boundaries,
  # and for powers that are integers, where terms of the tail series need
  # their limit
  a <- c(0.05, 0.3, 0.5, 1, 1.0001, 2, 3.5, 7, 20, 60, 200)
  log_A <- c(
    -30, -8, -2, -0.5, -2e-3, -1e-3 + 1e-9, -5e-5, -1e-7, 0, 1e-7, 5e-5,
    1e-3 - 1e-9, 2e-3, 0.5, 1, log(3) + 1e-9, 2, 5, 6.9, log(1000), 7.2, 10, 20,
    60
  )
  expected <- .loglogistic_ratio(a, log_A)
  for (i in seq_along(log_A)) {
    for (j in seq_along(a)) {
      expect_equal(
        loglogistic_moment_ratio(a[j], log_A[i], 1e-16), expected[i, j],
        tolerance = 1e-10,
        info = sprintf("a = %g, log A = %g", a[j], log_A[i])
      )
    }
  }
})

test_that("a shape below 0.01 uses the ODE path", {
  tiny <- c(3, 0.005)
  expect_identical(
    primarycensored_lcdf(2, 31L, tiny, 2, 0, Inf, 1L, numeric(0)),
    log(ll_ode_cdf(ll_ode_model(), tiny, 2, 2, 1L, numeric(0)))
  )
})

test_that("the vectorised and single uniform paths use the same rule", {
  # Below the minimum shape both use the numerical path
  params <- c(2, 0.009)
  expect_identical(check_for_analytical_params(31L, params, 1L, numeric(0)), 0L)
  vec <- primarycensored_lcdf_vectorized(1L, 3L, 31L, params, 1, 1L, numeric(0))
  single <- vapply(
    1:3, primarycensored_lcdf, numeric(1),
    31L, params, 1, 0, Inf, 1L, numeric(0)
  )
  expect_identical(vec, single)
  lpmf <- primarycensored_sone_lpmf_vectorized(
    2L, 0, Inf, 31L, params, 1, 1L, numeric(0)
  )
  expect_equal(
    lpmf,
    vapply(0:2, function(d) {
      primarycensored_lpmf(
        d, 31L, params, 1, d + 1, 0, Inf, 1L, numeric(0)
      )
    }, numeric(1)),
    tolerance = 1e-12
  )
})

test_that("primary terms and CDFs match a reference integral", {
  ode_model <- ll_ode_model()
  for (case in ll_cases) {
    cdf <- ll_case_cdf(case)
    for (pwindow in c(0.5, 1, 2, 7)) {
      d <- sort(c(
        1e-6, 1e-3, 0.3 * pwindow, pwindow - 1e-3, pwindow, pwindow + 1e-3,
        2, 3, 6, 12, 25, 100
      ))
      expected <- vapply(d, function(dd) {
        stats::integrate(
          cdf, max(dd - pwindow, 0), dd,
          rel.tol = 1e-13, abs.tol = 0, subdivisions = 2000L
        )$value / pwindow
      }, numeric(1))
      info <- ll_case_label(case, pwindow = pwindow)
      lcdf <- vapply(
        d, primarycensored_lcdf, numeric(1),
        31L, case$params, pwindow, 0, Inf, 1L, numeric(0)
      )
      expect_lt(max_rel_diff(exp(lcdf), expected), 1e-8, label = info)
      direct <- vapply(
        d, primarycensored_analytical_lcdf, numeric(1),
        31L, case$params, pwindow, 0, Inf, 1L, numeric(0)
      )
      expect_identical(direct, lcdf)
      ode <- ll_ode_cdf(ode_model, case$params, d, pwindow, 1L, numeric(0))
      expect_lt(max(abs(exp(lcdf) - ode)), 1e-4, label = info)
    }
  }
})

test_that("the uniform terms are the endpoint terms of the analytic CDF", {
  params <- c(5, 2)
  terms_d <- primarycensored_uniform_terms(6, 31L, params)
  terms_q <- primarycensored_uniform_terms(4, 31L, params)
  expect_equal(
    primarycensored_uniform_lcdf_from_terms(terms_d, terms_q, 2),
    primarycensored_lcdf(6, 31L, params, 2, 0, Inf, 1L, numeric(0)),
    tolerance = 1e-14
  )
  expect_identical(
    primarycensored_uniform_terms(0, 31L, params), c(-Inf, -Inf)
  )
  expect_identical(
    primarycensored_uniform_terms(-1, 31L, params), c(-Inf, -Inf)
  )
})

test_that("the uniform CDF matches the R implementation", {
  d <- c(0.3, 1, 2.5, 6, 15, 1e3)
  for (case in ll_cases) {
    obj <- do.call(
      new_pcens,
      c(
        list(
          pdist = pllogis_test, dprimary = dunif, primary_args = list()
        ),
        case$args
      )
    )
    expect_equal(
      exp(vapply(
        d, primarycensored_lcdf, numeric(1),
        31L, case$params, 3, 0, Inf, 1L, numeric(0)
      )),
      pcens_cdf(obj, d, 3),
      tolerance = 1e-9
    )
  }
})

test_that("primarycensored_lcdf_vectorized shares terms", {
  params <- c(5, 2)
  vec <- primarycensored_lcdf_vectorized(
    1L, 20L, 31L, params, 3, 1L, numeric(0)
  )
  single <- vapply(
    1:20, primarycensored_lcdf, numeric(1),
    31L, params, 3, 0, Inf, 1L, numeric(0)
  )
  expect_equal(vec, single, tolerance = 1e-9)
})

test_that("the vectorised PMF matches the per delay PMF with truncation", {
  params <- c(5, 2)
  for (bounds in list(c(0, Inf), c(0, 25), c(2, 25))) {
    vec <- primarycensored_sone_lpmf_vectorized(
      15L, bounds[1], bounds[2], 31L, params, 3, 1L, numeric(0)
    )
    single <- vapply(1:15, function(d) {
      primarycensored_lpmf(
        d - 1L, 31L, params, 3, d, bounds[1], bounds[2], 1L, numeric(0)
      )
    }, numeric(1))
    expect_equal(vec[1:15], single, tolerance = 1e-8)
  }
})

test_that("log-logistic uniform log CDFs have finite gradients matching
  finite differences", {
  model <- ll_gradient_model()
  points <- list(
    list(d = 0.3, pwindow = 2),
    list(d = 1, pwindow = 1),
    list(d = 6, pwindow = 3),
    list(d = 40, pwindow = 3),
    list(d = 3000, pwindow = 2)
  )
  for (case in ll_cases) {
    for (point in points) {
      res <- expect_ll_gradient(model, case, point, 1L)
      # The tilt does not enter, so its gradient is 0
      expect_identical(unname(res$gradient[3]), 0)
    }
  }
})

test_that("the vectorised log-logistic log PMF has finite gradients
  matching finite differences", {
  model <- ll_gradient_model()
  points <- list(
    list(d = 12, pwindow = 3, rho = 0, primary_id = 1L),
    list(d = 20, pwindow = 4, rho = 0, primary_id = 1L)
  )
  for (case in ll_cases[c(2, 3, 4)]) {
    for (point in points) {
      expect_ll_gradient(model, case, point, point$primary_id, TRUE)
    }
  }
})

test_that("pcd_cmdstan_model rejects the reserved distribution identifiers", {
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)))
  counts <- data.frame(
    delay = 1, delay_upper = 2, pwindow = 2, relative_obs_time = 25, n = 3
  )
  stan_data <- pcd_as_stan_data(
    counts,
    dist_id = pcd_stan_dist_id("loglogistic"), primary_id = 1,
    param_bounds = list(lower = c(0, 0), upper = c(Inf, Inf)),
    primary_param_bounds = list(lower = numeric(0), upper = numeric(0)),
    priors = list(location = c(5, 2), scale = c(5, 2)),
    primary_priors = list(location = numeric(0), scale = numeric(0))
  )
  model <- suppressMessages(suppressWarnings(pcd_cmdstan_model()))
  for (reserved in c(29L, 30L)) {
    stan_data$dist_id <- reserved
    # The rejection in transformed data stops every chain
    expect_message(
      suppressWarnings(model$sample(
        data = stan_data, seed = 1, chains = 1, iter_warmup = 1,
        iter_sampling = 1, refresh = 0
      )),
      "is reserved and is not a distribution"
    )
  }
  stan_data$dist_id <- 31L
  fit <- suppressMessages(suppressWarnings(model$sample(
    data = stan_data, seed = 1, chains = 1, iter_warmup = 50,
    iter_sampling = 50, refresh = 0, show_messages = FALSE
  )))
  expect_s3_class(fit, "CmdStanMCMC")
})

test_that("the Stan uniform CDF is accurate where the difference cancels", {
  for (case in uniform_conditioning_cases) {
    ref <- loglogistic_censored_reference(case[4], case[1], case[2], 0, case[3])
    actual <- exp(primarycensored_lcdf(
      case[4], 31L, c(case[2], case[1]), case[3], 0, Inf, 1L, numeric(0)
    ))
    # The survival is resolved to 1e-16 in absolute terms
    expect_lt(
      abs(actual - ref[["cdf"]]) / max(min(ref), 1e-9), 1e-6,
      label = paste(case, collapse = " ")
    )
  }
})

test_that("the Stan uniform PMF is positive and accurate at large q", {
  for (case in uniform_pmf_cases) {
    params <- c(case[2], case[1])
    ref <- vapply(
      case[4] + 0:1, loglogistic_censored_reference, numeric(2),
      shape = case[1], scale = case[2], rho = 0, pwindow = case[3]
    )
    expected <- -diff(ref["survival", ])
    actual <- exp(primarycensored_lpmf(
      as.integer(case[4]), 31L, params, case[3], case[4] + 1, 0, Inf, 1L,
      numeric(0)
    ))
    expect_gt(actual, 0)
    expect_lt(
      abs(actual - expected) / expected, 1e-6,
      label = paste(case, collapse = " ")
    )
  }
})

test_that("the vectorised uniform CDF uses the numerical CDF where the
  difference cancels", {
  params <- c(5, 0.5)
  n <- 1e5L
  vec <- primarycensored_lcdf_vectorized(
    n - 2L, n, 31L, params, 1, 1L, numeric(0)
  )
  single <- vapply(
    (n - 2L):n, primarycensored_lcdf, numeric(1),
    31L, params, 1, 0, Inf, 1L, numeric(0)
  )
  expect_identical(vec[(n - 2L):n], single)
  expect_true(all(diff(single) > 0))
})

test_that("the uniform conditioning check flags the ill conditioned terms", {
  params <- c(5, 2)
  ill <- function(d) {
    terms_d <- primarycensored_uniform_terms(d, 31L, params)
    terms_q <- primarycensored_uniform_terms(max(d - 1, 0), 31L, params)
    loglogistic_window_ill_conditioned(
      terms_d[1], terms_q[1], d,
      1, primarycensored_uniform_lcdf_from_terms(terms_d, terms_q, 1)
    )
  }
  expect_identical(ill(10), 0L)
  expect_identical(ill(1e6), 1L)
  # The PMF is small next to the tail in a heavy tail, so the rule also
  # follows the density of the CDF
  params <- c(5, 0.5)
  terms_d <- primarycensored_uniform_terms(1e4, 31L, params)
  terms_q <- primarycensored_uniform_terms(1e4 - 1, 31L, params)
  expect_identical(
    loglogistic_window_ill_conditioned(
      terms_d[1], terms_q[1], 1e4, 1,
      primarycensored_uniform_lcdf_from_terms(terms_d, terms_q, 1)
    ),
    1L
  )
})

test_that("the Stan numerical CDF matches the R numerical CDF", {
  cases <- list(
    c(5, 2, 1, 4), c(5, 2, 3, 6), c(5, 0.5, 1, 1e4), c(1, 0.05, 2, 1e6),
    c(16.9, 115.5, 0.00179, 19.8), c(8, 4.5, 2, 30)
  )
  for (case in cases) {
    family <- list(
      pdist = pllogis_test, args = list(shape = case[2], scale = case[1])
    )
    r_cdf <- .loglogistic_numeric_cdf(uniform_object(family), case[4], case[3])
    stan_cdf <- exp(loglogistic_numeric_lcdf(
      case[4], c(case[1], case[2]), case[3]
    ))
    smaller <- min(r_cdf, 1 - r_cdf)
    expect_lt(
      abs(stan_cdf - r_cdf) / max(smaller, 1e-14), 1e-6,
      label = paste(case, collapse = " ")
    )
  }
})

test_that("the log-logistic fallback has gradients matching differences", {
  model <- ll_gradient_model()
  # Points where the difference cancels
  uniform <- list(
    list(case = list(params = c(5, 0.5)), d = 1e5, pwindow = 1),
    list(case = list(params = c(1, 0.05)), d = 1e6, pwindow = 2),
    list(case = list(params = c(16.9, 115.5)), d = 19.8, pwindow = 0.00179)
  )
  for (point in uniform) {
    expect_ll_gradient(model, point$case, point, 1L, tolerance = 1e-3)
  }
})

test_that("Stan log-logistic CDFs match rprimarycensored samples", {
  withr::local_seed(3)
  n <- 20000
  pwindow <- 2
  for (case in ll_cases[c(1, 3, 5)]) {
    samples <- loglogistic_samples(n, ll_case_family(case), pwindow, 0)
    qs <- unname(stats::quantile(samples, c(0.05, 0.25, 0.5, 0.75, 0.95)))
    stan_cdf <- exp(vapply(
      qs, primarycensored_lcdf, numeric(1), 31L, case$params, pwindow, 0,
      Inf, 1L, numeric(0)
    ))
    expect_lt(
      max(abs(vapply(qs, function(q) mean(samples <= q), numeric(1)) -
        stan_cdf)),
      0.015,
      label = ll_case_label(case, pwindow = pwindow)
    )
  }
})
