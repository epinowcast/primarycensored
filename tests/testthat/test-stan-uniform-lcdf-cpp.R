skip_on_cran()

# The C++ functions from pcd_stan_header() must return the same values as
# their Stan bodies and the same gradients as Stan's autodiff through those
# bodies. See helper-stan-cpp.R for the probe models.

cpp_cases <- list(
  list(dist_id = 1L, params = list(
    c(1.6, 0.42), c(1.5, 0.5), c(0.2, 1.2), c(-1, 0.3), c(3.5, 0.05),
    c(4, 2)
  )),
  list(dist_id = 2L, params = list(
    c(1.4, 0.38), c(2, 0.5), c(0.5, 0.2), c(20, 4), c(1.2, 0.02),
    c(50, 0.5), c(5, 1.2)
  )),
  list(dist_id = 3L, params = list(
    c(1.5, 5), c(0.7, 2), c(3, 30), c(5, 0.8), c(0.8, 3)
  ))
)

case_info <- function(dist_id, params, setting) {
  sprintf(
    "dist %d, params (%s), setting (%s)",
    dist_id, toString(params), toString(setting)
  )
}

vectorized_models <- function() {
  compile_cpp_and_stan( # nolint: object_usage_linter.
    cpp_probe_code( # nolint: object_usage_linter.
      "  int start;",
      paste(
        "  out = primarycensored_uniform_lcdf_vectorized(",
        "start, N, dist_id, params, pwindow);"
      )
    ),
    "pcdcppvec",
    fallback = TRUE
  )
}

# start, n, pwindow
vectorized_settings <- list(c(1, 40, 1), c(5, 40, 7), c(1, 70, 3))

vectorized_data <- function(dist_id, setting, w) {
  list(
    dist_id = dist_id, start = setting[1], N = setting[2],
    pwindow = setting[3], w = w
  )
}

scalar_models <- function() {
  compile_cpp_and_stan( # nolint: object_usage_linter.
    cpp_probe_code( # nolint: object_usage_linter.
      "  vector[N] d;",
      paste(
        "  for (i in 1:N) out[i] = primarycensored_uniform_lcdf(",
        "d[i] | dist_id, params, pwindow);"
      )
    ),
    "pcdcppscalar"
  )
}

scalar_delays <- c(0, 0.5, 1, 1.5, 2.5, 4, 7, 10.25, 15.3, 25, 40, 64.5)

# pwindow
scalar_settings <- list(1, 1.5, 7, 0.3)

scalar_data <- function(dist_id, setting, w) {
  list(
    dist_id = dist_id, pwindow = setting, d = scalar_delays,
    N = length(scalar_delays), w = w
  )
}

# The vectorised log PMF, which uses the C++ vectorised log CDF
pmf_models <- function() {
  compile_cpp_and_stan( # nolint: object_usage_linter.
    cpp_probe_code( # nolint: object_usage_linter.
      "  int max_delay;\n  real L;\n  real D;",
      paste(
        "  out = primarycensored_sone_lpmf_vectorized(",
        "max_delay, L, D, dist_id, params, pwindow, 1, rep_array(0.0, 0));"
      )
    ),
    "pcdcpppmf"
  )
}

# max_delay, L, D, pwindow
pmf_settings <- list(
  c(14, 0, 15, 1), c(40, 0, 41, 1), c(20, -Inf, Inf, 1), c(20, 0, Inf, 3),
  c(20, 2, 30.5, 1), c(20, 2.5, Inf, 1), c(20, -1.5, 21, 7), c(30, 3, 31, 2)
)

# Whether L is so far into the upper tail that 1 - F(L) < 1e-4, using
# 1 - F(L) <= S_T(L - pwindow). The truncation normaliser is then a
# difference of log CDFs close to 0. It amplifies differences of about
# 1e-13 in the special functions by 1 / (1 - F(L)), past the gradient
# tolerance of 1e-8.
far_upper_L <- function(dist_id, params, setting) {
  at <- setting[2] - setting[4]
  if (!is.finite(at) || at <= 0) {
    return(FALSE)
  }
  surv <- switch(dist_id,
    stats::plnorm(at, params[1], params[2], lower.tail = FALSE),
    stats::pgamma(at, params[1], params[2], lower.tail = FALSE),
    stats::pweibull(at, params[1], params[2], lower.tail = FALSE)
  )
  surv < 1e-4
}

pmf_data <- function(dist_id, setting, w) {
  list(
    dist_id = dist_id, max_delay = setting[1], L = setting[2],
    D = setting[3], pwindow = setting[4], N = setting[1] + 1, w = w
  )
}

no_skip <- function(dist_id, params, setting) FALSE

probes <- list(
  vectorized = list(
    models = vectorized_models, settings = vectorized_settings,
    data = vectorized_data, n = function(setting) setting[2],
    skip = no_skip, upper_tail = FALSE
  ),
  scalar = list(
    models = scalar_models, settings = scalar_settings,
    data = scalar_data, n = function(setting) length(scalar_delays),
    skip = no_skip, upper_tail = FALSE
  ),
  pmf = list(
    models = pmf_models, settings = pmf_settings,
    data = pmf_data, n = function(setting) setting[1] + 1,
    skip = far_upper_L, upper_tail = TRUE
  )
)

# Calls `check(dist_id, params, setting)` for every case and setting of
# `probe` that it does not skip.
for_each_case <- function(probe, check, skip = no_skip) {
  for (case in cpp_cases) {
    for (params in case$params) {
      for (setting in probe$settings) {
        if (!skip(case$dist_id, params, setting)) {
          check(case$dist_id, params, setting)
        }
      }
    }
  }
}

for (probe_name in names(probes)) {
  probe <- probes[[probe_name]]

  test_that(paste("C++", probe_name, "values match the Stan body"), {
    models <- probe$models()
    for_each_case(probe, function(dist_id, params, setting) {
      probe_data <- probe$data(dist_id, setting, rep(1, probe$n(setting)))
      expect_cpp_values_match(
        models, probe_data, params, case_info(dist_id, params, setting)
      )
    })
  })

  test_that(paste("C++", probe_name, "gradient matches Stan autodiff"), {
    skip_if_no_log_prob()
    models <- probe$models()
    set.seed(1)
    for_each_case(probe, function(dist_id, params, setting) {
      probe_data <- probe$data(
        dist_id, setting, stats::runif(probe$n(setting), 0.5, 1.5)
      )
      expect_cpp_gradient_matches(
        models, probe_data, params, case_info(dist_id, params, setting),
        upper_tail = probe$upper_tail
      )
    }, skip = probe$skip)
  })
}

test_that("C++ vectorized gradient with the nested autodiff dP/da matches
  Stan", {
  skip_if_no_log_prob()
  models <- vectorized_models()
  set.seed(3)
  for (case in cpp_cases[2:3]) {
    for (params in case$params) {
      data <- vectorized_data(case$dist_id, c(1, 40, 1), stats::runif(40))
      expect_cpp_gradient_matches(
        models, data, params, case_info(case$dist_id, params, "fallback"),
        model = "fallback"
      )
    }
  }
})

test_that("C++ gradients match finite differences", {
  skip_if_no_log_prob()
  set.seed(2)
  h <- 1e-5
  for (probe in probes) {
    models <- probe$models()
    setting <- probe$settings[[1]]
    for (case in cpp_cases) {
      params <- case$params[[1]]
      data <- probe$data(
        case$dist_id, setting, stats::runif(probe$n(setting), 0.5, 1.5)
      )
      theta <- cpp_probe_theta(case$dist_id, params)
      cpp <- stan_log_prob(models$cpp, data, list(theta = theta))
      fd <- vapply(1:2, function(i) {
        shift <- replace(numeric(2), i, h)
        (stan_log_prob(models$cpp, data, list(theta = theta + shift))$lp -
          stan_log_prob(models$cpp, data, list(theta = theta - shift))$lp) /
          (2 * h)
      }, numeric(1))
      expect_equal(cpp$grad, fd, tolerance = 1e-6, info = toString(params))
    }
  }
})

# testdata/uniform-lpmf-reference.csv is computed at 50 digits by
# data-raw/uniform-lpmf-reference.py. Elements below -500, and past the mode
# below -10, are left out as in expect_cpp_gradient_matches(), as are
# elements whose reference is NA (not a number even at 50 digits).
test_that("log PMF and gradient with the C++ log CDF match a 50-digit
  reference", {
  skip_if_no_log_prob()
  models <- pmf_models()
  ref <- utils::read.csv(test_path("testdata", "uniform-lpmf-reference.csv"))
  for (case in split(ref, ref[c("dist_id", "p1", "p2", "max_delay")],
    drop = TRUE
  )) {
    params <- c(case$p1[1], case$p2[1])
    n <- case$max_delay[1] + 1
    data <- pmf_data(
      case$dist_id[1], c(case$max_delay[1], 0, n, 1), numeric(n)
    )
    theta <- list(theta = cpp_probe_theta(data$dist_id, params))
    past_mode <- seq_len(n) > which.max(case$lpmf)
    keep <- !is.na(case$lpmf) & !is.na(case$grad1) & !is.na(case$grad2) &
      case$lpmf > -500 & !(past_mode & case$lpmf < -10)
    cpp_lpmf <- stan_gq(models$cpp, data, theta, "out")
    cpp_grad <- t(vapply(which(keep), function(d) {
      data$w <- replace(numeric(n), d, 1)
      stan_log_prob(models$cpp, data, theta)$grad
    }, numeric(2)))
    info <- case_info(data$dist_id, params, n)
    expect_equal(
      cpp_lpmf[keep], case$lpmf[keep],
      tolerance = 1e-10, info = info
    )
    expect_equal(
      cpp_grad, as.matrix(case[keep, c("grad1", "grad2")]),
      tolerance = 1e-9, ignore_attr = TRUE, info = info
    )
  }
})

# Lognormal log CDF and its gradient with respect to (mu, log(sigma)) at d,
# by integrating the delay CDF and its derivatives over the primary window.
# This does not cancel in the lower tail.
lognormal_uniform_lcdf_ref <- function(d, mu, sigma, pwindow) {
  z <- function(t) (log(pmax(t, 1e-300)) - mu) / sigma
  cdf <- function(t) ifelse(t > 0, stats::plnorm(t, mu, sigma), 0)
  dmu <- function(t) ifelse(t > 0, -stats::dnorm(z(t)) / sigma, 0)
  dsigma <- function(t) ifelse(t > 0, -z(t) * stats::dnorm(z(t)), 0)
  over_window <- function(f) {
    stats::integrate(
      function(u) f(d - u), 0, pwindow,
      rel.tol = 1e-13, abs.tol = 0, subdivisions = 1000L
    )$value
  }
  den <- over_window(cdf)
  c(log(den / pwindow), over_window(dmu) / den, over_window(dsigma) / den)
}

test_that("C++ lognormal uniform lcdf gradient is accurate in the lower
  tail", {
  skip_if_no_log_prob()
  models <- vectorized_models()
  mu <- 3.5
  sigma <- 0.05
  # The log CDF runs from about -590 at d = 6 to about -2e-4 at d = 40
  for (d in c(6, 10, 20, 30, 40)) {
    data <- vectorized_data(1L, c(1, 40, 1), replace(numeric(40), d, 1))
    cpp <- stan_log_prob(
      models$cpp, data, list(theta = c(mu, log(sigma)))
    )
    ref <- lognormal_uniform_lcdf_ref(d, mu, sigma, 1)
    expect_equal(cpp$lp, ref[1], tolerance = 1e-10, info = paste("d =", d))
    expect_equal(cpp$grad, ref[2:3], tolerance = 1e-9, info = paste("d =", d))
  }
})

test_that("C++ functions reject invalid arguments as Stan does", {
  skip_if_no_log_prob()
  models <- vectorized_models()
  data <- vectorized_data(2L, c(1, 10, 1), rep(1, 10))
  upars <- list(theta = c(log(2), 800))
  for (model in models[c("stan", "cpp")]) {
    err <- tryCatch(
      stan_log_prob(model, data, upars),
      error = conditionMessage
    )
    expect_type(err, "character")
    expect_match(err, "rate must be positive finite")
  }
  pmf <- pmf_models()
  data <- pmf_data(1L, c(14, 0, 10, 1), rep(1, 15))
  upars <- list(theta = c(1.6, log(0.42)))
  for (model in pmf) {
    err <- tryCatch(
      stan_log_prob(model, data, upars),
      error = conditionMessage
    )
    expect_type(err, "character")
    expect_match(err, "D must be at least max_delay + 1", fixed = TRUE)
  }
})

test_that("C++ vectorized lcdf rejects a non-integer pwindow and bad start", {
  skip_if_no_log_prob()
  models <- vectorized_models()
  upars <- list(theta = c(1.6, log(0.42)))
  bad <- list(
    list(setting = c(1, 10, 1.5), message = "pwindow must be a positive"),
    list(setting = c(0, 10, 1), message = "need 1 <= start <= n")
  )
  for (case in bad) {
    data <- vectorized_data(1L, case$setting, rep(1, 10))
    err <- tryCatch(
      stan_log_prob(models$cpp, data, upars),
      error = conditionMessage
    )
    expect_type(err, "character")
    expect_match(err, case$message)
  }
})
