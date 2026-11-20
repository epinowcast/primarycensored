skip_on_cran()

# The C++ tilt transforms from pcd_stan_header() must return the same values
# as their Stan bodies and the same gradients as Stan's autodiff through
# those bodies, for the transform itself and for the exponentially tilted
# log CDFs and log PMFs that use it. See helper-stan-cpp.R for the shared
# helpers.

tilt_cpp_cases <- list(
  list(dist_id = 1L, params = list(c(1.5, 0.5), c(0.2, 1.2), c(3, 0.3))),
  list(dist_id = 4L, params = list(0.5, 2)),
  list(dist_id = 2L, params = list(
    c(2, 0.5), c(0.6, 1.3), c(20, 4), c(200, 20)
  )),
  list(dist_id = 18L, params = list(c(3, 2), c(-1, 3), c(8, 1))),
  list(dist_id = 3L, params = list(c(1.5, 5), c(0.7, 2), c(3, 30))),
  list(dist_id = 5L, params = list(c(1.5, 5, 2), c(0.8, 3, 0.6)))
)

# Tilts of both signs. With |rho| = 0.004 the small window form of
# primarycensored_exptilt_form() covers the first delays at pwindow 0.3
# and 1 for delays on the non-negative reals, and the scalar probe's
# pwindow of 2.5 puts its delays below 2.5 on the small delay form. The
# other tilts and settings use the direct form.
tilt_cpp_rhos <- c(0.2, -0.15, 0.004, -0.004, 1)

# Probe model with theta = (unconstrained delay parameters, rho). The first
# parameter of the lognormal (1) and normal (18) is on the identity scale,
# the other delay parameters on the log scale. The target is the weighted
# sum of the finite elements of `out`, and `out` is also a generated
# quantity.
tilt_probe_code <- function(data_code, out_code) {
  paste0(
    "data {\n  int dist_id;\n  int n_params;\n  int N;\n", data_code,
    "\n  vector[N] w;\n}\n",
    "parameters {\n  vector[n_params + 1] theta;\n}\n",
    "transformed parameters {\n  array[n_params] real params;\n",
    "  for (j in 1:n_params) {\n",
    "    params[j] = (dist_id == 1 || dist_id == 18) && j == 1 ? theta[j]\n",
    "      : exp(theta[j]);\n",
    "  }\n  real rho = theta[n_params + 1];\n}\n",
    "model {\n  vector[N] out;\n", out_code, "\n",
    "  for (i in 1:N) {\n",
    "    if (!is_inf(out[i]) && !is_nan(out[i])) {\n",
    "      target += w[i] * out[i];\n    }\n  }\n}\n",
    "generated quantities {\n  vector[N] out;\n", out_code, "\n}\n"
  )
}

tilt_theta <- function(dist_id, params, rho) {
  log_scale <- seq_along(params) > (dist_id %in% c(1L, 18L))
  c(ifelse(log_scale, log(abs(params)), params), rho)
}

tilt_rate <- function(dist_id, params) {
  switch(as.character(dist_id),
    "4" = params[1],
    "2" = params[2],
    Inf
  )
}

tilt_case_info <- function(dist_id, params, rho, setting) {
  sprintf(
    "dist %d, params (%s), rho %g, setting (%s)",
    dist_id, toString(params), rho, toString(setting)
  )
}

# The transform at -rho with the parameters, at 0 with the parameters (the
# delay CDF, with a data tilt), at -rho with the parameters as data and at
# -rho at a point t exp(rho / 10) that depends on rho, as the split point of
# the truncated logistic series does
pair_models <- function() {
  compile_cpp_and_stan( # nolint: object_usage_linter.
    tilt_probe_code(
      "  int M;\n  vector[M] t;\n  array[n_params] real fixed;",
      paste(
        "  for (i in 1:M) {",
        "    out[(12 * i - 11):(12 * i - 9)] = log_tilt_transform_pair(",
        "      t[i], dist_id, -rho, params);",
        "    out[(12 * i - 8):(12 * i - 6)] = log_tilt_transform_pair(",
        "      t[i], dist_id, 0, params);",
        "    out[(12 * i - 5):(12 * i - 3)] = log_tilt_transform_pair(",
        "      t[i], dist_id, -rho, fixed);",
        "    out[(12 * i - 2):(12 * i)] = log_tilt_transform_pair(",
        "      t[i] * exp(rho / 10), dist_id, -rho, params);",
        "  }",
        sep = "\n"
      )
    ),
    "pcdcpptiltpair"
  )
}

pair_points <- c(-2, 0, 0.5, 1, 3, 10, 40, 100)

pair_data <- function(dist_id, params, setting, w) {
  list(
    dist_id = dist_id, n_params = length(params), M = length(pair_points),
    t = pair_points, fixed = as.array(params),
    N = 12 * length(pair_points), w = w
  )
}

vectorized_models <- function() {
  compile_cpp_and_stan( # nolint: object_usage_linter.
    tilt_probe_code(
      "  int start;\n  int pw;",
      paste(
        "  out = primarycensored_exptilt_lcdf_vectorized(",
        "start, N, dist_id, params, pw, rho);"
      )
    ),
    "pcdcpptiltvec"
  )
}

# start, n, pwindow
vectorized_settings <- list(c(1, 30, 1), c(4, 40, 3))

vectorized_data <- function(dist_id, params, setting, w) {
  list(
    dist_id = dist_id, n_params = length(params), start = setting[1],
    N = setting[2], pw = setting[3], w = w
  )
}

# The scalar log CDF, whose q = d - pwindow is a local variable
scalar_models <- function() {
  compile_cpp_and_stan( # nolint: object_usage_linter.
    tilt_probe_code(
      "  vector[N] d;\n  real pwindow;",
      paste(
        "  for (i in 1:N) out[i] = primarycensored_exptilt_lcdf(",
        "d[i] | dist_id, params, pwindow, rho);"
      )
    ),
    "pcdcpptiltscalar"
  )
}

scalar_delays <- c(0.5, 1, 2.5, 7, 15.3, 30)

# pwindow
scalar_settings <- list(1, 2.5, 0.3)

scalar_data <- function(dist_id, params, setting, w) {
  list(
    dist_id = dist_id, n_params = length(params), d = scalar_delays,
    N = length(scalar_delays), pwindow = setting, w = w
  )
}

pmf_models <- function() {
  compile_cpp_and_stan( # nolint: object_usage_linter.
    tilt_probe_code(
      "  int max_delay;\n  real L;\n  real D;\n  real pwindow;",
      paste(
        "  out = primarycensored_sone_lpmf_vectorized(",
        "max_delay, L, D, dist_id, params, pwindow, 2, {rho});"
      )
    ),
    "pcdcpptiltpmf"
  )
}

# max_delay, L, D, pwindow
pmf_settings <- list(c(20, 0, 21, 1), c(20, 0, Inf, 2), c(15, 2, 30, 1))

pmf_data <- function(dist_id, params, setting, w) {
  list(
    dist_id = dist_id, n_params = length(params), max_delay = setting[1],
    L = setting[2], D = setting[3], pwindow = setting[4],
    N = setting[1] + 1, w = w
  )
}

tilt_probes <- list(
  pair = list(
    models = pair_models, settings = list(0), data = pair_data,
    n = function(setting) 12 * length(pair_points), upper_tail = FALSE
  ),
  vectorized = list(
    models = vectorized_models, settings = vectorized_settings,
    data = vectorized_data, n = function(setting) setting[2],
    upper_tail = FALSE
  ),
  scalar = list(
    models = scalar_models, settings = scalar_settings, data = scalar_data,
    n = function(setting) length(scalar_delays), upper_tail = FALSE
  ),
  pmf = list(
    models = pmf_models, settings = pmf_settings, data = pmf_data,
    n = function(setting) setting[1] + 1, upper_tail = TRUE
  )
)

# Calls `check(dist_id, params, rho, setting)` for every case, admissible
# tilt and setting of `probe`.
for_each_tilt_case <- function(probe, check) {
  for (case in tilt_cpp_cases) {
    for (params in case$params) {
      for (rho in tilt_cpp_rhos) {
        if (tilt_rate(case$dist_id, params) + rho <= 0) next
        for (setting in probe$settings) {
          check(case$dist_id, params, rho, setting)
        }
      }
    }
  }
}

for (probe_name in names(tilt_probes)) {
  probe <- tilt_probes[[probe_name]]

  test_that(paste("C++ tilt", probe_name, "values match the Stan body"), {
    models <- probe$models()
    for_each_tilt_case(probe, function(dist_id, params, rho, setting) {
      probe_data <- probe$data(
        dist_id, params, setting, rep(1, probe$n(setting))
      )
      expect_cpp_values_match( # nolint: object_usage_linter.
        models, probe_data, params,
        tilt_case_info(dist_id, params, rho, setting),
        theta = tilt_theta(dist_id, params, rho)
      )
    })
  })

  test_that(paste("C++ tilt", probe_name, "gradient matches Stan autodiff"), {
    skip_if_no_log_prob() # nolint: object_usage_linter.
    models <- probe$models()
    set.seed(1)
    for_each_tilt_case(probe, function(dist_id, params, rho, setting) {
      probe_data <- probe$data(
        dist_id, params, setting, stats::runif(probe$n(setting), 0.5, 1.5)
      )
      expect_cpp_gradient_matches( # nolint: object_usage_linter.
        models, probe_data, params,
        tilt_case_info(dist_id, params, rho, setting),
        upper_tail = probe$upper_tail,
        theta = tilt_theta(dist_id, params, rho)
      )
    })
  })
}

test_that("C++ tilt gradients match finite differences", {
  skip_if_no_log_prob() # nolint: object_usage_linter.
  set.seed(2)
  h <- 1e-5
  models <- vectorized_models()
  setting <- vectorized_settings[[1]]
  for (case in tilt_cpp_cases) {
    params <- case$params[[1]]
    data <- vectorized_data(
      case$dist_id, params, setting, stats::runif(setting[2], 0.5, 1.5)
    )
    theta <- tilt_theta(case$dist_id, params, 0.2)
    cpp <- stan_log_prob(models$cpp, data, list(theta = theta)) # nolint
    fd <- vapply(seq_along(theta), function(i) {
      shift <- replace(numeric(length(theta)), i, h)
      (stan_log_prob(models$cpp, data, list(theta = theta + shift))$lp - # nolint
        stan_log_prob(models$cpp, data, list(theta = theta - shift))$lp) / # nolint
        (2 * h)
    }, numeric(1))
    expect_equal(cpp$grad, fd, tolerance = 1e-6, info = toString(params))
  }
})

test_that("C++ tilt transform rejects an invalid shape as Stan does", {
  skip_if_no_log_prob() # nolint: object_usage_linter.
  models <- pair_models()
  data <- pair_data(2L, c(2, 0.5), 0, rep(1, 12 * length(pair_points)))
  upars <- list(theta = c(800, log(0.5), 0.1))
  for (model in models) {
    err <- tryCatch(
      stan_log_prob(model, data, upars), # nolint: object_usage_linter.
      error = conditionMessage
    )
    expect_type(err, "character")
    expect_match(err, "shape must be positive finite")
  }
})
