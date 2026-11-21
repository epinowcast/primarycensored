skip_on_cran()

# The truncated logistic (primary_id 3) and truncated Gumbel (primary_id 4)
# solutions sum tilt transforms from log_tilt_transform_pair(), which uses
# the C++ primarycensored_tilt_pair() with the header. Their log CDFs and
# log PMFs with the C++ functions (cpp = TRUE) must match the Stan bodies
# (cpp = FALSE) in value and in the gradient, including the gradient in the
# primary parameters. See helper-stan-cpp.R for the shared helpers.

# primary_id, delay cases and primary parameter cases. The locations are
# below, inside and above the window. Inside the window a location that is
# not an integer puts the split point d - location of the truncated
# logistic series at a point that depends on the location.
primary_cpp_cases <- list(
  list(
    primary_id = 3L,
    delays = list(
      list(dist_id = 4L, params = 2), list(dist_id = 4L, params = 25),
      list(dist_id = 2L, params = c(2, 0.5)),
      list(dist_id = 2L, params = c(20, 40)),
      list(dist_id = 18L, params = c(3, 2)),
      list(dist_id = 18L, params = c(-1, 3))
    ),
    primary = list(c(-0.5, 0.3), c(0.35, 0.3), c(0.5, 0.05), c(1.6, 0.4))
  ),
  list(
    primary_id = 4L,
    delays = list(
      list(dist_id = 18L, params = c(3, 2)),
      list(dist_id = 18L, params = c(8, 1)),
      list(dist_id = 18L, params = c(-1, 3))
    ),
    primary = list(c(-0.5, 0.5), c(0.4, 0.3), c(1.5, 0.5), c(3, 0.2))
  )
)

# Whether the series applies, from the Stan functions exposed in setup.R.
# Otherwise the numerical path is used.
primary_series_applies <- function(primary_id, dist_id, params, primary,
                                   pwindow) {
  if (primary_id == 3L) {
    check_for_tlogis_params( # nolint: object_usage_linter.
      dist_id, params, primary, pwindow
    ) == 1
  } else {
    check_for_gumbel_params( # nolint: object_usage_linter.
      dist_id, params, primary[1], primary[2]
    ) == 1
  }
}

# Probe model with theta = (unconstrained delay parameters, location,
# log scale). The first parameter of the normal (18) is on the identity
# scale, the other delay parameters on the log scale. The target is the
# weighted sum of the finite elements of `out`, and `out` is also a
# generated quantity.
primary_probe_code <- function(data_code, out_code) {
  paste0(
    "data {\n  int dist_id;\n  int n_params;\n  int primary_id;\n",
    "  real pwindow;\n  int N;\n", data_code, "\n  vector[N] w;\n}\n",
    "parameters {\n  vector[n_params + 2] theta;\n}\n",
    "transformed parameters {\n  array[n_params] real params;\n",
    "  for (j in 1:n_params) {\n",
    "    params[j] = dist_id == 18 && j == 1 ? theta[j] : exp(theta[j]);\n",
    "  }\n",
    "  array[2] real pp = {theta[n_params + 1], exp(theta[n_params + 2])};\n",
    "}\n",
    "model {\n  vector[N] out;\n", out_code, "\n",
    "  for (i in 1:N) {\n",
    "    if (!is_inf(out[i]) && !is_nan(out[i])) {\n",
    "      target += w[i] * out[i];\n    }\n  }\n}\n",
    "generated quantities {\n  vector[N] out;\n", out_code, "\n}\n"
  )
}

primary_theta <- function(dist_id, params, primary) {
  delay <- if (dist_id == 18L) c(params[1], log(params[-1])) else log(params)
  c(delay, primary[1], log(primary[2]))
}

primary_case_info <- function(primary_id, dist_id, params, primary,
                              setting) {
  sprintf(
    "primary %d (%s), dist %d (%s), setting (%s)",
    primary_id, toString(primary), dist_id, toString(params),
    toString(setting)
  )
}

# The log CDF through primarycensored_lcdf(), which chooses between the
# series and the numerical path
scalar_models <- function() {
  compile_cpp_and_stan( # nolint: object_usage_linter.
    primary_probe_code(
      "  vector[N] d;",
      paste(
        "  for (i in 1:N) out[i] = primarycensored_lcdf(",
        "d[i] | dist_id, params, pwindow,",
        "dist_id == 18 ? negative_infinity() : 0, positive_infinity(),",
        "primary_id, pp);"
      )
    ),
    "pcdcppprimscalar"
  )
}

scalar_delays <- c(0.5, 1, 2.5, 4.2, 9)

# pwindow
scalar_settings <- list(1, 2)

scalar_data <- function(case, setting, w) {
  list(
    d = scalar_delays, N = length(scalar_delays), pwindow = setting, w = w
  )
}

# The vectorised series, called directly, so only where the series applies
vectorized_models <- function() {
  compile_cpp_and_stan( # nolint: object_usage_linter.
    primary_probe_code(
      "  int start;",
      paste(
        "  if (primary_id == 3) {",
        "    out = primarycensored_tlogis_lcdf_vectorized(",
        "      start, N, dist_id, params, pwindow, pp[1], pp[2]);",
        "  } else {",
        "    out = primarycensored_gumbel_lcdf_vectorized(",
        "      start, N, dist_id, params, pwindow, pp[1], pp[2]);",
        "  }",
        sep = "\n"
      )
    ),
    "pcdcppprimvec"
  )
}

# start, n, pwindow
vectorized_settings <- list(c(1, 15, 1), c(3, 15, 2))

vectorized_data <- function(case, setting, w) {
  list(start = setting[1], N = setting[2], pwindow = setting[3], w = w)
}

# The vectorised log PMF, which chooses between the series and the
# numerical path
pmf_models <- function() {
  compile_cpp_and_stan( # nolint: object_usage_linter.
    primary_probe_code(
      "  int max_delay;\n  real L;\n  real D;",
      paste(
        "  out = primarycensored_sone_lpmf_vectorized(",
        "max_delay, L, D, dist_id, params, pwindow, primary_id, pp);"
      )
    ),
    "pcdcppprimpmf"
  )
}

# max_delay, L, D, pwindow
pmf_settings <- list(c(12, 0, 13, 1), c(12, 0, Inf, 2))

pmf_data <- function(case, setting, w) {
  list(
    max_delay = setting[1], L = setting[2], D = setting[3],
    pwindow = setting[4], N = setting[1] + 1, w = w
  )
}

primary_probes <- list(
  scalar = list(
    models = scalar_models, settings = scalar_settings, data = scalar_data,
    n = function(setting) length(scalar_delays), pwindow = identity,
    series_only = FALSE, upper_tail = FALSE
  ),
  vectorized = list(
    models = vectorized_models, settings = vectorized_settings,
    data = vectorized_data, n = function(setting) setting[2],
    pwindow = function(setting) setting[3], series_only = TRUE,
    upper_tail = FALSE
  ),
  pmf = list(
    models = pmf_models, settings = pmf_settings, data = pmf_data,
    n = function(setting) setting[1] + 1,
    pwindow = function(setting) setting[4], series_only = FALSE,
    upper_tail = TRUE
  )
)

# Calls `check(case, primary, setting, data_base)` for every case and
# setting of `probe`, skipping cases where the series does not apply for
# probes that call it directly.
for_each_primary_case <- function(probe, check) {
  for (group in primary_cpp_cases) {
    for (delay in group$delays) {
      for (primary in group$primary) {
        for (setting in probe$settings) {
          applies <- primary_series_applies(
            group$primary_id, delay$dist_id, delay$params, primary,
            probe$pwindow(setting)
          )
          if (probe$series_only && !applies) next
          base <- list(
            dist_id = delay$dist_id, n_params = length(delay$params),
            primary_id = group$primary_id
          )
          check(group$primary_id, delay, primary, setting, base)
        }
      }
    }
  }
}

test_that("the cases cover the series and the numerical path", {
  for (group in primary_cpp_cases) {
    applies <- unlist(lapply(group$delays, function(delay) {
      vapply(group$primary, function(primary) {
        vapply(c(1, 2), function(pwindow) {
          primary_series_applies(
            group$primary_id, delay$dist_id, delay$params, primary, pwindow
          )
        }, logical(1))
      }, logical(2))
    }))
    expect_true(any(applies), info = group$primary_id)
    expect_false(all(applies), info = group$primary_id)
  }
})

for (probe_name in names(primary_probes)) {
  probe <- primary_probes[[probe_name]]

  test_that(paste("C++", probe_name, "primary values match the Stan body"), {
    models <- probe$models()
    for_each_primary_case(
      probe, function(primary_id, delay, primary, setting, base) {
        probe_data <- c(
          base, probe$data(delay, setting, rep(1, probe$n(setting)))
        )
        expect_cpp_values_match( # nolint: object_usage_linter.
          models, probe_data, delay$params,
          primary_case_info(
            primary_id, delay$dist_id, delay$params, primary, setting
          ),
          theta = primary_theta(delay$dist_id, delay$params, primary)
        )
      }
    )
  })

  test_that(paste("C++", probe_name, "primary gradient matches Stan"), {
    skip_if_no_log_prob() # nolint: object_usage_linter.
    models <- probe$models()
    set.seed(1)
    for_each_primary_case(
      probe, function(primary_id, delay, primary, setting, base) {
        probe_data <- c(
          base,
          probe$data(
            delay, setting, stats::runif(probe$n(setting), 0.5, 1.5)
          )
        )
        expect_cpp_gradient_matches( # nolint: object_usage_linter.
          models, probe_data, delay$params,
          primary_case_info(
            primary_id, delay$dist_id, delay$params, primary, setting
          ),
          upper_tail = probe$upper_tail,
          theta = primary_theta(delay$dist_id, delay$params, primary)
        )
      }
    )
  })
}
