# Helpers for the log-logistic delay tests.
#
# `pllogis()` has the parameters of `flexsurv::pllogis()`, a shape and a
# scale, so the tests do not need flexsurv. It is named as the function
# exported by flexsurv and actuar is, which is how the package finds the
# analytic methods for it.
pllogis_test <- add_name_attribute(
  function(q, shape, scale = 1) {
    stats::plogis(shape * (log(pmax(q, 0)) - log(scale)))
  },
  "pllogis"
)

# Delay families covering a heavy tail, exponential-like and sharp shapes
loglogistic_families <- function() {
  list(
    list(
      label = "log-logistic shape 0.6", pdist = pllogis_test,
      args = list(shape = 0.6, scale = 2), positive = TRUE
    ),
    list(
      label = "log-logistic shape 1", pdist = pllogis_test,
      args = list(shape = 1, scale = 3), positive = TRUE
    ),
    list(
      label = "log-logistic shape 2", pdist = pllogis_test,
      args = list(shape = 2, scale = 5), positive = TRUE
    ),
    list(
      label = "log-logistic shape 4.5", pdist = pllogis_test,
      args = list(shape = 4.5, scale = 8), positive = TRUE
    ),
    list(
      label = "log-logistic shape 12", pdist = pllogis_test,
      args = list(shape = 12, scale = 6), positive = TRUE
    )
  )
}

# Reference for the partial moment M_a(A) scaled as in the series
# r_a(A) = M_a(A) (1 + A) / A^(a + 1), by integrating over s = y / A so the
# integrand does not overflow.
loglogistic_ratio_reference <- function(a, A) {
  integrand <- function(s) s^a / (1 + A * s)^2
  # The integrand peaks near 1 / A for small a, so split there
  breaks <- sort(unique(c(0, min(1, 1 / A), 1)))
  inner <- sum(vapply(seq_len(length(breaks) - 1L), function(i) {
    stats::integrate(
      integrand, breaks[i], breaks[i + 1L],
      rel.tol = 1e-13, abs.tol = 0, subdivisions = 5000L
    )$value
  }, numeric(1)))
  (1 + A) * inner
}

# A log-logistic delay with a uniform primary
uniform_object <- function(family) {
  do.call(
    new_pcens,
    c(
      list(
        pdist = family$pdist, dprimary = dunif, primary_args = list()
      ),
      family$args
    )
  )
}

# Independent reference for the primary event censored CDF and survival of a
# log-logistic delay. The CDF is the expectation of F(q - P) and the survival
# the expectation of S(q - P) over the primary event time P in (0, pwindow)
# with a uniform density for `rho = 0` and otherwise density
# rho exp(rho p) / (exp(rho pwindow) - 1). Each is integrated over p with
# breaks where the delay CDF moves, so the upper tail is accurate to 1e-13 of
# the survival rather than of the CDF.
loglogistic_censored_reference <- function(q, shape, scale, rho, pwindow) {
  primary_density <- function(p) {
    if (rho == 0) {
      rep(1 / pwindow, length(p))
    } else {
      rho * exp(rho * p) / expm1(rho * pwindow)
    }
  }
  lcdf <- function(u) {
    ifelse(u <= 0, -Inf, -log1p(exp(-shape * (log(pmax(u, 1e-300)) -
      log(scale)))))
  }
  breaks <- sort(unique(c(
    0, pwindow, if (q > 0 && q < pwindow) q,
    q - scale * exp(seq(-40, 40, by = 0.5) / shape)
  )))
  breaks <- breaks[breaks >= 0 & breaks <= pwindow]
  integral <- function(g) {
    sum(vapply(seq_len(length(breaks) - 1L), function(i) {
      stats::integrate(
        function(p) g(p) * primary_density(p), breaks[i], breaks[i + 1L],
        rel.tol = 1e-13, abs.tol = 0, subdivisions = 10000L,
        stop.on.error = FALSE
      )$value
    }, numeric(1)))
  }
  cdf <- integral(function(p) exp(lcdf(q - p)))
  survival <- integral(function(p) {
    u <- q - p
    ifelse(u <= 0, 1, -expm1(lcdf(u)))
  })
  c(cdf = cdf, survival = survival)
}

# Cases where the difference G_1(q) - G_1(q - w) of the uniform solution
# loses precision, each (shape, scale, pwindow, q)
uniform_conditioning_cases <- list(
  c(0.01, 1, 1, 1e17),
  c(0.01, 1, 1, 1e8),
  c(0.02, 1, 1, 1e10),
  c(0.3, 5, 1, 1e5),
  c(0.5, 5, 1, 1e4),
  c(0.5, 5, 1, 1e5),
  c(0.5, 5, 1, 1e6),
  c(1, 5, 1, 1e3),
  c(2, 5, 1, 1e3),
  c(115.5, 16.9, 0.00179, 19.8)
)

# Cases of the uniform PMF at a large delay, where the PMF is small next to
# the smaller tail in a heavy tail, each (shape, scale, pwindow, d)
uniform_pmf_cases <- list(
  c(0.05, 1, 1, 1e6),
  c(0.3, 5, 1, 1e5),
  c(0.5, 5, 1, 1e4),
  c(0.5, 5, 1, 1e5),
  c(0.2, 7.94, 2, 2e5),
  c(1, 5, 1, 1e3),
  c(2, 5, 1, 500)
)

# Inverse CDF sampler for the log-logistic
rllogis_test <- function(n, shape, scale = 1) {
  u <- stats::runif(n)
  scale * (u / (1 - u))^(1 / shape)
}

# Samples from rprimarycensored() with a uniform (`rho = 0`) or a tilted
# primary event window
loglogistic_samples <- function(n, family, pwindow, rho) {
  do.call(
    rprimarycensored,
    c(
      list(
        n = n, rdist = rllogis_test, pwindow = pwindow, swindow = 0,
        rprimary = if (rho == 0) stats::runif else rexpgrowth,
        rprimary_args = if (rho == 0) list() else list(r = rho)
      ),
      family$args
    )
  )
}

# Parameters of the Stan tests, `c(scale, shape)`, with the R arguments
ll_cases <- list(
  list(params = c(2, 0.6), args = list(shape = 0.6, scale = 2)),
  list(params = c(3, 1), args = list(shape = 1, scale = 3)),
  list(params = c(5, 2), args = list(shape = 2, scale = 5)),
  list(params = c(8, 4.5), args = list(shape = 4.5, scale = 8)),
  list(params = c(6, 12), args = list(shape = 12, scale = 6))
)

ll_case_label <- function(case, ...) {
  paste0(
    "params ", toString(case$params), ", ",
    paste(names(list(...)), unlist(list(...)), sep = " = ", collapse = ", ")
  )
}

ll_case_cdf <- function(case) {
  function(x) do.call(pllogis_test, c(list(x), case$args))
}

ll_case_family <- function(case) {
  list(pdist = pllogis_test, args = case$args)
}

# The ODE branch of primarycensored_cdf() for given delays, as the CDF of a
# fixed parameter run of a model that calls the ODE solver directly
ll_ode_model <- function() {
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
    "  int n;\n",
    "  array[n] real d;\n",
    "  array[2] real params;\n",
    "  real pwindow;\n",
    "  int primary_id;\n",
    "  int n_primary;\n",
    "  array[n_primary] real primary_params;\n",
    "}\n",
    "generated quantities {\n",
    "  array[n] real cdf;\n",
    "  for (i in 1:n) {\n",
    "    real lower_bound = fmax(d[i] - pwindow, 0);\n",
    "    array[2 + n_primary] real theta =\n",
    "      append_array(params, primary_params);\n",
    "    array[4] int ids = {31, primary_id, 2, n_primary};\n",
    "    cdf[i] = ode_rk45(\n",
    "      primarycensored_ode, rep_vector(0.0, 1), lower_bound, {d[i]},\n",
    "      theta, {d[i], pwindow}, ids\n",
    "    )[1, 1];\n",
    "  }\n",
    "}\n"
  )
  path <- file.path(tempdir(), "pcd_loglogistic_ode.stan")
  writeLines(code, path)
  suppressMessages(suppressWarnings(cmdstanr::cmdstan_model(path)))
}

ll_ode_cdf <- function(model, params, d, pwindow, primary_id, primary_params) {
  fit <- model$sample(
    data = list(
      n = length(d), d = as.array(d), params = as.array(params),
      pwindow = pwindow, primary_id = primary_id,
      n_primary = length(primary_params),
      primary_params = as.array(primary_params)
    ),
    fixed_param = TRUE, chains = 1, iter_sampling = 1, refresh = 0,
    show_messages = FALSE, sig_figs = 18
  )
  as.numeric(fit$draws("cdf", format = "matrix"))
}

# Gradients are only observable from a compiled model, so this builds a
# minimal one whose target is the log CDF or the vectorised log PMF, and runs
# `stan_gradient_at()` from helper-stan-gradient.R.
ll_gradient_model <- function() {
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
    "  int primary_id;\n",
    "  int vectorised;\n",
    "  real d;\n",
    "  real pwindow;\n",
    "}\n",
    "parameters {\n",
    "  real<lower=0> scale;\n",
    "  real<lower=0> shape;\n",
    "  real rho;\n",
    "}\n",
    "model {\n",
    "  array[2] real params = {scale, shape};\n",
    "  array[primary_id == 2 ? 1 : 0] real primary_params;\n",
    "  if (primary_id == 2) primary_params[1] = rho;\n",
    "  if (vectorised) {\n",
    "    target += sum(primarycensored_sone_lpmf_vectorized(\n",
    "      to_int(d), 0.0, positive_infinity(), 31, params, pwindow,\n",
    "      primary_id, primary_params\n",
    "    ));\n",
    "  } else {\n",
    "    target += primarycensored_lcdf(\n",
    "      d | 31, params, pwindow, 0.0, positive_infinity(), primary_id,\n",
    "      primary_params\n",
    "    );\n",
    "  }\n",
    "}\n"
  )
  path <- file.path(tempdir(), "pcd_loglogistic_gradient.stan")
  writeLines(code, path)
  suppressMessages(suppressWarnings(cmdstanr::cmdstan_model(path)))
}

# Checks the gradient is finite and matches the finite difference gradient one
# component at a time, relative to its size with a floor for tiny ones.
expect_ll_gradient <- function(model, case, point, primary_id,
                               vectorised = FALSE, tolerance = 1e-4) {
  rho <- if (is.null(point$rho)) 0 else point$rho
  label <- ll_case_label(
    case, d = point$d, pwindow = point$pwindow, r = rho,
    vectorised = vectorised
  )
  res <- stan_gradient_at( # nolint: object_usage_linter.
    model,
    data = list(
      primary_id = primary_id, vectorised = as.integer(vectorised),
      d = point$d, pwindow = point$pwindow
    ),
    init = list(scale = case$params[1], shape = case$params[2], rho = rho)
  )
  expect_false(res$gradient_not_finite, info = label)
  expect_false(res$rejected, info = label)
  expect_true(all(is.finite(res$gradient)), info = label)
  allowed <- tolerance * pmax(abs(res$finite_diff), 1e-2)
  expect_true(
    all(abs(res$gradient - res$finite_diff) <= allowed),
    info = paste0(
      label, ": gradient ", toString(signif(res$gradient, 5)),
      ", finite difference ", toString(signif(res$finite_diff, 5))
    )
  )
  invisible(res)
}
