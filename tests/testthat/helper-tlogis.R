# Reference CDF by integrating the delay CDF against the window density to a
# relative 1e-12, split at the delay kink and around the location.
tlogis_reference <- function(q, pwindow, location, scale, cdf,
                             positive = TRUE) {
  vapply(q, function(qq) {
    integrand <- function(p) {
      cdf(qq - p) * dtlogis(p, 0, pwindow, location, scale)
    }
    centre <- min(max(location, 0), pwindow)
    spike <- centre + scale * c(-40, -20, -10, -5, -2, 0, 2, 5, 10, 20, 40)
    breaks <- sort(unique(c(
      0, pwindow,
      if (positive && qq > 0 && qq < pwindow) qq,
      spike[spike > 0 & spike < pwindow]
    )))
    sum(vapply(seq_len(length(breaks) - 1L), function(i) {
      stats::integrate(
        integrand, breaks[i], breaks[i + 1L],
        rel.tol = 1e-12, abs.tol = 0, subdivisions = 2000L
      )$value
    }, numeric(1)))
  }, numeric(1))
}

# Delay CDF of a family from `exptilt_families()` as a function of the
# quantile alone.
tlogis_object <- function(family, location, scale) {
  do.call(
    new_pcens,
    c(
      list(
        pdist = family$pdist, dprimary = dtlogis,
        primary_args = list(location = location, scale = scale)
      ),
      family$args
    )
  )
}

# Label for a failed expectation.
tlogis_label <- function(family, pwindow, location, scale) {
  sprintf(
    "%s, pwindow = %g, location = %g, scale = %g",
    family$label, pwindow, location, scale
  )
}

tlogis_stan_cases <- list(
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
    dist_id = 2L, params = c(2.5, 1000), pdist = pgamma,
    args = list(shape = 2.5, rate = 1000)
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

tlogis_case_label <- function(case, ...) {
  paste0(
    "dist ", case$dist_id, " params ", toString(case$params), ", ",
    paste(names(list(...)), unlist(list(...)), sep = " = ", collapse = ", ")
  )
}

tlogis_case_cdf <- function(case) {
  function(x) do.call(case$pdist, c(list(x), case$args))
}

# Internal lower bound used by primarycensored_lcdf for each support
tlogis_case_lower <- function(case) {
  if (case$dist_id == 18L) -Inf else 0
}

tlogis_case_obj <- function(case, location, scale) {
  tlogis_object(list(pdist = case$pdist, args = case$args), location, scale)
}

# Gradients are only observable from a compiled model, so this builds a
# minimal one whose target is the log CDF or the vectorised log PMF, and runs
# `stan_gradient_at()` from helper-stan-gradient.R.
tlogis_gradient_model <- function() {
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
    "  real location;\n",
    "  real<lower=0> scale;\n",
    "}\n",
    "model {\n",
    "  array[2] real all_params = {p1, p2};\n",
    "  array[n_params] real params = all_params[1:n_params];\n",
    "  if (vectorised) {\n",
    "    target += sum(primarycensored_sone_lpmf_vectorized(\n",
    "      to_int(d), L, positive_infinity(), dist_id, params, pwindow, 3,\n",
    "      {location, scale}\n",
    "    ));\n",
    "  } else {\n",
    "    target += primarycensored_lcdf(\n",
    "      d | dist_id, params, pwindow, L, positive_infinity(), 3,\n",
    "      {location, scale}\n",
    "    );\n",
    "  }\n",
    "}\n"
  )
  path <- file.path(tempdir(), "pcd_tlogis_cdf_gradient.stan")
  writeLines(code, path)
  suppressMessages(suppressWarnings(cmdstanr::cmdstan_model(path)))
}

tlogis_gradient_at <- function(model, case, d, pwindow, location, scale,
                               vectorised = FALSE) {
  init <- list(
    p1 = case$params[1],
    p2 = if (length(case$params) > 1) case$params[2] else 1,
    location = location, scale = scale
  )
  stan_gradient_at( # nolint: object_usage_linter.
    model,
    data = list(
      dist_id = case$dist_id, n_params = length(case$params),
      vectorised = as.integer(vectorised), d = d, pwindow = pwindow,
      L = tlogis_case_lower(case)
    ),
    init = init
  )
}
