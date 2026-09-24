skip_on_cran()

analytical_delays <- list(
  list(dist_id = 1, params = c(1.5, 0.5)), # Lognormal
  list(dist_id = 1, params = c(3.5, 0.05)), # Narrow lognormal
  list(dist_id = 2, params = c(2, 0.5)), # Gamma
  list(dist_id = 2, params = c(1.2, 0.02)), # Wide gamma
  list(dist_id = 3, params = c(1.5, 5)), # Weibull
  list(dist_id = 5, params = c(1.5, 3, 2)) # Generalised gamma
)

vectorized_lpmf <- function(delay, max_delay, L, D, pwindow,
                            primary_id = 1, primary_params = numeric(0)) {
  primarycensored_sone_lpmf_vectorized( # nolint: object_usage_linter.
    max_delay, L, D, delay$dist_id, delay$params, pwindow, primary_id,
    primary_params
  )
}

test_that(
  "primarycensored_sone_lpmf_vectorized matches primarycensored_lpmf for
   analytical delays",
  {
    settings <- list(
      list(max_delay = 20, L = 0, D = 21),
      list(max_delay = 40, L = 0, D = Inf),
      list(max_delay = 15, L = 2, D = 30),
      list(max_delay = 20, L = -Inf, D = 21),
      list(max_delay = 0, L = 0, D = 1)
    )
    for (delay in analytical_delays) {
      for (s in settings) {
        for (pwindow in c(1, 2, 3, 7, 1.5)) {
          vec <- vectorized_lpmf(delay, s$max_delay, s$L, s$D, pwindow)
          scalar <- vapply(0:s$max_delay, function(d) {
            primarycensored_lpmf(
              d, delay$dist_id, delay$params, pwindow, d + 1, s$L, s$D,
              1, numeric(0)
            )
          }, numeric(1))
          expect_equal(vec, scalar, tolerance = 1e-12)
        }
      }
    }
  }
)

test_that(
  "primarycensored_sone_pmf_vectorized sums to the mass in [L, max_delay + 1]",
  {
    for (delay in analytical_delays) {
      for (pwindow in c(1, 2)) {
        full <- exp(vectorized_lpmf(delay, 20, 0, 21, pwindow))
        expect_equal(sum(full), 1, tolerance = 1e-10)
        for (D in c(30, Inf)) {
          part <- exp(vectorized_lpmf(delay, 20, 2, D, pwindow))
          mass <- exp(primarycensored_lcdf(
            21, delay$dist_id, delay$params, pwindow, 2, D, 1, numeric(0)
          ))
          expect_equal(sum(part), mass, tolerance = 1e-10)
        }
      }
    }
  }
)

test_that(
  "primarycensored_sone_lpmf_vectorized is continuous in pwindow at an
   integer window",
  {
    for (delay in analytical_delays) {
      for (pwindow in c(1, 2)) {
        at <- exp(vectorized_lpmf(delay, 20, 0, 21, pwindow))
        near <- exp(vectorized_lpmf(delay, 20, 0, 21, pwindow + 1e-9))
        expect_equal(at, near, tolerance = 1e-6)
      }
    }
  }
)

test_that("primarycensored_sone_pmf_vectorized has finite gradients", {
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)))
  functions <- pcd_load_stan_functions(
    wrap_in_block = TRUE, write_to_file = FALSE
  )
  code <- paste0(
    functions, "\n",
    "data {\n",
    "  int max_delay;\n",
    "  int dist_id;\n",
    "  int n_params;\n",
    "  real pwindow;\n",
    "}\n",
    "parameters {\n",
    "  array[n_params] real params;\n",
    "}\n",
    "model {\n",
    "  array[0] real primary_params;\n",
    "  target += dot_product(\n",
    "    linspaced_vector(max_delay + 1, 1, 2),\n",
    "    primarycensored_sone_pmf_vectorized(\n",
    "      max_delay, 0, max_delay + 1, dist_id, params, pwindow, 1,\n",
    "      primary_params\n",
    "    )\n",
    "  );\n",
    "}\n"
  )
  path <- file.path(tempdir(), "pcd_sone_vectorized_gradient.stan")
  writeLines(code, path)
  model <- suppressMessages(suppressWarnings(cmdstanr::cmdstan_model(path)))
  for (delay in analytical_delays) {
    for (pwindow in c(1, 1.5)) {
      data_file <- tempfile(fileext = ".json")
      cmdstanr::write_stan_json(
        list(
          max_delay = 20, dist_id = delay$dist_id,
          n_params = length(delay$params), pwindow = pwindow
        ),
        data_file
      )
      params_file <- tempfile(fileext = ".json")
      cmdstanr::write_stan_json(
        list(params = as.array(delay$params)), params_file
      )
      out_file <- tempfile(fileext = ".csv")
      cmdstan_log <- suppressWarnings(system2(
        model$exe_file(),
        c(
          "log_prob", "jacobian=0",
          paste0("constrained_params=", params_file),
          "data", paste0("file=", data_file),
          "output", paste0("file=", out_file)
        ),
        stdout = TRUE, stderr = TRUE
      ))
      expect_null(
        attr(cmdstan_log, "status"),
        label = paste(cmdstan_log, collapse = "\n")
      )
      csv <- readLines(out_file)
      grad <- utils::read.csv(text = csv[!startsWith(csv, "#")])
      expect_true(
        all(is.finite(unlist(grad))),
        label = paste("dist", delay$dist_id, "pwindow", pwindow)
      )
    }
  }
})
