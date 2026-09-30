skip_on_cran()

# The uniform primary analytical solutions in Stan are built from terms at d
# and at q = max(d - pwindow, 0), see primarycensored_uniform_terms(). These
# tests check the log CDF against tight numerical integration and the
# gradients against finite differences, across parameters, windows, tails,
# `d <= pwindow` and `d` near 0. The reference and cases are in
# helper-uniform-reference.R.

stan_uniform_lcdf <- function(d, case, args, pwindow) {
  lcdf <- primarycensored_analytical_lcdf # nolint: object_usage_linter.
  vapply(
    d, lcdf, numeric(1),
    case$stan_id, case$stan_params(args), pwindow, 0, Inf, 1L, numeric(0)
  )
}

test_that("Stan uniform primary analytical log CDFs match tight numerical
  integration across parameters, windows and tails", {
  d <- c(1e-8, 1e-3, 0.05, 0.3, 0.75, 1, 1.5, 2.5, 6, 12, 40)
  for (case in unif_cases()) {
    for (args in case$grid) {
      for (pwindow in c(0.5, 1, 3)) {
        info <- paste(
          case$name, toString(unlist(args)), "pwindow", pwindow
        )
        reference <- do.call(
          unif_reference, c(list(case$pdist, d, pwindow), args)
        )
        result <- exp(stan_uniform_lcdf(d, case, args, pwindow))
        expect_close(result, reference, info = info)
      }
    }
  }
})

test_that("Stan uniform primary analytical log CDFs are -Inf at and below 0
  and 0 far into the upper tail", {
  for (case in unif_cases()) {
    args <- case$grid[[1]]
    expect_identical(
      stan_uniform_lcdf(c(-1, 0), case, args, 1), c(-Inf, -Inf)
    )
    expect_equal(stan_uniform_lcdf(1e4, case, args, 1), 0, tolerance = 1e-10)
  }
})

gradient_model <- function() {
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
    "  int dist_id;\n",
    "  int n_params;\n",
    "  real pwindow;\n",
    "  array[n] real d;\n",
    "}\n",
    "parameters {\n",
    "  array[n_params] real params;\n",
    "}\n",
    "model {\n",
    "  for (i in 1:n) {\n",
    "    target += primarycensored_analytical_lcdf(\n",
    "      d[i] | dist_id, params, pwindow, 0, positive_infinity(), 1,\n",
    "      {0.0}[1:0]\n",
    "    );\n",
    "  }\n",
    "}\n"
  )
  path <- file.path(tempdir(), "pcd_uniform_terms_gradient.stan")
  writeLines(code, path)
  suppressMessages(suppressWarnings(cmdstanr::cmdstan_model(path)))
}

test_that("Stan uniform primary analytical log CDFs have gradients that match
  finite differences", {
  model <- gradient_model()
  # Delays at the lower tail, central mass and upper tail of each case, and
  # for d <= pwindow, where q = 0 and only the terms at d contribute.
  # Weibull and generalised gamma shapes are 0.8 and above. The autodiff
  # gradient of gamma_p in a = 1 + 1 / shape is inaccurate for shapes below
  # about 0.2, see #393.
  cases <- list(
    list(id = 2L, params = c(3, 0.5), d = c(0.4, 1.5, 4, 9, 20)),
    list(id = 2L, params = c(0.6, 0.2), d = c(0.3, 2, 8, 40)),
    list(id = 1L, params = c(1.5, 0.6), d = c(0.9, 2, 4.5, 9, 20)),
    list(id = 1L, params = c(0.3, 1.2), d = c(0.2, 1, 3, 30)),
    list(id = 3L, params = c(1.6, 6), d = c(0.4, 2, 5, 9, 20)),
    list(id = 3L, params = c(0.8, 2), d = c(0.3, 1.5, 6, 25)),
    list(id = 5L, params = c(1.5, 4, 1.2), d = c(0.4, 2, 5, 9, 20)),
    list(id = 5L, params = c(0.8, 2, 0.5), d = c(0.3, 1.5, 6, 25))
  )
  for (case in cases) {
    for (pwindow in c(1, 2.5)) {
      info <- paste(
        "dist", case$id, toString(case$params), "pwindow", pwindow
      )
      res <- stan_gradient_at( # nolint: object_usage_linter.
        model,
        data = list(
          n = length(case$d), dist_id = case$id,
          n_params = length(case$params), pwindow = pwindow,
          d = case$d
        ),
        init = list(params = case$params)
      )
      expect_false(res$rejected, info = info)
      expect_false(res$gradient_not_finite, info = info)
      expect_length(res$gradient, length(case$params))
      expect_true(all(is.finite(res$gradient)), info = info)
      expect_equal(
        res$gradient, res$finite_diff,
        tolerance = 1e-5, info = info
      )
    }
  }
})
