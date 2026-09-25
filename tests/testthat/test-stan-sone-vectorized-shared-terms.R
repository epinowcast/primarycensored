skip_on_cran()

# With a uniform primary and an integer pwindow,
# primarycensored_sone_lpmf_vectorized() uses
# primarycensored_analytical_lcdf_vectorized(), which computes the analytical
# uniform primary terms once per integer delay and shares them. These tests
# check it against per-delay Stan calls and the R functions.

vectorized_dists <- list(
  list(dist_id = 1L, params = list(
    c(1.5, 0.5), c(0.2, 1.2), c(-1, 0.3), c(3.5, 0.05), c(4, 2)
  )),
  list(dist_id = 2L, params = list(
    c(2, 0.5), c(0.5, 0.2), c(20, 4), c(1.2, 0.02), c(50, 0.5)
  )),
  list(dist_id = 3L, params = list(
    c(1.5, 5), c(0.7, 2), c(3, 30), c(5, 0.8)
  )),
  list(dist_id = 5L, params = list(c(1.5, 3, 2), c(0.8, 1, 0.5)))
)

per_delay_lcdf <- function(delays, dist_id, params, pwindow,
                           primary_id = 1L, primary_params = numeric(0)) {
  vapply(
    delays, primarycensored_lcdf, numeric(1), # nolint: object_usage_linter.
    dist_id, params, pwindow, 0, Inf, primary_id, primary_params
  )
}

test_that("check_for_analytical_vectorized needs uniform terms and an
  integer pwindow", {
  for (dist_id in c(1L, 2L, 3L, 5L)) {
    expect_identical(check_for_analytical_vectorized(dist_id, 1L, 1), 1L)
    expect_identical(check_for_analytical_vectorized(dist_id, 1L, 7), 1L)
    expect_identical(check_for_analytical_vectorized(dist_id, 1L, 1.5), 0L)
    expect_identical(check_for_analytical_vectorized(dist_id, 1L, 0.5), 0L)
    expect_identical(check_for_analytical_vectorized(dist_id, 2L, 1), 0L)
  }
  for (dist_id in c(4L, 18L, 26L, 27L, 28L)) {
    expect_identical(check_for_analytical_vectorized(dist_id, 1L, 1), 0L)
  }
})

test_that("uniform primary terms are -Inf for t <= 0", {
  for (dist in vectorized_dists) {
    for (t in c(0, -0.5, -3)) {
      expect_identical(
        primarycensored_uniform_terms(t, dist$dist_id, dist$params[[1]]),
        c(-Inf, -Inf)
      )
    }
  }
})

test_that(
  "primarycensored_uniform_lcdf_from_terms is -Inf when all terms underflow",
  {
    expect_identical(
      primarycensored_uniform_lcdf_from_terms(
        c(-Inf, -Inf), c(-Inf, -Inf), 2
      ),
      -Inf
    )
  }
)

test_that("primarycensored_analytical_lcdf_vectorized matches
  primarycensored_lcdf at each delay", {
  n <- 41L
  for (dist in vectorized_dists) {
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
            per_delay_lcdf(start:n, dist$dist_id, params, pwindow),
            info = info
          )
        }
      }
    }
  }
})

test_that("primarycensored_lcdf_vectorized matches primarycensored_lcdf
  off the analytical path", {
  cases <- list(
    list(
      dist_id = 2L, params = c(2, 0.5), pwindow = 1.5, primary_id = 1L,
      primary_params = numeric(0)
    ),
    list(
      dist_id = 1L, params = c(1.5, 0.5), pwindow = 1, primary_id = 2L,
      primary_params = 0.2
    ),
    list(
      dist_id = 4L, params = 0.3, pwindow = 1, primary_id = 1L,
      primary_params = numeric(0)
    )
  )
  for (case in cases) {
    vectorised <- primarycensored_lcdf_vectorized(
      1L, 10L, case$dist_id, case$params, case$pwindow, case$primary_id,
      case$primary_params
    )
    expect_identical(
      vectorised,
      per_delay_lcdf(
        1:10, case$dist_id, case$params, case$pwindow, case$primary_id,
        case$primary_params
      ),
      info = paste("dist", case$dist_id, "primary", case$primary_id)
    )
  }
})

test_that("primarycensored_sone_lpmf_vectorized matches primarycensored_lpmf
  under truncation", {
  dists <- list(
    list(dist_id = 1L, params = c(1.5, 0.5)),
    list(dist_id = 1L, params = c(0.2, 1.2)),
    list(dist_id = 2L, params = c(2, 0.5)),
    list(dist_id = 2L, params = c(20, 4)),
    list(dist_id = 3L, params = c(1.5, 5)),
    list(dist_id = 3L, params = c(0.7, 2)),
    list(dist_id = 5L, params = c(1.5, 3, 2))
  )
  settings <- list(
    list(max_delay = 0, L = 0, D = 1),
    list(max_delay = 1, L = 0, D = Inf),
    list(max_delay = 20, L = 0, D = 21),
    list(max_delay = 20, L = -Inf, D = 21),
    list(max_delay = 20, L = 3, D = 21),
    list(max_delay = 20, L = 2.5, D = 25.5),
    list(max_delay = 15, L = 0, D = 30),
    list(max_delay = 20, L = 4, D = Inf)
  )
  for (dist in dists) {
    for (s in settings) {
      for (pwindow in c(1, 2, 3, 7)) {
        info <- paste(
          "dist", dist$dist_id, "params", toString(dist$params),
          "max_delay", s$max_delay, "L", s$L, "D", s$D, "pwindow", pwindow
        )
        vectorised <- primarycensored_sone_lpmf_vectorized(
          s$max_delay, s$L, s$D, dist$dist_id, dist$params, pwindow,
          1L, numeric(0)
        )
        # A bin that L cuts is left out: primarycensored_lpmf() only takes
        # integer lower bounds
        delays <- 0:s$max_delay
        full <- delays >= s$L
        below <- delays + 1 <= s$L
        per_delay <- vapply(delays[full], function(d) {
          primarycensored_lpmf(
            d, dist$dist_id, dist$params, pwindow, d + 1, s$L, s$D,
            1L, numeric(0)
          )
        }, numeric(1))
        expect_equal(
          vectorised[full], per_delay,
          tolerance = 1e-10, info = info
        )
        expect_identical(vectorised[below], rep(-Inf, sum(below)))
        if (s$D == s$max_delay + 1) {
          expect_equal(
            sum(exp(vectorised)), 1,
            tolerance = 1e-10, info = info
          )
        }
      }
    }
  }
})

test_that("primarycensored_sone_pmf_vectorized matches R dprimarycensored
  for integer pwindow", {
  cases <- list(
    list(
      dist_id = 1L, params = c(1.5, 0.5), pdist = plnorm,
      args = list(meanlog = 1.5, sdlog = 0.5)
    ),
    list(
      dist_id = 2L, params = c(2, 0.5), pdist = pgamma,
      args = list(shape = 2, rate = 0.5)
    ),
    list(
      dist_id = 3L, params = c(1.5, 5), pdist = pweibull,
      args = list(shape = 1.5, scale = 5)
    )
  )
  if (requireNamespace("flexsurv", quietly = TRUE)) {
    cases[[length(cases) + 1]] <- list(
      dist_id = 5L, params = c(1.5, 3, 2), pdist = flexsurv::pgengamma.orig,
      args = list(shape = 1.5, scale = 3, k = 2)
    )
  }
  max_delay <- 15
  for (case in cases) {
    for (pwindow in c(1, 2, 3, 7)) {
      for (D in c(max_delay + 1, 30, Inf)) {
        info <- paste("dist", case$dist_id, "pwindow", pwindow, "D", D)
        stan_pmf <- primarycensored_sone_pmf_vectorized(
          max_delay, 0, D, case$dist_id, case$params, pwindow, 1L,
          numeric(0)
        )
        r_pmf <- do.call(dprimarycensored, c(
          list(
            0:max_delay, case$pdist,
            pwindow = pwindow, swindow = 1, D = D
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
vectorized_gradient_model <- function() {
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
    "  int max_delay;\n",
    "  real L;\n",
    "  real D;\n",
    "  int dist_id;\n",
    "  real pwindow;\n",
    "  vector[max_delay + 1] w;\n",
    "}\n",
    "parameters {\n",
    "  array[2] real<lower=0> params;\n",
    "}\n",
    "model {\n",
    "  target += dot_product(w, primarycensored_sone_pmf_vectorized(\n",
    "    max_delay, L, D, dist_id, params, pwindow, 1, {0.0}[1:0]\n",
    "  ));\n",
    "}\n"
  )
  path <- file.path(tempdir(), "pcd_sone_vectorized_gradient.stan")
  writeLines(code, path)
  suppressMessages(suppressWarnings(cmdstanr::cmdstan_model(path)))
}

vectorized_gradient_at <- function(model, stan_data, params) {
  stan_gradient_at( # nolint: object_usage_linter.
    model,
    data = stan_data, init = list(params = params)
  )
}

test_that("primarycensored_sone_pmf_vectorized gradients match finite
  differences with shared terms", {
  model <- vectorized_gradient_model()
  cases <- list(
    list(dist_id = 1L, params = c(1.5, 0.5)),
    list(dist_id = 2L, params = c(2, 0.5)),
    list(dist_id = 3L, params = c(1.5, 5))
  )
  max_delay <- 15
  for (case in cases) {
    for (pwindow in c(1, 3)) {
      for (s in list(list(L = 0, D = max_delay + 1), list(L = 2, D = Inf))) {
        info <- paste(
          "dist", case$dist_id, "pwindow", pwindow, "L", s$L, "D", s$D
        )
        res <- vectorized_gradient_at(model, list(
          max_delay = max_delay, L = s$L, D = s$D, dist_id = case$dist_id,
          pwindow = pwindow, w = seq(0.5, 1.5, length.out = max_delay + 1)
        ), case$params)
        expect_length(res$gradient, 2)
        expect_true(all(is.finite(res$gradient)), info = info)
        expect_equal(
          res$gradient, res$finite_diff,
          tolerance = 1e-4, info = info
        )
      }
    }
  }
})
