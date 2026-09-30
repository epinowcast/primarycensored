skip_on_cran()

# The truncated logistic primary is primary_id 3 with
# primary_params = c(location, scale). The Stan functions should agree with
# the R functions dtlogis(), ptlogis() and rtlogis().
stan_tlogis_cases <- list(
  list(min = 0, max = 1, location = -0.5, scale = 0.2, label = "m below"),
  list(min = 0, max = 1, location = 0, scale = 0.1, label = "m at min"),
  list(min = 0, max = 1, location = 0.5, scale = 0.2, label = "m inside"),
  list(min = 0, max = 2, location = 1, scale = 1, label = "m mid, wide"),
  list(min = 0, max = 1, location = 1.5, scale = 0.1, label = "m above"),
  list(min = 0, max = 2, location = 40, scale = 5, label = "flat, far above"),
  list(min = 0, max = 1, location = -6, scale = 0.02, label = "far below"),
  list(min = 10, max = 20, location = 14, scale = 2, label = "non-zero min")
)

for (tc in stan_tlogis_cases) {
  x <- seq(tc$min, tc$max, length.out = 11)
  inner <- x[x > tc$min & x < tc$max]

  test_that(paste("Stan tlogis_lpdf matches R dtlogis:", tc$label), {
    stan_lpdf <- vapply(
      x, tlogis_lpdf, numeric(1), tc$min, tc$max, tc$location, tc$scale
    )
    expect_equal(
      stan_lpdf,
      dtlogis(x, tc$min, tc$max, tc$location, tc$scale, log = TRUE),
      tolerance = 1e-9
    )
    expect_identical(
      tlogis_lpdf(tc$min - 0.1, tc$min, tc$max, tc$location, tc$scale), -Inf
    )
    expect_identical(
      tlogis_lpdf(tc$max + 0.1, tc$min, tc$max, tc$location, tc$scale), -Inf
    )
  })

  test_that(paste("Stan tlogis_lcdf and tlogis_cdf match ptlogis:", tc$label), {
    stan_lcdf <- vapply(
      inner, tlogis_lcdf, numeric(1), tc$min, tc$max, tc$location, tc$scale
    )
    expect_equal(
      stan_lcdf,
      ptlogis(inner, tc$min, tc$max, tc$location, tc$scale, log.p = TRUE),
      tolerance = 1e-9
    )
    stan_cdf <- vapply(
      x, tlogis_cdf, numeric(1), tc$min, tc$max, tc$location, tc$scale
    )
    expect_equal(
      stan_cdf,
      ptlogis(x, tc$min, tc$max, tc$location, tc$scale),
      tolerance = 1e-9
    )
    expect_identical(
      tlogis_lcdf(tc$min - 0.1, tc$min, tc$max, tc$location, tc$scale), -Inf
    )
    expect_identical(
      tlogis_lcdf(tc$max + 0.1, tc$min, tc$max, tc$location, tc$scale), 0
    )
  })

  test_that(paste("Stan tlogis_rng matches the R CDF:", tc$label), {
    n <- 10000
    samples <- replicate(
      n, tlogis_rng(tc$min, tc$max, tc$location, tc$scale)
    )
    expect_true(all(samples >= tc$min & samples <= tc$max))
    ks <- suppressWarnings(ks.test(
      samples,
      function(q) ptlogis(q, tc$min, tc$max, tc$location, tc$scale)
    ))
    expect_gt(ks$p.value, 1e-3)
  })

  test_that(paste("primary_lpdf and primary_lcdf use tlogis:", tc$label), {
    params <- c(tc$location, tc$scale)
    pwindow <- tc$max - tc$min
    expect_equal(
      vapply(
        x, primary_lpdf, numeric(1), 3L, params, tc$min, tc$max
      ),
      dtlogis(x, tc$min, tc$max, tc$location, tc$scale, log = TRUE),
      tolerance = 1e-9
    )
    # primary_lcdf is on the window [0, pwindow], with the location on the
    # same scale
    p <- x - tc$min
    p <- p[p > 0 & p < pwindow]
    expect_equal(
      vapply(p, primary_lcdf, numeric(1), 3L, params, pwindow),
      ptlogis(p, 0, pwindow, tc$location, tc$scale, log.p = TRUE),
      tolerance = 1e-9
    )
  })
}

test_that("primary_lcdf is on [0, pwindow] at and beyond the ends", {
  params <- c(0.5, 0.2)
  expect_identical(primary_lcdf(0, 3L, params, 2), -Inf)
  expect_identical(primary_lcdf(-1, 3L, params, 2), -Inf)
  expect_identical(primary_lcdf(2, 3L, params, 2), 0)
  expect_identical(primary_lcdf(3, 3L, params, 2), 0)
})

test_that("tlogis_lcdf and tlogis_lpdf have finite gradients in the location
  and scale", {
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)))
  functions <- pcd_load_stan_functions(
    wrap_in_block = TRUE, write_to_file = FALSE
  )
  code <- paste0(
    functions, "\n",
    "data {\n",
    "  real x;\n",
    "  real xmin;\n",
    "  real xmax;\n",
    "  int use_pdf;\n",
    "}\n",
    "parameters {\n",
    "  real location;\n",
    "  real<lower=0> scale;\n",
    "}\n",
    "model {\n",
    "  if (use_pdf) {\n",
    "    target += tlogis_lpdf(x | xmin, xmax, location, scale);\n",
    "  } else {\n",
    "    target += tlogis_lcdf(x | xmin, xmax, location, scale);\n",
    "  }\n",
    "}\n"
  )
  path <- file.path(tempdir(), "pcd_tlogis_gradient.stan")
  writeLines(code, path)
  model <- suppressMessages(suppressWarnings(cmdstanr::cmdstan_model(path)))
  points <- list(
    list(x = 0.3, xmin = 0, xmax = 1, location = 0.5, scale = 0.2),
    list(x = 0.8, xmin = 0, xmax = 1, location = -0.5, scale = 0.2),
    list(x = 0.2, xmin = 0, xmax = 1, location = 1.5, scale = 0.1),
    list(x = 1.5, xmin = 0, xmax = 2, location = 40, scale = 5),
    list(x = 0.5, xmin = 0, xmax = 1, location = 7, scale = 0.05),
    list(x = 0.5, xmin = 0, xmax = 1, location = -7, scale = 0.05)
  )
  for (point in points) {
    for (use_pdf in c(0L, 1L)) {
      res <- stan_gradient_at(
        model,
        data = c(point[c("x", "xmin", "xmax")], list(use_pdf = use_pdf)),
        init = point[c("location", "scale")]
      )
      label <- paste0(
        "x ", point$x, ", location ", point$location, ", scale ", point$scale,
        ", pdf ", use_pdf
      )
      expect_false(res$gradient_not_finite, info = label)
      expect_false(res$rejected, info = label)
      expect_length(res$gradient, 2)
      expect_true(all(is.finite(res$gradient)), info = label)
      allowed <- 1e-5 * pmax(abs(res$finite_diff), 1e-2)
      expect_true(
        all(abs(res$gradient - res$finite_diff) <= allowed),
        info = paste0(
          label, ": gradient ", toString(signif(res$gradient, 5)),
          ", finite difference ", toString(signif(res$finite_diff, 5))
        )
      )
    }
  }
})

test_that("the ODE path integrates against the truncated logistic primary", {
  # primarycensored_numeric_cdf is the ODE CDF the analytical solutions are
  # tested against. It has relative and absolute tolerances of 1e-6.
  cases <- list(
    list(dist_id = 2L, params = c(2.5, 0.4), pdist = pgamma,
         args = list(shape = 2.5, rate = 0.4), pos = TRUE),
    list(dist_id = 4L, params = 0.3, pdist = pexp,
         args = list(rate = 0.3), pos = TRUE),
    list(dist_id = 18L, params = c(3, 2), pdist = pnorm,
         args = list(mean = 3, sd = 2), pos = FALSE)
  )
  for (case in cases) {
    for (window in list(c(1, 0.5, 0.2), c(2, 1.5, 0.3), c(1, -0.5, 1))) {
      pwindow <- window[1]
      primary_params <- window[2:3]
      obj <- do.call(
        new_pcens,
        c(
          list(
            pdist = case$pdist, dprimary = dtlogis,
            primary_args = list(
              location = primary_params[1], scale = primary_params[2]
            )
          ),
          case$args
        )
      )
      d <- c(0.6, 1.5, 3, 6, 12)
      r_numeric <- pcens_cdf(obj, d, pwindow, use_numeric = TRUE)
      stan_ode <- vapply(
        d, primarycensored_numeric_cdf, numeric(1),
        case$dist_id, case$params, pwindow, 3L, primary_params
      )
      expect_lt(
        max(abs(r_numeric - stan_ode)), 1e-5,
        label = paste(case$dist_id, toString(window))
      )
    }
  }
})
