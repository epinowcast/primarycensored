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

# Compared in absolute terms as the ODE error is absolute
tlogis_narrow_cases <- list(
  list(dist_id = 4L, params = 1.5, pdist = pexp, args = list(rate = 1.5)),
  list(
    dist_id = 2L, params = c(2, 1), pdist = pgamma,
    args = list(shape = 2, rate = 1)
  ),
  list(
    dist_id = 1L, params = c(1, 0.5), pdist = plnorm,
    args = list(meanlog = 1, sdlog = 0.5)
  ),
  list(
    dist_id = 1L, params = c(0, 0.5), pdist = plnorm,
    args = list(meanlog = 0, sdlog = 0.5)
  ),
  list(
    dist_id = 3L, params = c(2, 2), pdist = pweibull,
    args = list(shape = 2, scale = 2)
  ),
  list(
    dist_id = 18L, params = c(3, 2), pdist = pnorm,
    args = list(mean = 3, sd = 2)
  )
)

test_that("the ODE path resolves a narrow truncated logistic primary", {
  pwindow <- 2
  d <- c(0.5, 1, 3, 10)
  for (case in tlogis_narrow_cases) {
    cdf <- tlogis_case_cdf(case)
    for (location in c(-0.5, 0.7, 2.5)) {
      for (scale in c(0.02, 0.005, 0.001)) {
        window <- c(location, scale)
        expected <- tlogis_reference(
          d, pwindow, location, scale, cdf, case$dist_id != 18L
        )
        ode <- vapply(
          d, primarycensored_numeric_cdf, numeric(1),
          case$dist_id, case$params, pwindow, 3L, window
        )
        expect_lt(
          max(abs(ode - expected)), 1e-7,
          label = tlogis_case_label(
            case,
            location = location, scale = scale
          )
        )
      }
    }
  }
})

test_that("the log CDF and log PMF are correct for a narrow primary", {
  pwindow <- 2
  for (case in tlogis_narrow_cases) {
    cdf <- tlogis_case_cdf(case)
    lower <- tlogis_case_lower(case)
    for (location in c(0.7, 1)) {
      for (scale in c(0.02, 0.005)) {
        window <- c(location, scale)
        label <- tlogis_case_label(case, location = location, scale = scale)
        d <- c(0.5, 1, 3, 10)
        expected <- tlogis_reference(
          d, pwindow, location, scale, cdf, case$dist_id != 18L
        )
        lcdf <- vapply(
          d, primarycensored_lcdf, numeric(1),
          case$dist_id, case$params, pwindow, lower, Inf, 3L, window
        )
        expect_lt(max(abs(exp(lcdf) - expected)), 1e-7, label = label)
        # The PMF over integer delays
        ref <- tlogis_reference(
          0:6, pwindow, location, scale, cdf, case$dist_id != 18L
        )
        pmf <- exp(primarycensored_sone_lpmf_vectorized(
          5, lower, Inf, case$dist_id, case$params, pwindow, 3L, window
        ))
        expect_lt(
          max(abs(pmf - diff(ref)[1:6])), 1e-7,
          label = paste("pmf:", label)
        )
      }
    }
  }
})

test_that("non-parametric delays are analytic with a truncated logistic
  primary", {
  boundaries <- c(0, 1, 3, 6, 10)
  pmf <- c(0.2, 0.3, 0.35, 0.15)
  params <- c(boundaries, pmf)
  window <- c(0.5, 0.3)
  obj <- new_pcens(
    pdist = pdiscretestep, dprimary = dtlogis,
    primary_args = list(location = 0.5, scale = 0.3),
    boundaries = boundaries, pmf = pmf
  )
  d <- c(0.5, 1.5, 3, 5, 8, 12)
  expect_identical(check_for_analytical(26L, 3L), 1L)
  analytic <- vapply(
    d, primarycensored_cdf, numeric(1), 26L, params, 2, -Inf, Inf, 3L, window
  )
  expect_equal(analytic, pcens_cdf(obj, d, 2), tolerance = 1e-9)
  # The step CDF has kinks, so the ODE path is accurate to about 1e-4
  ode <- vapply(
    d, primarycensored_numeric_cdf, numeric(1), 26L, params, 2, 3L, window
  )
  expect_lt(max(abs(analytic - ode)), 1e-4)
})

test_that("the ODE path matches Stan random draws", {
  set.seed(20260930)
  n <- 5000
  pwindow <- 2
  probs <- seq(0.1, 0.9, by = 0.2)
  for (case in tlogis_stan_cases[c(1, 4, 7)]) {
    rdist <- switch(as.character(case$dist_id),
      "4" = function(n) rexp(n, case$args$rate),
      "2" = function(n) rgamma(n, case$args$shape, case$args$rate),
      "18" = function(n) rnorm(n, case$args$mean, case$args$sd)
    )
    for (window in list(c(-0.5, 0.3), c(1, 0.4), c(4, 0.5))) {
      primary <- vapply(
        seq_len(n),
        function(i) tlogis_rng(0, pwindow, window[1], window[2]),
        numeric(1)
      )
      q <- unname(quantile(primary + rdist(n), probs))
      cdf <- vapply(
        q, primarycensored_numeric_cdf, numeric(1),
        case$dist_id, case$params, pwindow, 3L, window
      )
      expect_lt(
        max(abs(cdf - probs)), 0.03,
        label = tlogis_case_label(
          case,
          pwindow = pwindow, location = window[1], scale = window[2]
        )
      )
    }
  }
})

test_that("the ODE path with a truncated logistic primary has finite
  gradients matching finite differences", {
  model <- tlogis_gradient_model()
  # The integral is split at times that depend on the location and scale,
  # and the gradients through them must cancel
  cases <- list(tlogis_stan_cases[[4]], tlogis_stan_cases[[2]])
  points <- list(
    list(d = 3, pwindow = 2, location = 0.7, scale = 0.05),
    list(d = 1, pwindow = 2, location = 0.7, scale = 0.02),
    list(d = 3, pwindow = 2, location = -0.5, scale = 0.05),
    list(d = 4, pwindow = 2, location = 2.5, scale = 0.05)
  )
  for (case in cases) {
    for (point in points) {
      label <- tlogis_case_label(
        case,
        d = point$d, pwindow = point$pwindow, location = point$location,
        scale = point$scale
      )
      res <- tlogis_gradient_at(
        model, case, point$d, point$pwindow, point$location, point$scale
      )
      expect_false(res$gradient_not_finite, info = label)
      expect_false(res$rejected, info = label)
      expect_true(all(is.finite(res$gradient)), info = label)
      allowed <- 1e-2 * pmax(abs(res$finite_diff), 1e-2)
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
