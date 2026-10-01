skip_on_cran()

# Stan solutions for a truncated Gumbel primary (primary_id 4 with
# primary_params = [mu, beta]). Normal delays have the series, and the other
# delays use the numerical path.

gumbel_stan_cases <- list(
  list(
    dist_id = 4L, params = 60, pdist = pexp, args = list(rate = 60),
    d = c(0.005, 0.02, 0.1, 0.5, 1, 2.5)
  ),
  list(
    dist_id = 2L, params = c(3, 80), pdist = pgamma,
    args = list(shape = 3, rate = 80), d = c(0.01, 0.05, 0.1, 0.5, 1, 2.5)
  ),
  list(
    dist_id = 18L, params = c(3, 2), pdist = pnorm,
    args = list(mean = 3, sd = 2), d = c(-5, -1, 0.3, 1, 2.5, 5, 8, 12)
  ),
  list(
    dist_id = 18L, params = c(-1, 3), pdist = pnorm,
    args = list(mean = -1, sd = 3), d = c(-8, -1, 0.3, 1, 2.5, 4, 8)
  )
)

gumbel_series_cases <- gumbel_stan_cases[3:4]

gumbel_case_lower <- function(case) {
  if (case$dist_id == 18L) -Inf else 0
}

gumbel_case_cdf <- function(case) {
  function(x) do.call(case$pdist, c(list(x), case$args))
}

gumbel_case_label <- function(case, ...) {
  paste0(
    "dist ", case$dist_id, " params ", toString(case$params), ", ",
    paste(names(list(...)), unlist(list(...)), sep = " = ", collapse = ", ")
  )
}

# The delay object to call the R implementation with
gumbel_case_object <- function(case, mu, beta) {
  gumbel_object(
    list(pdist = case$pdist, args = case$args), mu, beta
  )
}

test_that("Stan tgumbel functions match the R functions", {
  settings <- list(
    list(min = 0, max = 1, mu = 0.5, beta = 0.3),
    list(min = 0, max = 2, mu = -0.5, beta = 0.2),
    list(min = 0, max = 2, mu = 1.5, beta = 0.1),
    list(min = 0, max = 1, mu = 0, beta = 1),
    list(min = 5, max = 7, mu = 6, beta = 0.5)
  )
  for (s in settings) {
    x <- seq(s$min, s$max, length.out = 9)
    info <- toString(unlist(s))
    stan_lpdf <- vapply(
      x, tgumbel_lpdf, numeric(1), s$min, s$max, s$mu, s$beta
    )
    expect_equal(
      stan_lpdf, dtgumbel(x, s$min, s$max, s$mu, s$beta, log = TRUE),
      tolerance = 1e-10, info = info
    )
    inner <- x[x > s$min & x < s$max]
    stan_lcdf <- vapply(
      inner, tgumbel_lcdf, numeric(1), s$min, s$max, s$mu, s$beta
    )
    expect_equal(
      stan_lcdf,
      ptgumbel(inner, s$min, s$max, s$mu, s$beta, log.p = TRUE),
      tolerance = 1e-10, info = info
    )
    expect_identical(tgumbel_lcdf(s$min - 1, s$min, s$max, s$mu, s$beta), -Inf)
    expect_identical(tgumbel_lcdf(s$max + 1, s$min, s$max, s$mu, s$beta), 0)
    expect_identical(
      tgumbel_lpdf(s$max + 1, s$min, s$max, s$mu, s$beta), -Inf
    )
  }
})

test_that("Stan tgumbel_rng matches the theoretical CDF", {
  for (s in list(
    list(min = 0, max = 1, mu = 0.5, beta = 0.3),
    list(min = 0, max = 2, mu = -0.5, beta = 0.2),
    list(min = 0, max = 2, mu = 1.5, beta = 0.1)
  )) {
    samples <- replicate(5000, tgumbel_rng(s$min, s$max, s$mu, s$beta))
    expect_true(all(samples >= s$min & samples <= s$max))
    ks <- suppressWarnings(stats::ks.test(
      samples,
      function(x) ptgumbel(x, s$min, s$max, s$mu, s$beta)
    ))
    expect_gt(ks$p.value, 0.001)
  }
})

test_that("gumbel_n_terms matches the R number of terms", {
  for (log_s0 in c(-10, -3, -1, 0, 0.5, 1, 2, log(15))) {
    expect_identical(
      gumbel_n_terms(log_s0), as.integer(.gumbel_n_terms(log_s0))
    )
  }
})

test_that("check_for_gumbel is structural", {
  expect_identical(check_for_gumbel(18L, 4L), 1L)
  expect_identical(check_for_gumbel(18L, 1L), 0L)
  expect_identical(check_for_gumbel(18L, 2L), 0L)
  expect_identical(check_for_analytical_params(18L, c(3, 2), 4L, c(0, 1)), 1L)
  for (dist_id in c(1L, 2L, 3L, 4L, 5L, 12L)) {
    expect_identical(check_for_gumbel(dist_id, 4L), 0L)
    expect_identical(check_for_analytical(dist_id, 4L), 0L)
  }
  # Non-parametric delays are not analytic for this primary
  expect_identical(check_for_analytical(26L, 4L), 0L)
})

test_that("check_for_gumbel_params needs a bounded window value and a small
  rounding estimate", {
  expect_identical(check_for_gumbel_params(18L, c(3, 2), 0, 1), 1L)
  expect_identical(check_for_gumbel_params(18L, c(3, 2), 1, 1), 1L)
  # exp(mu / beta) above 15 is not used
  expect_identical(check_for_gumbel_params(18L, c(3, 2), 2.8, 1), 0L)
  expect_identical(check_for_gumbel_params(18L, c(3, 2), 0.5, 0.1), 0L)
  expect_identical(check_for_gumbel_params(18L, c(3, 2), -5, 0.1), 1L)
  expect_identical(check_for_gumbel_params(18L, c(3, 2), 0, 0), 0L)
  expect_identical(check_for_gumbel_params(18L, c(3, 2), 0, -1), 0L)
  # A wide delay and a large mu / beta make the rounding estimate large
  expect_identical(check_for_gumbel_params(18L, c(0, 4), 1.15, 0.5), 0L)
  expect_identical(check_for_gumbel_params(18L, c(0, 0.5), 1.15, 0.5), 1L)
  # Only the normal delay has a series
  expect_identical(check_for_gumbel_params(4L, 80, 0, 1), 0L)
  expect_identical(check_for_gumbel_params(2L, c(2, 80), 0, 1), 0L)
  expect_identical(
    check_for_analytical_params(18L, c(3, 2), 4L, c(0, 1)), 1L
  )
  expect_identical(
    check_for_analytical_params(18L, c(3, 2), 4L, c(0.5, 0.1)), 0L
  )
  expect_identical(
    check_for_analytical_params(4L, 80, 4L, c(0, 1)), 0L
  )
  expect_identical(gumbel_screen_tolerance(), .gumbel_screen_tol)
  # Unchanged for the other solutions
  expect_identical(
    check_for_analytical_params(2L, c(2, 0.4), 1L, numeric(0)), 1L
  )
  expect_identical(
    check_for_analytical_params(2L, c(2, 0.4), 2L, 0.3), 1L
  )
})

test_that("Stan gumbel terms match the R transforms", {
  for (case in gumbel_series_cases) {
    obj <- gumbel_case_object(case, -0.5, 1)
    n_terms <- gumbel_n_terms(-0.5)
    ts <- c(0.01, 0.4, 1, 3.5)
    for (t in ts) {
      terms <- as.vector(primarycensored_gumbel_terms(
        t, case$dist_id, 1, n_terms, case$params
      ))
      expect_length(terms, 2 * (n_terms + 1))
      lower <- terms[seq(1, by = 2, length.out = n_terms + 1)]
      upper <- terms[seq(2, by = 2, length.out = n_terms + 1)]
      r_lower <- vapply(
        0:n_terms, function(n) .pcens_tilt_transform(obj, t, n), numeric(1)
      )
      r_upper <- vapply(
        0:n_terms,
        function(n) .pcens_tilt_transform(obj, t, n, upper = TRUE),
        numeric(1)
      )
      expect_equal(
        lower, r_lower,
        tolerance = 1e-10, info = gumbel_case_label(case, t = t)
      )
      expect_equal(
        upper, r_upper,
        tolerance = 1e-10, info = gumbel_case_label(case, t = t)
      )
    }
  }
})

mus <- c(-0.5, 0, 0.5, 1, 1.5)
betas <- c(0.1, 0.2, 1)
pwindows <- c(1, 2)

test_that("primarycensored_gumbel_lcdf matches the R series and the
  reference integral where the series applies", {
  for (case in gumbel_series_cases) {
    cdf <- gumbel_case_cdf(case)
    positive <- case$dist_id != 18L
    for (pwindow in pwindows) {
      for (beta in betas) {
        for (mu in mus) {
          if (!check_for_gumbel_params(case$dist_id, case$params, mu, beta)) {
            next
          }
          info <- gumbel_case_label(
            case,
            pwindow = pwindow, mu = mu, beta = beta
          )
          obj <- gumbel_case_object(case, mu, beta)
          n_terms <- .gumbel_n_terms(mu / beta)
          fit <- .gumbel_lcdf(
            obj, case$d, pwindow, mu, beta, n_terms,
            gumbel_case_lower(case)
          )
          accepted <- fit$error <= .gumbel_tol
          if (!any(accepted)) next
          d <- case$d[accepted]
          stan <- vapply(
            d, primarycensored_gumbel_lcdf, numeric(1),
            case$dist_id, case$params, pwindow, mu, beta
          )
          expect_equal(
            stan, fit$log_cdf[accepted],
            tolerance = 1e-9, info = info
          )
          expect_lt(
            max_rel_diff(
              exp(stan), gumbel_reference(d, pwindow, mu, beta, cdf, positive)
            ),
            1e-8,
            label = info
          )
        }
      }
    }
  }
})

test_that("the Stan series is accurate at extreme and narrow windows", {
  case <- gumbel_stan_cases[[3]]
  cdf <- gumbel_case_cdf(case)
  settings <- list(
    c(mu = -5, beta = 0.1, w = 1), c(mu = -50, beta = 1, w = 1),
    c(mu = 0, beta = 1, w = 0.05), c(mu = 0, beta = 5, w = 0.5),
    c(mu = 3, beta = 5, w = 2), c(mu = -1, beta = 0.3, w = 0.1)
  )
  for (s in settings) {
    expect_identical(
      check_for_gumbel_params(18L, case$params, s[["mu"]], s[["beta"]]), 1L
    )
    stan <- vapply(
      case$d, primarycensored_gumbel_lcdf, numeric(1),
      18L, case$params, s[["w"]], s[["mu"]], s[["beta"]]
    )
    expect_lt(
      max_rel_diff(
        exp(stan),
        gumbel_reference(
          case$d, s[["w"]], s[["mu"]], s[["beta"]], cdf, FALSE
        )
      ),
      1e-10,
      label = toString(s)
    )
  }
})

test_that("primarycensored_lcdf and primarycensored_cdf use the series and
  agree with the ODE path", {
  for (case in gumbel_stan_cases) {
    lower <- gumbel_case_lower(case)
    positive <- case$dist_id != 18L
    cdf <- gumbel_case_cdf(case)
    for (pwindow in pwindows) {
      for (beta in betas) {
        for (mu in mus) {
          info <- gumbel_case_label(
            case,
            pwindow = pwindow, mu = mu, beta = beta
          )
          analytic <- check_for_analytical_params(
            case$dist_id, case$params, 4L, c(mu, beta)
          ) == 1L
          lcdf <- vapply(
            case$d, primarycensored_lcdf, numeric(1),
            case$dist_id, case$params, pwindow, lower, Inf, 4L, c(mu, beta)
          )
          plain <- vapply(
            case$d, primarycensored_cdf, numeric(1),
            case$dist_id, case$params, pwindow, lower, Inf, 4L, c(mu, beta)
          )
          expect_equal(exp(lcdf), plain, tolerance = 1e-9, info = info)
          # Relative error where the reference is above 1e-8
          ode <- vapply(
            case$d, primarycensored_gumbel_numeric_cdf, numeric(1),
            case$dist_id, case$params, pwindow, mu, beta
          )
          reference <- gumbel_reference(
            case$d, pwindow, mu, beta, cdf, positive
          )
          expect_lt(
            gumbel_error(ode, reference, floor = 1e-8), 1e-7, label = info
          )
          accepted <- rep(FALSE, length(case$d))
          if (analytic) {
            fit <- .gumbel_lcdf(
              gumbel_case_object(case, mu, beta), case$d, pwindow, mu, beta,
              .gumbel_n_terms(mu / beta), lower
            )
            accepted <- fit$error <= .gumbel_tol
          }
          expect_lt(
            gumbel_error(plain, reference, floor = 1e-8), 1e-7, label = info
          )
          expect_lt(
            max(0, abs(plain - reference)[accepted]), 1e-8,
            label = info
          )
        }
      }
    }
  }
})

test_that("the numerical path is used where the series does not apply", {
  # exp(mu / beta) is about 148 here
  expect_identical(
    check_for_analytical_params(18L, c(3, 2), 4L, c(0.5, 0.1)), 0L
  )
  for (d in c(-1, 3, 8)) {
    ode <- primarycensored_gumbel_numeric_cdf(d, 18L, c(3, 2), 1, 0.5, 0.1)
    expect_identical(
      primarycensored_cdf(d, 18L, c(3, 2), 1, -Inf, Inf, 4L, c(0.5, 0.1)),
      ode
    )
    expect_equal(
      primarycensored_lcdf(d, 18L, c(3, 2), 1, -Inf, Inf, 4L, c(0.5, 0.1)),
      log(ode),
      tolerance = 1e-12
    )
  }
  # A rate that is too small for the tilts needed
  ode <- primarycensored_gumbel_numeric_cdf(2, 2L, c(2, 0.5), 2, 0, 1)
  expect_identical(
    primarycensored_cdf(2, 2L, c(2, 0.5), 2, 0, Inf, 4L, c(0, 1)),
    ode
  )
})

test_that("points where the series loses accuracy use the ODE", {
  # A window much narrower than the scale loses precision in the bracket
  obj <- gumbel_case_object(gumbel_stan_cases[[3]], 0, 1)
  d <- c(1, 3)
  n_terms <- gumbel_n_terms(0)
  fit <- .gumbel_lcdf(obj, d, 1e-7, 0, 1, n_terms, -Inf)
  expect_true(all(fit$error > .gumbel_tol))
  ode <- vapply(
    d, primarycensored_gumbel_numeric_cdf, numeric(1),
    18L, c(3, 2), 1e-7, 0, 1
  )
  expect_identical(
    vapply(
      d, primarycensored_gumbel_lcdf, numeric(1),
      18L, c(3, 2), 1e-7, 0, 1
    ),
    log(ode)
  )
})

# Delays long relative to the window, for which the series does not apply,
# with the matching family of gumbel_spike_families()
gumbel_spike_stan_cases <- list(
  list(dist_id = 18L, params = c(3, 2), family = 1L),
  list(dist_id = 4L, params = 1, family = 2L),
  list(dist_id = 2L, params = c(3, 1), family = 3L),
  list(dist_id = 1L, params = c(1, 0.5), family = 4L),
  list(dist_id = 3L, params = c(2, 3), family = 5L)
)

test_that("the ODE path resolves a narrow spike of the window density", {
  # mu / beta of 15 to 50, where integrating over the window misses the spike
  families <- gumbel_spike_families()
  for (case in gumbel_spike_stan_cases) {
    family <- families[[case$family]]
    cdf <- gumbel_cdf(family)
    for (s in gumbel_spike_settings()) {
      label <- gumbel_label(family, s[["w"]], s[["mu"]], s[["beta"]])
      reference <- gumbel_reference(
        family$q, s[["w"]], s[["mu"]], s[["beta"]], cdf, family$positive
      )
      primary <- c(s[["mu"]], s[["beta"]])
      ode <- vapply(
        family$q, primarycensored_gumbel_numeric_cdf, numeric(1),
        case$dist_id, case$params, s[["w"]], primary[1], primary[2]
      )
      expect_lt(
        gumbel_error(ode, reference, floor = 1e-8), 1e-7, label = label
      )
      plain <- vapply(
        family$q, primarycensored_cdf, numeric(1),
        case$dist_id, case$params, s[["w"]],
        if (family$positive) 0 else -Inf, Inf, 4L, primary
      )
      expect_lt(
        gumbel_error(plain, reference, floor = 1e-8), 1e-7, label = label
      )
      lcdf <- vapply(
        family$q, primarycensored_lcdf, numeric(1),
        case$dist_id, case$params, s[["w"]],
        if (family$positive) 0 else -Inf, Inf, 4L, primary
      )
      expect_false(anyNA(lcdf), label = label)
      expect_true(all(lcdf <= 0), label = label)
      keep <- reference > 1e-12
      expect_equal(
        exp(lcdf[keep]), reference[keep], tolerance = 1e-7, info = label
      )
    }
  }
})

test_that("the ODE path is accurate for a large mu over beta", {
  family <- gumbel_spike_families()[[1]]
  cdf <- gumbel_cdf(family)
  for (ratio in c(15, 20, 50, 200)) {
    beta <- 0.1
    mu <- ratio * beta
    for (pwindow in c(0.4, 1, 3)) {
      reference <- gumbel_reference(family$q, pwindow, mu, beta, cdf, FALSE)
      ode <- vapply(
        family$q, primarycensored_gumbel_numeric_cdf, numeric(1),
        18L, c(3, 2), pwindow, mu, beta
      )
      expect_lt(
        gumbel_error(ode, reference, floor = 1e-8), 1e-7,
        label = sprintf("mu over beta %g, pwindow %g", ratio, pwindow)
      )
    }
  }
})

test_that("the ODE path handles a location far below the window", {
  # s(pwindow) underflows, so the reference integrates dtgumbel()
  x <- c(1e-6, 0.3, 1, 3)
  for (s in list(
    c(mu = -50, beta = 0.02, w = 1), c(mu = -200, beta = 0.05, w = 2)
  )) {
    primary <- c(s[["mu"]], s[["beta"]])
    ode <- vapply(
      x, primarycensored_gumbel_numeric_cdf, numeric(1),
      4L, 1, s[["w"]], primary[1], primary[2]
    )
    reference <- vapply(x, function(d) {
      stats::integrate(
        function(z) {
          pexp(d - z, 1) * dtgumbel(z, 0, s[["w"]], s[["mu"]], s[["beta"]])
        },
        0, min(d, s[["w"]]), rel.tol = 1e-12
      )$value
    }, numeric(1))
    expect_true(all(ode >= 0 & ode <= 1), info = toString(s))
    expect_lt(
      gumbel_error(ode, reference, floor = 1e-8), 1e-7,
      label = toString(s)
    )
  }
})

test_that("the ODE log CDF is accurate far in the lower tail", {
  # The CDF is 1e-85 to 1e-225 here
  points <- list(
    c(mu = 2, beta = 0.1, w = 1),
    c(mu = 0.5, beta = 0.2, w = 2),
    c(mu = -3, beta = 0.5, w = 1),
    c(mu = 1.5, beta = 0.05, w = 7)
  )
  for (pt in points) {
    for (d in c(-35, -60)) {
      lcdf <- primarycensored_gumbel_numeric_lcdf(
        d, 18L, c(3, 2), pt[["w"]], pt[["mu"]], pt[["beta"]]
      )
      reference <- gumbel_reference(
        d, pt[["w"]], pt[["mu"]], pt[["beta"]],
        function(x) pnorm(x, 3, 2), FALSE
      )
      expect_equal(
        lcdf, log(reference), tolerance = 1e-8,
        info = paste(toString(pt), "d", d)
      )
    }
    # Beyond the double precision CDF the log CDF is finite and decreasing
    lcdf <- vapply(
      c(-60, -90, -120), primarycensored_gumbel_numeric_lcdf, numeric(1),
      18L, c(3, 2), pt[["w"]], pt[["mu"]], pt[["beta"]]
    )
    expect_true(all(is.finite(lcdf)), info = toString(pt))
    expect_true(all(diff(lcdf) < 0), info = toString(pt))
  }
})

test_that("a log CDF from the ODE is never NaN or above zero", {
  # Solver rounding must not give a NaN or positive log CDF
  points <- list(
    list(d = 31, id = 18L, par = c(3, 2), w = 1, mu = 2.7, beta = 0.05),
    list(d = 37, id = 18L, par = c(3, 2), w = 7, mu = 1.5, beta = 0.05),
    list(d = 0.01, id = 18L, par = c(1, 0.1), w = 1, mu = 0.5, beta = 0.2),
    list(d = 13, id = 18L, par = c(5, 1), w = 1, mu = 0.3, beta = 0.111),
    list(d = 12, id = 18L, par = c(5, 1), w = 1, mu = 0.3, beta = 0.111),
    list(d = 37, id = 4L, par = 1, w = 7, mu = 1.5, beta = 0.05),
    list(d = 37, id = 2L, par = c(3, 1), w = 7, mu = 1.5, beta = 0.05)
  )
  for (pt in points) {
    lower <- if (pt$id == 18L) -Inf else 0
    lcdf <- primarycensored_lcdf(
      pt$d, pt$id, pt$par, pt$w, lower, Inf, 4L, c(pt$mu, pt$beta)
    )
    info <- toString(unlist(pt))
    expect_false(is.nan(lcdf), info = info)
    expect_lte(lcdf, 0)
    cdf <- gumbel_reference(
      pt$d, pt$w, pt$mu, pt$beta,
      function(x) {
        if (pt$id == 18L) {
          pnorm(x, pt$par[1], pt$par[2])
        } else if (pt$id == 4L) {
          pexp(x, pt$par)
        } else {
          pgamma(x, pt$par[1], pt$par[2])
        }
      },
      pt$id != 18L
    )
    expect_equal(exp(lcdf), cdf, tolerance = 1e-7, info = info)
  }
})

test_that("the series is used where its estimate is below the tolerance", {
  expect_identical(gumbel_error_tolerance(), .gumbel_tol)
  case <- gumbel_stan_cases[[3]]
  used <- 0
  for (mu in c(0.2, 0.25, 0.27, 0.3)) {
    beta <- 0.1
    n_terms <- gumbel_n_terms(mu / beta)
    for (d in c(-2, 0.5, 2, 5, 9, 14)) {
      fit <- as.vector(primarycensored_gumbel_lcdf_from_terms(
        primarycensored_gumbel_terms(d, 18L, beta, n_terms, case$params),
        primarycensored_gumbel_terms(d - 1, 18L, beta, n_terms, case$params),
        d, 1, mu, beta, n_terms
      ))
      lcdf <- primarycensored_gumbel_lcdf(d, 18L, case$params, 1, mu, beta)
      info <- sprintf("mu %g, d %g", mu, d)
      if (fit[2] <= gumbel_error_tolerance()) {
        expect_identical(lcdf, fit[1], info = info)
        used <- used + 1
      } else {
        expect_identical(
          lcdf,
          primarycensored_gumbel_numeric_lcdf(
            d, 18L, case$params, 1, mu, beta
          ),
          info = info
        )
      }
      expect_equal(
        exp(lcdf),
        gumbel_reference(d, 1, mu, beta, gumbel_case_cdf(case), FALSE),
        tolerance = 1e-8, info = info
      )
    }
  }
  expect_gt(used, 0)
})

test_that("Stan tgumbel is accurate for a spike and a far location", {
  # s(xmax) = exp(25), so the density needs the difference of s
  t <- c(0, 1e-13, 1e-12, 5e-12)
  x <- 1 - t
  expect_equal(
    vapply(x, tgumbel_lpdf, numeric(1), 0, 1, 1.5, 0.02),
    dtgumbel(x, 0, 1, 1.5, 0.02, log = TRUE),
    tolerance = 1e-12
  )
  expect_equal(
    tgumbel_lpdf(1, 0, 1, 1.5, 0.02), 25 - log(0.02),
    tolerance = 1e-13
  )
  # The normalisation underflows for a location far below the window
  x <- c(0, 0.05, 0.5, 1)
  for (beta in c(0.02, 0.05)) {
    expect_equal(
      vapply(x, tgumbel_lpdf, numeric(1), 0, 1, -50, beta),
      dexpgrowth(x, 0, 1, r = -1 / beta, log = TRUE),
      tolerance = 1e-10
    )
    inner <- x[2:3]
    expect_equal(
      vapply(inner, tgumbel_lcdf, numeric(1), 0, 1, -50, beta),
      pexpgrowth(inner, 0, 1, r = -1 / beta, log.p = TRUE),
      tolerance = 1e-10
    )
  }
  draws <- replicate(50, tgumbel_rng(0, 1, 1.5, 0.02))
  expect_true(all(draws <= 1 & draws > 1 - 1e-8))
  draws <- replicate(50, tgumbel_rng(0, 1, -50, 0.02))
  expect_true(all(is.finite(draws) & draws >= 0 & draws <= 1))
})

test_that("the analytical function rejects an inadmissible primary", {
  expect_error(
    primarycensored_analytical_lcdf(
      2, 18L, c(3, 2), 1, -Inf, Inf, 4L, c(0.5, 0.1)
    ),
    "truncated Gumbel"
  )
})

test_that("the analytical lcdf applies truncation for a normal delay", {
  case <- gumbel_stan_cases[[3]]
  pwindow <- 2
  mu <- 0
  beta <- 1
  cdf_at <- function(x) {
    exp(primarycensored_lcdf(
      x, 18L, case$params, pwindow, -Inf, Inf, 4L, c(mu, beta)
    ))
  }
  d <- 4
  L <- 1
  D <- 9
  expected <- (cdf_at(d) - cdf_at(L)) / (cdf_at(D) - cdf_at(L))
  expect_equal(
    exp(primarycensored_lcdf(
      d, 18L, case$params, pwindow, L, D, 4L, c(mu, beta)
    )),
    expected,
    tolerance = 1e-9
  )
})

test_that("the vectorised Gumbel solution needs an integer pwindow", {
  expect_identical(check_for_analytical_vectorized(18L, 4L, 1), 1L)
  expect_identical(check_for_analytical_vectorized(18L, 4L, 7), 1L)
  expect_identical(check_for_analytical_vectorized(18L, 4L, 1.5), 0L)
  expect_identical(check_for_analytical_vectorized(18L, 4L, 0.5), 0L)
  expect_identical(check_for_analytical_vectorized(18L, 3L, 1), 0L)
  for (dist_id in c(1L, 2L, 3L, 4L, 26L)) {
    expect_identical(check_for_analytical_vectorized(dist_id, 4L, 1), 0L)
  }
})

per_delay_gumbel_lcdf <- function(delays, case, pwindow, mu, beta) {
  vapply(
    delays, primarycensored_lcdf, numeric(1), # nolint: object_usage_linter.
    case$dist_id, case$params, pwindow, gumbel_case_lower(case), Inf, 4L,
    c(mu, beta)
  )
}

test_that("the vectorised Gumbel CDF matches the per delay CDF", {
  n <- 21L
  for (case in gumbel_series_cases) {
    for (pwindow in c(1, 2, 5)) {
      for (beta in c(0.2, 1)) {
        for (mu in c(-1, 0, 1)) {
          if (!check_for_gumbel_params(case$dist_id, case$params, mu, beta)) {
            next
          }
          for (start in c(1L, 5L)) {
            vectorised <- primarycensored_analytical_lcdf_vectorized(
              start, n, case$dist_id, case$params, pwindow, 4L, c(mu, beta)
            )
            expect_length(vectorised, n)
            expect_identical(
              vectorised[start:n],
              per_delay_gumbel_lcdf(start:n, case, pwindow, mu, beta),
              info = gumbel_case_label(
                case,
                pwindow = pwindow, mu = mu, beta = beta, start = start
              )
            )
          }
        }
      }
    }
  }
})

test_that("primarycensored_lcdf_vectorized uses the shared Gumbel terms", {
  for (case in gumbel_series_cases) {
    expect_identical(
      primarycensored_lcdf_vectorized(
        1L, 15L, case$dist_id, case$params, 3, 4L, c(-1, 1)
      ),
      primarycensored_analytical_lcdf_vectorized(
        1L, 15L, case$dist_id, case$params, 3, 4L, c(-1, 1)
      )
    )
  }
  # A non-integer window and inadmissible parameters use the per delay path
  case <- gumbel_stan_cases[[3]]
  expect_identical(
    primarycensored_lcdf_vectorized(
      1L, 10L, case$dist_id, case$params, 1.5, 4L, c(0, 1)
    ),
    per_delay_gumbel_lcdf(1:10, case, 1.5, 0, 1)
  )
  expect_identical(
    primarycensored_lcdf_vectorized(
      1L, 10L, case$dist_id, case$params, 2, 4L, c(0.5, 0.1)
    ),
    per_delay_gumbel_lcdf(1:10, case, 2, 0.5, 0.1)
  )
})

test_that("the vectorised PMF matches the per delay PMF with truncation", {
  settings <- list(
    list(max_delay = 10, L = 0, D = 11),
    list(max_delay = 10, L = 0, D = Inf),
    list(max_delay = 20, L = 2, D = 30),
    list(max_delay = 20, L = -Inf, D = Inf)
  )
  for (case in gumbel_stan_cases[c(1, 3)]) {
    lower_support <- case$dist_id == 18L
    for (setting in settings) {
      if (lower_support && is.finite(setting$L)) next
      for (pwindow in c(1, 3)) {
        info <- gumbel_case_label(
          case,
          pwindow = pwindow, L = setting$L, D = setting$D
        )
        vectorised <- primarycensored_sone_lpmf_vectorized(
          setting$max_delay, setting$L, setting$D, case$dist_id, case$params,
          pwindow, 4L, c(-1, 1)
        )
        expected <- vapply(
          seq_len(setting$max_delay + 1L) - 1L,
          function(d) {
            primarycensored_lpmf(
              d, case$dist_id, case$params, pwindow, d + 1, setting$L,
              setting$D, 4L, c(-1, 1)
            )
          },
          numeric(1)
        )
        expect_equal(vectorised, expected, tolerance = 1e-9, info = info)
      }
    }
  }
})

# Gradients are only observable from a compiled model, so this builds a
# minimal one whose target is the log CDF or the vectorised log PMF, and runs
# `stan_gradient_at()` from helper-stan-gradient.R.
gumbel_gradient_model <- function() {
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
    "  real mu;\n",
    "  real<lower=0> beta;\n",
    "}\n",
    "model {\n",
    "  array[2] real all_params = {p1, p2};\n",
    "  array[n_params] real params = all_params[1:n_params];\n",
    "  if (vectorised) {\n",
    "    target += sum(primarycensored_sone_lpmf_vectorized(\n",
    "      to_int(d), L, positive_infinity(), dist_id, params, pwindow, 4,\n",
    "      {mu, beta}\n",
    "    ));\n",
    "  } else {\n",
    "    target += primarycensored_lcdf(\n",
    "      d | dist_id, params, pwindow, L, positive_infinity(), 4,\n",
    "      {mu, beta}\n",
    "    );\n",
    "  }\n",
    "}\n"
  )
  path <- file.path(tempdir(), "pcd_gumbel_gradient.stan")
  writeLines(code, path)
  suppressMessages(suppressWarnings(cmdstanr::cmdstan_model(path)))
}

gumbel_gradient_at <- function(model, case, d, pwindow, mu, beta,
                               vectorised = FALSE) {
  init <- list(
    p1 = case$params[1],
    p2 = if (length(case$params) > 1) case$params[2] else 1,
    mu = mu, beta = beta
  )
  stan_gradient_at( # nolint: object_usage_linter.
    model,
    data = list(
      dist_id = case$dist_id, n_params = length(case$params),
      vectorised = as.integer(vectorised), d = d, pwindow = pwindow,
      L = gumbel_case_lower(case)
    ),
    init = init
  )
}

# Stan's gamma shape gradient has a relative error of about 1e-3
expect_gumbel_gradient_close <- function(res, case, label, scale = 1,
                                         slack = 0) {
  tolerance <- rep(1e-4, 4) * scale
  if (case$dist_id == 2L) {
    tolerance[1] <- 2e-2
  }
  allowed <- tolerance * pmax(abs(res$finite_diff), 1e-2) + slack
  testthat::expect_true(
    all(abs(res$gradient - res$finite_diff) <= allowed),
    info = paste0(
      label, ": gradient ", toString(signif(res$gradient, 5)),
      ", finite difference ", toString(signif(res$finite_diff, 5))
    )
  )
}

test_that("Gumbel log CDFs have finite gradients matching finite
  differences", {
  model <- gumbel_gradient_model()
  points <- list(
    list(d = 0.3, pwindow = 2, mu = -0.5, beta = 0.5),
    list(d = 1, pwindow = 1, mu = 0, beta = 1),
    list(d = 2.5, pwindow = 2, mu = 0.5, beta = 1),
    list(d = 6, pwindow = 3, mu = -1, beta = 0.3),
    list(d = 20, pwindow = 3, mu = 1, beta = 1),
    list(d = -3, pwindow = 2, mu = 0, beta = 0.5),
    list(d = 1.5, pwindow = 1, mu = -5, beta = 0.1),
    list(d = 2, pwindow = 2, mu = -30, beta = 1),
    list(d = 3, pwindow = 1, mu = 0, beta = 8),
    # The ODE is used
    list(d = 4, pwindow = 2, mu = 0.5, beta = 0.1),
    # Narrow spikes of the window density
    list(d = 1, pwindow = 1, mu = 2, beta = 0.1),
    list(d = 5, pwindow = 1, mu = 1.5, beta = 0.05),
    list(d = 7, pwindow = 7, mu = 1.5, beta = 0.05),
    list(d = 3, pwindow = 2, mu = 3, beta = 0.06),
    # Around where the series is replaced by the ODE
    list(d = 2, pwindow = 1, mu = 0.25, beta = 0.1),
    list(d = 2, pwindow = 1, mu = 0.27, beta = 0.1),
    list(d = 5, pwindow = 1, mu = 0.3, beta = 0.111),
    list(d = 11, pwindow = 1, mu = 0.3, beta = 0.111)
  )
  cases <- gumbel_stan_cases[c(1, 3, 4)]
  for (case in cases) {
    for (point in points) {
      if (case$dist_id != 18L && point$d <= 0) next
      label <- gumbel_case_label(
        case,
        d = point$d, pwindow = point$pwindow, mu = point$mu, beta = point$beta
      )
      res <- gumbel_gradient_at(
        model, case, point$d, point$pwindow, point$mu, point$beta
      )
      expect_false(res$gradient_not_finite, info = label)
      expect_false(res$rejected, info = label)
      expect_length(res$gradient, 4)
      expect_true(all(is.finite(res$gradient)), info = label)
      expect_gumbel_gradient_close(res, case, label)
    }
  }
})

test_that("Gumbel log CDF gradients match finite differences across the
  acceptance grid", {
  model <- gumbel_gradient_model()
  cases <- list(
    list(dist_id = 18L, params = c(0.5, 1), d = c(0.02, 1, 3)),
    list(dist_id = 18L, params = c(3, 2), d = c(0.02, 1, 3, 6))
  )
  for (case in cases) {
    for (mu in c(-0.5, 0, 0.5, 1, 1.5)) {
      for (beta in c(0.1, 0.2, 1)) {
        for (pwindow in c(1, 2)) {
          for (d in case$d) {
            label <- gumbel_case_label(
              case, d = d, pwindow = pwindow, mu = mu, beta = beta
            )
            # Finite differences of a CDF below 1e-260 are not meaningful
            lcdf <- primarycensored_lcdf(
              d, case$dist_id, case$params, pwindow,
              gumbel_case_lower(case), Inf, 4L, c(mu, beta)
            )
            if (lcdf < -600) next
            res <- gumbel_gradient_at(model, case, d, pwindow, mu, beta)
            expect_false(res$gradient_not_finite, info = label)
            expect_false(res$rejected, info = label)
            expect_length(res$gradient, 4)
            # Finite differences have an absolute error of about 3e-5, and
            # a relative error of about 1e-4 from the solver on the
            # numerical path
            expect_gumbel_gradient_close(
              res, case, label, scale = 10, slack = 1e-4
            )
          }
        }
      }
    }
  }
})

test_that("Gumbel gradients are finite for delays whose CDF has an
  infinite slope at 0", {
  # A gamma or Weibull with shape below 1 has an unbounded density at 0
  model <- gumbel_gradient_model()
  cases <- list(
    list(dist_id = 2L, params = c(0.5, 1)),
    list(dist_id = 2L, params = c(0.8, 1)),
    list(dist_id = 2L, params = c(0.3, 2)),
    list(dist_id = 2L, params = c(0.5, 200)),
    list(dist_id = 3L, params = c(0.5, 1))
  )
  points <- list(
    list(d = 0.5, pwindow = 2, mu = 0, beta = 0.5),
    list(d = 0.2, pwindow = 1, mu = -1, beta = 0.3),
    list(d = 1, pwindow = 1, mu = 0.5, beta = 0.3)
  )
  for (case in cases) {
    for (point in points) {
      label <- gumbel_case_label(
        case,
        d = point$d, pwindow = point$pwindow, mu = point$mu, beta = point$beta
      )
      elapsed <- system.time({
        res <- gumbel_gradient_at(
          model, case, point$d, point$pwindow, point$mu, point$beta
        )
      })[["elapsed"]]
      expect_false(res$gradient_not_finite, info = label)
      expect_false(res$rejected, info = label)
      expect_length(res$gradient, 4)
      expect_true(all(is.finite(res$gradient)), info = label)
      expect_lt(elapsed, 5, label = label)
      expect_gumbel_gradient_close(res, case, label, scale = 5)
    }
  }
})

test_that("the vectorised Gumbel log PMF has finite gradients for delays
  whose CDF has an infinite slope at 0", {
  model <- gumbel_gradient_model()
  cases <- list(
    list(dist_id = 2L, params = c(0.5, 1)),
    list(dist_id = 3L, params = c(0.5, 1))
  )
  points <- list(
    list(d = 10, pwindow = 1, mu = 0.5, beta = 0.3),
    list(d = 4, pwindow = 2, mu = 0, beta = 0.5)
  )
  for (case in cases) {
    for (point in points) {
      label <- gumbel_case_label(
        case,
        d = point$d, pwindow = point$pwindow, mu = point$mu, beta = point$beta
      )
      elapsed <- system.time({
        res <- gumbel_gradient_at(
          model, case, point$d, point$pwindow, point$mu, point$beta,
          vectorised = TRUE
        )
      })[["elapsed"]]
      expect_false(res$gradient_not_finite, info = label)
      expect_false(res$rejected, info = label)
      expect_true(all(is.finite(res$gradient)), info = label)
      expect_lt(elapsed, 5, label = label)
      expect_gumbel_gradient_close(res, case, label, scale = 5)
    }
  }
})

test_that("lognormal and Weibull gradients are finite for a narrow window", {
  model <- gumbel_gradient_model()
  cases <- list(
    list(dist_id = 1L, params = c(1, 0.5)),
    list(dist_id = 3L, params = c(2, 3)),
    list(dist_id = 2L, params = c(4, 2))
  )
  points <- list(
    list(d = 1.9, pwindow = 2, mu = 2.1, beta = 0.1),
    list(d = 1, pwindow = 2, mu = 2, beta = 0.5),
    list(d = 1, pwindow = 2, mu = 2.2, beta = 0.2),
    list(d = 4, pwindow = 2, mu = 3, beta = 0.1)
  )
  for (case in cases) {
    for (point in points) {
      label <- gumbel_case_label(
        case,
        d = point$d, pwindow = point$pwindow, mu = point$mu, beta = point$beta
      )
      elapsed <- system.time({
        res <- gumbel_gradient_at(
          model, case, point$d, point$pwindow, point$mu, point$beta
        )
      })[["elapsed"]]
      expect_false(res$gradient_not_finite, info = label)
      expect_false(res$rejected, info = label)
      expect_true(all(is.finite(res$gradient)), info = label)
      expect_lt(elapsed, 5, label = label)
      expect_gumbel_gradient_close(res, case, label, scale = 5)
    }
  }
})

test_that("vectorised lognormal and Weibull gradients are finite", {
  model <- gumbel_gradient_model()
  cases <- list(
    list(dist_id = 1L, params = c(1, 0.5)),
    list(dist_id = 3L, params = c(2, 3))
  )
  for (case in cases) {
    label <- gumbel_case_label(case, d = 6, pwindow = 2, mu = 2.1, beta = 0.3)
    elapsed <- system.time({
      res <- gumbel_gradient_at(
        model, case, 6, 2, 2.1, 0.3, vectorised = TRUE
      )
    })[["elapsed"]]
    expect_false(res$gradient_not_finite, info = label)
    expect_false(res$rejected, info = label)
    expect_true(all(is.finite(res$gradient)), info = label)
    expect_lt(elapsed, 5, label = label)
    expect_gumbel_gradient_close(res, case, label, scale = 5)
  }
})

test_that("a log CDF below -1000 is -inf in the numerical path", {
  # The window mass is at 2, so a delay of at most 1 needs u of exp(20)
  for (case in list(
    list(id = 2L, par = c(4, 2)), list(id = 1L, par = c(1, 0.5)),
    list(id = 4L, par = 1)
  )) {
    expect_identical(
      primarycensored_gumbel_numeric_lcdf(
        1, case$id, case$par, 2, 3, 0.1
      ),
      -Inf
    )
    expect_identical(
      primarycensored_lcdf(1, case$id, case$par, 2, 0, Inf, 4L, c(3, 0.1)),
      -Inf
    )
  }
})

test_that("the Stan numerical path keeps the lower tail on the log scale", {
  # References are an independent log scale integral and R
  cases <- list(
    list(id = 4L, par = 1, d = 1 / 3, w = 1, mu = 1.5, beta = 0.3,
         ref = -48.6807),
    list(id = 4L, par = 1, d = 1, w = 3, mu = 1.5, beta = 0.1,
         ref = -155.72),
    list(id = 2L, par = c(5, 1), d = 2, w = 4, mu = 3, beta = 0.3,
         ref = -51.1984),
    list(id = 18L, par = c(1, 0.1), d = 1, w = 3, mu = 0.9, beta = 0.1,
         ref = -32.9165),
    list(id = 18L, par = c(0, 1), d = -10, w = 10, mu = 45, beta = 10,
         ref = -108.60),
    list(id = 18L, par = c(0, 1), d = -3, w = 10, mu = 45, beta = 10,
         ref = -51.684)
  )
  for (pt in cases) {
    info <- toString(unlist(pt[c("id", "par", "d", "w", "mu", "beta")]))
    lower <- if (pt$id == 18L) -Inf else 0
    lcdf <- primarycensored_lcdf(
      pt$d, pt$id, pt$par, pt$w, lower, Inf, 4L, c(pt$mu, pt$beta)
    )
    expect_equal(lcdf, pt$ref, tolerance = 1e-4, info = info)
    family <- switch(
      as.character(pt$id),
      "4" = list(pdist = pexp, args = list(rate = pt$par)),
      "2" = list(
        pdist = pgamma, args = list(shape = pt$par[1], rate = pt$par[2])
      ),
      "18" = list(
        pdist = pnorm, args = list(mean = pt$par[1], sd = pt$par[2])
      )
    )
    r_cdf <- pcens_cdf(
      gumbel_object(family, pt$mu, pt$beta), pt$d, pt$w, use_numeric = TRUE
    )
    expect_equal(lcdf, log(r_cdf), tolerance = 1e-8, info = info)
    numeric_lcdf <- primarycensored_gumbel_numeric_lcdf(
      pt$d, pt$id, pt$par, pt$w, pt$mu, pt$beta
    )
    expect_equal(numeric_lcdf, log(r_cdf), tolerance = 1e-8, info = info)
  }
  # The log PMF of a day scale gamma delay is finite
  lpmf <- primarycensored_lpmf(
    1, 2L, c(5, 1), 4, 2, 0, Inf, 4L, c(3, 0.3)
  )
  expect_true(is.finite(lpmf))
})

test_that("the vectorised Gumbel log PMF has finite gradients matching
  finite differences", {
  model <- gumbel_gradient_model()
  points <- list(
    list(d = 12, pwindow = 3, mu = -0.5, beta = 0.5),
    list(d = 8, pwindow = 2, mu = 0, beta = 1),
    list(d = 10, pwindow = 1, mu = 1, beta = 1)
  )
  # Only normal delays, as the PMF of a short exponential delay is a
  # difference of CDFs that are 1 in double precision. Finite differences of
  # a sum of PMFs have an error of about 5e-4.
  for (case in gumbel_stan_cases[c(3, 4)]) {
    for (point in points) {
      label <- gumbel_case_label(
        case,
        d = point$d, pwindow = point$pwindow, mu = point$mu, beta = point$beta
      )
      res <- gumbel_gradient_at(
        model, case, point$d, point$pwindow, point$mu, point$beta,
        vectorised = TRUE
      )
      expect_false(res$gradient_not_finite, info = label)
      expect_false(res$rejected, info = label)
      expect_true(all(is.finite(res$gradient)), info = label)
      expect_gumbel_gradient_close(res, case, label, scale = 5)
    }
  }
})

# The compiled model is shared by the tests that sample from it
gumbel_pcens_model_cache <- new.env()
gumbel_pcens_model <- function() {
  if (is.null(gumbel_pcens_model_cache$model)) {
    gumbel_pcens_model_cache$model <- suppressMessages(
      suppressWarnings(pcd_cmdstan_model())
    )
  }
  gumbel_pcens_model_cache$model
}

test_that(
  "pcd_cmdstan_model recovers true values for a normal delay with a
   truncated Gumbel primary",
  {
    testthat::skip_if_not_installed("cmdstanr")
    testthat::skip_if_not_installed("dplyr")
    testthat::skip_if(
      is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))
    )
    set.seed(321)
    n <- 2000
    true_mean <- 1
    true_sd <- 2

    simulated_delays <- rprimarycensored(
      n = n,
      rdist = rnorm,
      mean = true_mean,
      sd = true_sd,
      pwindow = 2,
      L = -5,
      D = 8,
      rprimary = rtgumbel,
      rprimary_args = list(mu = -0.5, beta = 1)
    )
    delay_counts <- dplyr::summarise(
      data.frame(
        delay = simulated_delays,
        delay_upper = simulated_delays + 1,
        pwindow = 2,
        start_relative_obs_time = -5,
        relative_obs_time = 8
      ),
      n = dplyr::n(),
      .by = c(
        pwindow, start_relative_obs_time, relative_obs_time,
        delay, delay_upper
      )
    )

    # The normal with a truncated Gumbel primary has a series solution, so
    # this fit uses it with the vectorised shared terms. The delay mean and
    # the window are confounded, so the window has informative priors.
    stan_data <- pcd_as_stan_data(
      delay_counts,
      dist_id = pcd_stan_dist_id("normal", "delay"),
      primary_id = pcd_stan_dist_id("tgumbel", "primary"),
      param_bounds = list(lower = c(-Inf, 0.01), upper = c(Inf, Inf)),
      primary_param_bounds = list(lower = c(-Inf, 0.05), upper = c(Inf, Inf)),
      priors = list(location = c(0, 0), scale = c(5, 2.5)),
      primary_priors = list(location = c(-0.5, 1), scale = c(0.1, 0.1))
    )

    model <- gumbel_pcens_model()
    fit <- suppressMessages(suppressWarnings(model$sample(
      data = stan_data,
      seed = 321,
      chains = 2,
      parallel_chains = 2,
      refresh = 0,
      show_messages = FALSE,
      iter_warmup = 500,
      iter_sampling = 500
    )))

    posterior <- fit$draws(c("params[1]", "params[2]"), format = "df")
    expect_equal(mean(posterior$`params[1]`), true_mean, tolerance = 0.1)
    expect_equal(mean(posterior$`params[2]`), true_sd, tolerance = 0.1)
    ci_mean <- quantile(posterior$`params[1]`, c(0.05, 0.95))
    ci_sd <- quantile(posterior$`params[2]`, c(0.05, 0.95))
    expect_gt(true_mean, ci_mean[1])
    expect_lt(true_mean, ci_mean[2])
    expect_gt(true_sd, ci_sd[1])
    expect_lt(true_sd, ci_sd[2])
  }
)
