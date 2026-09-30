skip_on_cran()

# Stan solutions for exponential, gamma and normal delays with a truncated
# Gumbel primary (primary_id 4 with primary_params = [mu, beta]). These tests
# check the Stan primary functions against R, the series against the R
# implementation, a reference integral and the ODE path, the dispatch and the
# fallback, the shared endpoint vectorised form, and gradients.

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
  for (dist_id in c(2L, 4L, 18L)) {
    expect_identical(check_for_gumbel(dist_id, 4L), 1L)
    expect_identical(check_for_gumbel(dist_id, 1L), 0L)
    expect_identical(check_for_gumbel(dist_id, 2L), 0L)
    expect_identical(check_for_analytical(dist_id, 4L), 1L)
  }
  for (dist_id in c(1L, 3L, 5L, 12L)) {
    expect_identical(check_for_gumbel(dist_id, 4L), 0L)
    expect_identical(check_for_analytical(dist_id, 4L), 0L)
  }
  # Non-parametric delays are not analytic for this primary
  expect_identical(check_for_analytical(26L, 4L), 0L)
})

test_that("check_for_gumbel_params needs a bounded window value and the
  tilted delay", {
  # exp(mu / beta) above 15 is not used
  expect_identical(check_for_gumbel_params(18L, c(3, 2), 0, 1), 1L)
  expect_identical(check_for_gumbel_params(18L, c(3, 2), 2.7, 1), 1L)
  expect_identical(check_for_gumbel_params(18L, c(3, 2), 2.8, 1), 0L)
  expect_identical(check_for_gumbel_params(18L, c(3, 2), 0.5, 0.1), 0L)
  expect_identical(check_for_gumbel_params(18L, c(3, 2), -5, 0.1), 1L)
  expect_identical(check_for_gumbel_params(18L, c(3, 2), 0, 0), 0L)
  expect_identical(check_for_gumbel_params(18L, c(3, 2), 0, -1), 0L)
  # The delay needs a rate above n_terms / beta
  n_terms <- gumbel_n_terms(0)
  expect_identical(
    check_for_gumbel_params(4L, n_terms + 1, 0, 1), 1L
  )
  expect_identical(check_for_gumbel_params(4L, n_terms - 1, 0, 1), 0L)
  expect_identical(check_for_gumbel_params(4L, 1, 0, 1), 0L)
  expect_identical(check_for_gumbel_params(2L, c(2, 1), 0, 1), 0L)
  expect_identical(check_for_gumbel_params(2L, c(2, 80), 0, 1), 1L)
  expect_identical(
    check_for_analytical_params(18L, c(3, 2), 4L, c(0, 1)), 1L
  )
  expect_identical(
    check_for_analytical_params(18L, c(3, 2), 4L, c(0.5, 0.1)), 0L
  )
  expect_identical(
    check_for_analytical_params(4L, 1, 4L, c(0, 1)), 0L
  )
  # Unchanged for the other solutions
  expect_identical(
    check_for_analytical_params(2L, c(2, 0.4), 1L, numeric(0)), 1L
  )
  expect_identical(
    check_for_analytical_params(2L, c(2, 0.4), 2L, 0.3), 1L
  )
})

test_that("Stan gumbel terms match the R transforms", {
  for (case in gumbel_stan_cases) {
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
  for (case in gumbel_stan_cases) {
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

test_that("Stan and R accept the same points of the series", {
  case <- gumbel_stan_cases[[3]]
  for (beta in betas) {
    for (mu in mus) {
      if (!check_for_gumbel_params(case$dist_id, case$params, mu, beta)) next
      n_terms <- gumbel_n_terms(mu / beta)
      obj <- gumbel_case_object(case, mu, beta)
      fit <- .gumbel_lcdf(obj, case$d, 1, mu, beta, n_terms, -Inf)
      stan_error <- vapply(case$d, function(d) {
        as.vector(primarycensored_gumbel_lcdf_from_terms(
          primarycensored_gumbel_terms(d, 18L, beta, n_terms, case$params),
          primarycensored_gumbel_terms(d - 1, 18L, beta, n_terms, case$params),
          d, 1, mu, beta, n_terms
        ))[2]
      }, numeric(1))
      info <- gumbel_case_label(case, mu = mu, beta = beta)
      expect_identical(
        stan_error <= gumbel_error_tolerance(), fit$error <= .gumbel_tol,
        info = info
      )
      expect_equal(
        log(stan_error), log(fit$error),
        tolerance = 1e-3, info = info
      )
    }
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
          # The ODE path is used where the series does not apply, and is
          # accurate to 1e-8 of the CDF, relative where it is above 1e-12
          ode <- vapply(
            case$d, primarycensored_numeric_cdf, numeric(1),
            case$dist_id, case$params, pwindow, 4L, c(mu, beta)
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
          # The series is accurate to 1e-8 and so is the ODE path, at every
          # mu / beta, where the density can be a narrow spike
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
    ode <- primarycensored_numeric_cdf(d, 18L, c(3, 2), 1, 4L, c(0.5, 0.1))
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
  ode <- primarycensored_numeric_cdf(2, 2L, c(2, 0.5), 2, 4L, c(0, 1))
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
    d, primarycensored_numeric_cdf, numeric(1),
    18L, c(3, 2), 1e-7, 4L, c(0, 1)
  )
  expect_identical(
    vapply(
      d, primarycensored_gumbel_lcdf, numeric(1),
      18L, c(3, 2), 1e-7, 0, 1
    ),
    log(ode)
  )
})

# The ODE path is accurate to about 1e-15 absolute for a small CDF, so it is
# compared with a floor of 1e-8 on the reference, see gumbel_error().

# Delays long relative to the window, for which the series does not apply,
# with the matching family of gumbel_spike_families()
gumbel_spike_stan_cases <- list(
  list(dist_id = 18L, params = c(3, 2), family = 1L),
  list(dist_id = 4L, params = 1, family = 2L),
  list(dist_id = 2L, params = c(3, 1), family = 3L)
)

test_that("the ODE path resolves a narrow spike of the window density", {
  # mu / beta of 15 to 50 at the end of the window and inside it, where
  # integrating over the window returns 0 or more than 1 without an error
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
        family$q, primarycensored_numeric_cdf, numeric(1),
        case$dist_id, case$params, s[["w"]], 4L, primary
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
        family$q, primarycensored_numeric_cdf, numeric(1),
        18L, c(3, 2), pwindow, 4L, c(mu, beta)
      )
      expect_lt(
        gumbel_error(ode, reference, floor = 1e-8), 1e-7,
        label = sprintf("mu over beta %g, pwindow %g", ratio, pwindow)
      )
    }
  }
})

test_that("the ODE path handles a location far below the window", {
  # s(pwindow) underflows and the density is the exponentially decaying one,
  # so the reference integrates dtgumbel() over the window
  x <- c(1e-6, 0.3, 1, 3)
  for (s in list(
    c(mu = -50, beta = 0.02, w = 1), c(mu = -200, beta = 0.05, w = 2)
  )) {
    primary <- c(s[["mu"]], s[["beta"]])
    ode <- vapply(
      x, primarycensored_numeric_cdf, numeric(1),
      4L, 1, s[["w"]], 4L, primary
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
  # The CDF is 1e-85 to 1e-225 here, far below the absolute solver
  # tolerance, so the integral is scaled by the largest delay CDF
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
  # The ODE CDF can be 0 or negative by a rounding error of the solver, or
  # above 1, where its log must not be NaN or positive
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

test_that("the series is not used for the gamma where it amplifies the
  shape gradient error", {
  case <- gumbel_stan_cases[[2]]
  beta <- 1
  d <- 0.1
  fit_at <- function(mu) {
    n_terms <- gumbel_n_terms(mu / beta)
    as.vector(primarycensored_gumbel_lcdf_from_terms(
      primarycensored_gumbel_terms(d, 2L, beta, n_terms, case$params),
      primarycensored_gumbel_terms(d - 1, 2L, beta, n_terms, case$params),
      d, 1, mu, beta, n_terms
    ))
  }
  # The value is accurate, but the ratio of the absolute terms to the sum
  # multiplies the error of the derivative of the gamma CDF in the shape
  fit <- fit_at(1.5)
  expect_lt(fit[2], gumbel_error_tolerance())
  expect_gt(
    fit[3] * gumbel_gamma_shape_gradient_error(), gumbel_gradient_tolerance()
  )
  expect_identical(gumbel_series_accepted(2L, fit), 0L)
  expect_identical(gumbel_series_accepted(18L, fit), 1L)
  expect_identical(
    primarycensored_gumbel_lcdf(d, 2L, case$params, 1, 1.5, beta),
    primarycensored_gumbel_numeric_lcdf(d, 2L, case$params, 1, 1.5, beta)
  )
  # A ratio of 32 is accepted
  fit <- fit_at(-8)
  expect_lt(fit[3], 200)
  expect_identical(gumbel_series_accepted(2L, fit), 1L)
  expect_identical(
    primarycensored_gumbel_lcdf(d, 2L, case$params, 1, -8, beta), fit[1]
  )
})

test_that("a series rejected by the estimate is replaced by an accurate ODE", {
  # The estimate of the series here is 1.9e-5, and the ODE is accurate to
  # 1e-8 of the CDF where the series loses precision
  lcdf <- primarycensored_gumbel_lcdf(1, 18L, c(5, 1), 1, 0.3, 0.111)
  expected <- log(gumbel_reference(
    1, 1, 0.3, 0.111, function(x) pnorm(x, 5, 1), FALSE
  ))
  expect_equal(lcdf, expected, tolerance = 1e-8)
  # Where the series was rejected and the ODE used to be off by 23%
  lcdf <- primarycensored_gumbel_lcdf(-2, 18L, c(3, 2), 2, 0.5, 0.2)
  expect_equal(
    exp(lcdf),
    gumbel_reference(-2, 2, 0.5, 0.2, function(x) pnorm(x, 3, 2), FALSE),
    tolerance = 1e-7
  )
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

test_that("check_for_gumbel_vectorized needs an integer pwindow", {
  for (dist_id in c(2L, 4L, 18L)) {
    expect_identical(check_for_gumbel_vectorized(dist_id, 4L, 1), 1L)
    expect_identical(check_for_gumbel_vectorized(dist_id, 4L, 7), 1L)
    expect_identical(check_for_gumbel_vectorized(dist_id, 4L, 1.5), 0L)
    expect_identical(check_for_gumbel_vectorized(dist_id, 4L, 0.5), 0L)
    expect_identical(check_for_gumbel_vectorized(dist_id, 1L, 1), 0L)
  }
  for (dist_id in c(1L, 3L, 26L)) {
    expect_identical(check_for_gumbel_vectorized(dist_id, 4L, 1), 0L)
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
  for (case in gumbel_stan_cases) {
    for (pwindow in c(1, 2, 5)) {
      for (beta in c(0.2, 1)) {
        for (mu in c(-1, 0, 1)) {
          if (!check_for_gumbel_params(case$dist_id, case$params, mu, beta)) {
            next
          }
          for (start in c(1L, 5L)) {
            vectorised <- primarycensored_gumbel_lcdf_vectorized(
              start, n, case$dist_id, case$params, pwindow, mu, beta
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

test_that("primarycensored_lcdf_vectorized uses the Gumbel shared terms", {
  for (case in gumbel_stan_cases) {
    expect_identical(
      primarycensored_lcdf_vectorized(
        1L, 15L, case$dist_id, case$params, 3, 4L, c(-1, 1)
      ),
      primarycensored_gumbel_lcdf_vectorized(
        1L, 15L, case$dist_id, case$params, 3, -1, 1
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
      if (lower_support && setting$L == 0) next
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
            if (lower_support && is.finite(setting$L) && d < setting$L) {
              return(-Inf)
            }
            primarycensored_lpmf(
              d, case$dist_id, case$params, pwindow, d + 1, setting$L,
              setting$D, 4L, c(-1, 1)
            )
          },
          numeric(1)
        )
        if (lower_support) {
          # Delays below zero are not part of the vectorised interval
          next
        }
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

# Gradient of the delay parameters and the primary parameters. The gradient
# of the parameters that an exponential delay does not have is zero, and the
# gamma shape has a gradient with a relative error of about 1e-3 in Stan.
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
    # Points where the series is not accurate and the ODE is used, which
    # has gradients with the accuracy of the solver
    list(d = 4, pwindow = 2, mu = 0.5, beta = 0.1),
    # Narrow spikes of the window density, at the end of the window and
    # inside it, where the ODE used to give gradients of the wrong size
    list(d = 1, pwindow = 1, mu = 2, beta = 0.1),
    list(d = 5, pwindow = 1, mu = 1.5, beta = 0.05),
    list(d = 7, pwindow = 7, mu = 1.5, beta = 0.05),
    list(d = 3, pwindow = 2, mu = 3, beta = 0.06),
    # Around where the series is replaced by the ODE, mu / beta of about
    # 2.5 to 2.7, where the log CDF must be smooth
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

# A model whose target is the log of the standard normal CDF, for the
# gradient of the Stan function primarycensored_log_std_normal_cdf
log_std_normal_gradient_model <- function() {
  testthat::skip_if_not_installed("cmdstanr")
  testthat::skip_if(
    is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))
  )
  functions <- pcd_load_stan_functions(
    wrap_in_block = TRUE, write_to_file = FALSE
  )
  code <- paste0(
    functions, "\n",
    "parameters {\n  real z;\n}\n",
    "model {\n  target += primarycensored_log_std_normal_cdf(z);\n}\n"
  )
  path <- file.path(tempdir(), "pcd_log_std_normal_gradient.stan")
  writeLines(code, path)
  suppressMessages(suppressWarnings(cmdstanr::cmdstan_model(path)))
}

test_that("the Stan log standard normal CDF is accurate with an exact
  gradient in the far lower tail", {
  z <- c(-5, -20, -36.9, -37.1, -45, -60, -100, -300)
  expect_equal(
    vapply(z, primarycensored_log_std_normal_cdf, numeric(1)),
    pnorm(z, log.p = TRUE),
    tolerance = 1e-13
  )
  model <- log_std_normal_gradient_model()
  for (z0 in z) {
    res <- stan_gradient_at(model, data = list(), init = list(z = z0))
    # The derivative is the density over the CDF, which the derivative of
    # std_normal_lcdf() misses by a relative 4e-4 at -60 and 8e-3 at -100
    exact <- exp(dnorm(z0, log = TRUE) - pnorm(z0, log.p = TRUE))
    expect_false(res$gradient_not_finite, info = as.character(z0))
    # CmdStan prints the gradient to 6 significant figures
    expect_equal(
      res$gradient, exact, tolerance = 1e-5, info = as.character(z0)
    )
  }
})

test_that("Gumbel log CDF gradients match finite differences across the
  acceptance grid", {
  # The series path gave the right log CDF but a gradient that was wrong by
  # up to 6 times where the tilts push the normal transform below -37
  model <- gumbel_gradient_model()
  cases <- list(
    list(dist_id = 18L, params = c(0.5, 1), d = c(0.02, 1, 3)),
    list(dist_id = 18L, params = c(3, 2), d = c(0.02, 1, 3, 6)),
    list(dist_id = 4L, params = 60, d = c(0.02, 0.1, 0.5, 2.5)),
    list(dist_id = 2L, params = c(3, 80), d = c(0.02, 0.1, 0.5, 2.5))
  )
  for (case in cases) {
    for (mu in c(-0.5, 0, 0.5, 1, 1.5)) {
      for (beta in c(0.1, 0.2, 1)) {
        for (pwindow in c(1, 2)) {
          for (d in case$d) {
            label <- gumbel_case_label(
              case, d = d, pwindow = pwindow, mu = mu, beta = beta
            )
            # A CDF below 1e-260 is not a point to fit, finite differences
            # of it are not meaningful
            lcdf <- primarycensored_lcdf(
              d, case$dist_id, case$params, pwindow,
              gumbel_case_lower(case), Inf, 4L, c(mu, beta)
            )
            if (lcdf < -600) next
            res <- gumbel_gradient_at(model, case, d, pwindow, mu, beta)
            expect_false(res$gradient_not_finite, info = label)
            expect_false(res$rejected, info = label)
            expect_length(res$gradient, 4)
            # The finite differences of the log CDF have an absolute error of
            # about 3e-5, from rounding in the log CDF of 3e-11
            expect_gumbel_gradient_close(res, case, label, slack = 1e-4)
          }
        }
      }
    }
  }
})

test_that("Gumbel gradients are finite for delays whose CDF has an
  infinite slope at 0", {
  # A gamma or Weibull with shape below 1 has an unbounded density at 0. The
  # series is not admissible for these in daily units, so mu below pwindow
  # with d at or below pwindow takes the numerical z branch, where the upper
  # limit is d. Its sensitivity used to be integrated and is singular, which
  # ended in the solver step limit after about 15 seconds
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
      elapsed <- system.time(
        res <- gumbel_gradient_at(
          model, case, point$d, point$pwindow, point$mu, point$beta
        )
      )[["elapsed"]]
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
      elapsed <- system.time(
        res <- gumbel_gradient_at(
          model, case, point$d, point$pwindow, point$mu, point$beta,
          vectorised = TRUE
        )
      )[["elapsed"]]
      expect_false(res$gradient_not_finite, info = label)
      expect_false(res$rejected, info = label)
      expect_true(all(is.finite(res$gradient)), info = label)
      expect_lt(elapsed, 5, label = label)
      expect_gumbel_gradient_close(res, case, label, scale = 5)
    }
  }
})

test_that("the Stan numerical path keeps the lower tail on the log scale", {
  # The kink of the delay CDF beyond the integration range gave -Inf, the
  # density over the integration range below the solver tolerance gave -Inf
  # or a wrong value, and the mass beyond a window quantile of 1 - 2^-53 was
  # cut. The references are an independent log scale integral, and R
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
  # The PMF of a delay as short as the exponential is the difference of CDFs
  # that are 1 to double precision, so only the normal delays are used. The
  # gradient of the upper tail log PMF is wrong where the log CDF is close to
  # 0, see #404. The finite differences of a sum of PMFs in the upper tail
  # have an error of about 5e-4, larger than the 1e-4 of a single log CDF.
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

test_that("the model rejects the reserved primary identifier 3", {
  testthat::skip_if_not_installed("cmdstanr")
  testthat::skip_if_not_installed("dplyr")
  testthat::skip_if(
    is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))
  )
  delay_counts <- data.frame(
    delay = 1:3, delay_upper = 2:4, n = 5L, pwindow = 1L,
    start_relative_obs_time = 0, relative_obs_time = 10
  )
  stan_data <- pcd_as_stan_data(
    delay_counts,
    dist_id = pcd_stan_dist_id("lognormal", "delay"),
    primary_id = 1,
    param_bounds = list(lower = c(-Inf, 0.01), upper = c(Inf, Inf)),
    primary_param_bounds = list(lower = numeric(0), upper = numeric(0)),
    priors = list(location = c(0, 1), scale = c(5, 2.5)),
    primary_priors = list(location = numeric(0), scale = numeric(0))
  )
  model <- gumbel_pcens_model()
  # The identifier 3 is reserved and has no primary distribution
  stan_data$primary_id <- 3L
  messages <- character()
  withCallingHandlers(
    suppressWarnings(tryCatch(
      model$sample(
        data = stan_data, chains = 1, iter_warmup = 1, iter_sampling = 1,
        refresh = 0, show_messages = TRUE
      ),
      error = function(e) NULL
    )),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m)) # nolint
      invokeRestart("muffleMessage")
    }
  )
  expect_true(
    any(grepl("primary_id 3 is reserved", messages, fixed = TRUE))
  )
})

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
