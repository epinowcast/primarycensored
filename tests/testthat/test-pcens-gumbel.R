families <- gumbel_families()
mus <- c(-0.5, 0, 0.5, 1, 1.5)
betas <- c(0.1, 0.2, 1)
pwindows <- c(1, 2)
# Relative tolerance of the series where it is accepted, and absolute
# tolerance of the public CDF everywhere
series_tol <- 1e-8
public_tol <- 1e-5

test_that("normal delays dispatch to the series method", {
  expect_s3_class(
    gumbel_object(families[[4]], 0, 1), "pcens_pnorm_dtgumbel"
  )
  expect_false(is.null(
    utils::getS3method("pcens_cdf", "pcens_pnorm_dtgumbel", optional = TRUE)
  ))
  for (cls in c("pcens_pexp_dtgumbel", "pcens_pgamma_dtgumbel")) {
    expect_null(utils::getS3method("pcens_cdf", cls, optional = TRUE))
  }
})

test_that("the number of terms meets the truncation bound", {
  expect_gte(.gumbel_n_terms(-10), 2)
  for (log_s0 in c(-5, -1, 0, 1, 2.5)) {
    n <- .gumbel_n_terms(log_s0)
    bound <- n * log_s0 - lgamma(n + 1) - max(log_s0, 0)
    expect_lt(bound, log(1e-18))
    bound_before <- (n - 1) * log_s0 - lgamma(n) - max(log_s0, 0)
    expect_gte(bound_before, log(1e-18))
  }
})

test_that("the series agrees with numerical integration where accepted", {
  for (family in gumbel_series_families()) {
    cdf <- gumbel_cdf(family)
    for (pwindow in pwindows) {
      for (beta in betas) {
        for (mu in mus) {
          label <- gumbel_label(family, pwindow, mu, beta)
          obj <- gumbel_object(family, mu, beta)
          if (!.gumbel_available(obj, mu, beta)) next
          lower <- .pcens_tilt_lower(obj)
          q <- family$q[!is.finite(lower) | family$q > lower]
          n_terms <- .gumbel_n_terms(mu / beta)
          fit <- .gumbel_lcdf(obj, q, pwindow, mu, beta, n_terms, lower)
          reference <- gumbel_reference(
            q, pwindow, mu, beta, cdf, family$positive
          )
          accepted <- fit$error <= .gumbel_tol
          if (!any(accepted)) next
          expect_lt(
            max_rel_diff(
              exp(fit$log_cdf)[accepted], reference[accepted]
            ),
            series_tol,
            label = label
          )
          # The estimate bounds the actual error
          actual <- abs(exp(fit$log_cdf)[accepted] / reference[accepted] - 1)
          expect_true(
            all(actual <= 10 * fit$error[accepted] + 1e-12),
            label = paste("estimate bounds error,", label)
          )
        }
      }
    }
  }
})

test_that("the public CDF agrees with numerical integration everywhere", {
  for (family in families) {
    cdf <- gumbel_cdf(family)
    for (pwindow in pwindows) {
      for (beta in betas) {
        for (mu in mus) {
          label <- gumbel_label(family, pwindow, mu, beta)
          obj <- gumbel_object(family, mu, beta)
          expected <- gumbel_reference(
            family$q, pwindow, mu, beta, cdf, family$positive
          )
          expect_lt(
            max(abs(pcens_cdf(obj, family$q, pwindow) - expected)),
            public_tol,
            label = label
          )
        }
      }
    }
  }
})

test_that("the series is accurate at extreme and narrow windows", {
  # Very negative mu / beta, a wide scale and a narrow window
  family <- families[[4]]
  cdf <- gumbel_cdf(family)
  settings <- list(
    c(mu = -5, beta = 0.1, w = 1), c(mu = -2, beta = 0.1, w = 2),
    c(mu = -50, beta = 1, w = 1), c(mu = 0, beta = 1, w = 0.05),
    c(mu = 0, beta = 5, w = 0.5), c(mu = 3, beta = 5, w = 2),
    c(mu = 0.3, beta = 10, w = 1), c(mu = -1, beta = 0.3, w = 0.1)
  )
  for (s in settings) {
    obj <- gumbel_object(family, s[["mu"]], s[["beta"]])
    expect_true(.gumbel_available(obj, s[["mu"]], s[["beta"]]))
    fit <- .gumbel_lcdf(
      obj, family$q, s[["w"]], s[["mu"]], s[["beta"]],
      .gumbel_n_terms(s[["mu"]] / s[["beta"]]), -Inf
    )
    expect_true(all(fit$error <= .gumbel_tol), info = toString(s))
    expect_lt(
      max_rel_diff(
        exp(fit$log_cdf),
        gumbel_reference(
          family$q, s[["w"]], s[["mu"]], s[["beta"]], cdf, FALSE
        )
      ),
      1e-10,
      label = toString(s)
    )
  }
})

test_that("use_numeric returns the numerical result", {
  obj <- gumbel_object(families[[4]], -0.5, 0.5)
  q <- c(-1, 3, 8)
  expect_identical(
    pcens_cdf(obj, q, 1, use_numeric = TRUE),
    pcens_cdf.default(obj, q, 1)
  )
})

test_that("a large window value falls back to the numerical path", {
  # exp(mu / beta) is about 148 here
  obj <- gumbel_object(families[[4]], 0.5, 0.1)
  expect_false(.gumbel_available(obj, 0.5, 0.1))
  q <- c(-1, 3, 8)
  expect_identical(pcens_cdf(obj, q, 1), pcens_cdf.default(obj, q, 1))
})

test_that("exponential and gamma delays use the numerical path", {
  q <- c(0.5, 3, 8)
  for (family in families[1:2]) {
    obj <- gumbel_object(family, 0, 1)
    expect_identical(pcens_cdf(obj, q, 2), pcens_cdf.default(obj, q, 2))
  }
})

test_that("the series is skipped where its rounding estimate is large", {
  # A wide normal delay with a large mu / beta loses the series precision
  wide <- new_pcens(
    pnorm, dtgumbel, primary_args = list(mu = 1.15, beta = 0.5),
    mean = 0, sd = 4
  )
  expect_false(.gumbel_available(wide, 1.15, 0.5))
  narrow <- gumbel_object(families[[4]], 0.5, 0.5)
  expect_true(.gumbel_available(narrow, 0.5, 0.5))
  q <- c(-5, 0, 3, 8)
  expect_equal(
    pcens_cdf(wide, q, 1),
    gumbel_reference(q, 1, 1.15, 0.5, function(x) pnorm(x, 0, 4), FALSE),
    tolerance = 1e-7
  )
})

test_that("a narrow window relative to the scale uses the numerical path", {
  # The bracket is a small difference of large terms when w << beta
  obj <- gumbel_object(families[[4]], 0, 1)
  fit <- .gumbel_lcdf(
    obj, c(1, 3), 1e-7, 0, 1, .gumbel_n_terms(0), -Inf
  )
  expect_true(all(fit$error > .gumbel_tol))
  q <- c(1, 3)
  expect_equal(
    pcens_cdf(obj, q, 1e-7),
    gumbel_reference(q, 1e-7, 0, 1, gumbel_cdf(families[[4]]), FALSE),
    tolerance = 1e-6
  )
})

test_that("q at and below the support and at infinity are handled", {
  for (family in families[c(1, 4)]) {
    obj <- gumbel_object(family, -0.5, 1)
    result <- pcens_cdf(obj, c(-Inf, Inf, NA), 1)
    expect_identical(result[1:2], c(0, 1))
    expect_true(is.na(result[3]))
  }
  exp_obj <- gumbel_object(families[[1]], -0.5, 1)
  expect_identical(pcens_cdf(exp_obj, c(-3, 0), 1), c(0, 0))
})

test_that("small q for delays on the non-negative reals is accurate", {
  family <- families[[1]]
  cdf <- gumbel_cdf(family)
  q <- c(1e-6, 1e-4, 0.01, 0.5, 0.99, 1.01)
  for (mu in c(-1, 0.5)) {
    obj <- gumbel_object(family, mu, 1)
    expect_equal(
      pcens_cdf(obj, q, 1),
      gumbel_reference(q, 1, mu, 1, cdf, TRUE),
      tolerance = 1e-8
    )
  }
})

test_that("tail probabilities are accurate on the relative scale", {
  family <- families[[4]]
  cdf <- gumbel_cdf(family)
  obj <- gumbel_object(family, -0.5, 1)
  # Far into the lower tail, where the CDF is small
  q <- c(-12, -9, -6)
  expected <- gumbel_reference(q, 2, -0.5, 1, cdf, FALSE)
  expect_true(all(expected < 1e-3 & expected > 0))
  expect_lt(max_rel_diff(pcens_cdf(obj, q, 2), expected), 1e-8)
  # Far into the upper tail the complement is accurate
  q <- c(9, 12, 15)
  expected <- gumbel_reference(q, 2, -0.5, 1, cdf, FALSE)
  expect_equal(1 - pcens_cdf(obj, q, 2), 1 - expected, tolerance = 1e-10)
})

test_that("pwindow must be a single positive window for the series", {
  obj <- gumbel_object(families[[4]], -0.5, 1)
  q <- c(0, 3)
  expect_equal(
    pcens_cdf(obj, q, 1),
    gumbel_reference(q, 1, -0.5, 1, gumbel_cdf(families[[4]]), FALSE),
    tolerance = 1e-8
  )
  expect_identical(pcens_cdf(obj, q, 0), pnorm(q, 3, 2))
})

test_that("the primary arguments are checked", {
  obj <- new_pcens(pnorm, dtgumbel, primary_args = list(mu = 0), mean = 3)
  expect_error(pcens_cdf(obj, 1, 1), "mu and beta")
  obj <- new_pcens(
    pnorm, dtgumbel,
    primary_args = list(mu = 0, beta = -1), mean = 3
  )
  expect_error(pcens_cdf(obj, 1, 1), "beta")
  obj <- new_pcens(
    pnorm, dtgumbel,
    primary_args = list(mu = c(0, 1), beta = 1), mean = 3
  )
  expect_error(pcens_cdf(obj, 1, 1), "single")
})

test_that("the PMF shares endpoints and matches the reference", {
  family <- families[[4]]
  cdf <- gumbel_cdf(family)
  for (pwindow in pwindows) {
    obj <- gumbel_object(family, 0, 0.5)
    x <- -3:10
    expected <- diff(gumbel_reference(-3:11, pwindow, 0, 0.5, cdf, FALSE))
    expect_equal(
      pcens_pmf(obj, x, pwindow), expected,
      tolerance = 1e-7
    )
    expect_equal(sum(pcens_pmf(obj, -40:40, pwindow)), 1, tolerance = 1e-9)
  }
  exp_family <- families[[1]]
  obj <- gumbel_object(exp_family, -1, 1)
  x <- c(0, 0.05, 0.3)
  pmf <- pcens_pmf(obj, x, 2, swindow = 0.1)
  expected <- gumbel_reference(
    x + 0.1, 2, -1, 1, gumbel_cdf(exp_family), TRUE
  ) - gumbel_reference(x, 2, -1, 1, gumbel_cdf(exp_family), TRUE)
  expect_equal(pmf, expected, tolerance = 1e-7)
})

test_that("a narrow spike of the window density is resolved", {
  # The window density is a spike when mu / beta is large
  for (family in gumbel_spike_families()) {
    cdf <- gumbel_cdf(family)
    for (s in gumbel_spike_settings()) {
      label <- gumbel_label(family, s[["w"]], s[["mu"]], s[["beta"]])
      obj <- gumbel_object(family, s[["mu"]], s[["beta"]])
      expected <- gumbel_reference(
        family$q, s[["w"]], s[["mu"]], s[["beta"]], cdf, family$positive
      )
      expect_lt(
        gumbel_error(pcens_cdf(obj, family$q, s[["w"]]), expected),
        1e-7,
        label = label
      )
      expect_lt(
        gumbel_error(
          pcens_cdf(obj, family$q, s[["w"]], use_numeric = TRUE), expected
        ),
        1e-7,
        label = paste("use_numeric,", label)
      )
    }
  }
})

test_that("the numerical path is accurate for large mu over beta", {
  family <- gumbel_spike_families()[[1]]
  cdf <- gumbel_cdf(family)
  for (ratio in c(15, 20, 50, 200)) {
    beta <- 0.1
    mu <- ratio * beta
    obj <- gumbel_object(family, mu, beta)
    # Above, inside and below the window
    for (pwindow in c(0.4, 1, 3)) {
      expected <- gumbel_reference(
        family$q, pwindow, mu, beta, cdf, FALSE
      )
      expect_lt(
        gumbel_error(pcens_cdf(obj, family$q, pwindow), expected),
        1e-7,
        label = sprintf("mu over beta %g, pwindow %g", ratio, pwindow)
      )
    }
  }
})

test_that("a point mass at the window end is the limit of a large mu", {
  # mu / beta of 1000 puts the mass within 1e-300 of pwindow
  obj <- gumbel_object(gumbel_spike_families()[[1]], 100, 0.1)
  q <- c(-3, 1, 5)
  expect_equal(pcens_cdf(obj, q, 1), pnorm(q - 1, 3, 2), tolerance = 1e-8)
})

test_that("dtgumbel passes the density check for a narrow spike", {
  expect_null(
    check_dprimary(dtgumbel, pwindow = 1, list(mu = 2, beta = 0.1))
  )
  expect_error(
    check_dprimary(dtgumbel, pwindow = 1, list(mu = 2, beta = -1)),
    "beta"
  )
})

test_that("pprimarycensored works for a narrow spike", {
  p <- pprimarycensored(
    c(1, 5, 20), pnorm,
    pwindow = 1, dprimary = dtgumbel,
    primary_args = list(mu = 2, beta = 0.1), mean = 3, sd = 2
  )
  expect_equal(
    p, c(0.0668075, 0.6914633, 1), tolerance = 1e-6
  )
})

test_that("delays without a solution use the numerical path for a spike", {
  # The lognormal has no solution, so it uses pcens_cdf.default()
  cdf <- function(x) plnorm(x, 1, 0.5)
  q <- c(0.5, 1, 3, 8, 20)
  for (s in gumbel_spike_settings()) {
    obj <- new_pcens(
      plnorm, dtgumbel,
      primary_args = list(mu = s[["mu"]], beta = s[["beta"]]),
      meanlog = 1, sdlog = 0.5
    )
    expected <- gumbel_reference(
      q, s[["w"]], s[["mu"]], s[["beta"]], cdf, TRUE
    )
    expect_lt(
      gumbel_error(pcens_cdf(obj, q, s[["w"]]), expected), 1e-7,
      label = toString(s)
    )
  }
})

test_that("the numerical path keeps a small delay at the end of the window", {
  # A spike within rounding of the window end at q = pwindow
  obj <- gumbel_object(gumbel_spike_families()[[3]], 2.7, 0.02)
  cdf <- gumbel_cdf(gumbel_spike_families()[[3]])
  expected <- gumbel_reference(1, 1, 2.7, 0.02, cdf, TRUE)
  expect_gt(expected, 0)
  expect_equal(pcens_cdf(obj, 1, 1), expected, tolerance = 1e-7)
})

test_that("the numerical path agrees with the reference where the delay
  starts well after a spike at the window end", {
  for (case in gumbel_late_cases()) {
    cdf <- gumbel_case_cdf(case)
    for (mu in c(1.5, 2, 3)) {
      for (beta in c(0.05, 0.1)) {
        obj <- gumbel_object(case, mu, beta)
        actual <- pcens_cdf(obj, c(1, 2, 3), 1)
        expected <- gumbel_reference(c(1, 2, 3), 1, mu, beta, cdf, TRUE)
        expect_lt(
          max(abs(actual - expected)), 1e-9,
          label = gumbel_case_label(case, mu = mu, beta = beta)
        )
        numeric <- pcens_cdf(obj, c(1, 2, 3), 1, use_numeric = TRUE)
        expect_equal(actual, numeric, tolerance = 1e-12)
      }
    }
  }
})

test_that("rprimarycensored samples match the analytic CDF", {
  args <- list(
    list(rdist = rexp, pdist = pexp, rate = 0.4),
    list(rdist = rgamma, pdist = pgamma, shape = 3, rate = 1),
    list(rdist = rnorm, pdist = pnorm, mean = 5, sd = 1),
    list(rdist = rlnorm, pdist = plnorm, meanlog = 1, sdlog = 0.5),
    list(rdist = rweibull, pdist = pweibull, shape = 2, scale = 3),
    list(rdist = rlogis, pdist = plogis, location = 2, scale = 1)
  )
  # Location before, inside and after the window
  primaries <- list(
    list(mu = -1, beta = 0.5),
    list(mu = 0.5, beta = 0.3),
    list(mu = 3, beta = 1)
  )
  for (primary_args in primaries) {
    for (a in args) {
      set.seed(10)
      dist_args <- a[setdiff(names(a), c("rdist", "pdist"))]
      draws <- do.call(
        rprimarycensored,
        c(
          list(
            n = 2000, rdist = a$rdist, pwindow = 2, swindow = 0,
            rprimary = rtgumbel, rprimary_args = primary_args
          ),
          dist_args
        )
      )
      ks <- suppressWarnings(stats::ks.test(
        draws,
        function(x) {
          do.call(
            pprimarycensored,
            c(
              list(
                q = x, pdist = a$pdist, pwindow = 2, dprimary = dtgumbel,
                primary_args = primary_args
              ),
              dist_args
            )
          )
        }
      ))
      expect_gt(
        ks$p.value, 0.001,
        label = paste(deparse(a$rdist), "mu", primary_args$mu)
      )
    }
  }
})
