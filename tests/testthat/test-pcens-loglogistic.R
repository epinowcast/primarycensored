# Log-logistic delays with a uniform primary. The CDF is
# (G_1(q) - G_1(q - w)) / w with G_1(t) = int_0^t F(u) du, which needs only
# the partial moment ratio r_{1/shape}.

families <- loglogistic_families()

test_that("log-logistic delays dispatch to the uniform primary method", {
  obj <- uniform_object(families[[2]])
  expect_s3_class(obj, "pcens_pllogis_dunif")
  expect_false(
    is.null(utils::getS3method(
      "pcens_cdf", "pcens_pllogis_dunif",
      optional = TRUE
    ))
  )
})

test_that("the partial moment ratio matches numerical integration", {
  a <- c(0, 0.25, 0.5, 1, 1.0001, 2, 2.5, 3, 7.3, 20, 50)
  A <- c(1e-8, 1e-3, 0.3, 0.9999, 1.0001, 1.5, 2.9, 3.0001, 5, 30, 1e3, 1e5)
  ratios <- .loglogistic_ratio(a, log(A))
  expect_identical(dim(ratios), c(length(A), length(a)))
  for (i in seq_along(A)) {
    for (j in seq_along(a)) {
      expect_equal(
        ratios[i, j], loglogistic_ratio_reference(a[j], A[i]),
        tolerance = 1e-10,
        info = sprintf("a = %g, A = %g", a[j], A[i])
      )
    }
  }
})

test_that("the partial moment ratio handles very large and small A", {
  # r_a(A) for A far into the upper tail follows M_a(Inf) / A^a for a < 1
  # and the logarithmic forms at a = 1
  ratios <- .loglogistic_ratio(c(0.5, 1, 2), log(c(1e10, 1e30)))
  expect_true(all(is.finite(ratios)))
  expect_equal(ratios[1, 2], (log1p(1e10) - 1 + 1 / (1 + 1e10)) *
    (1 + 1e10) / 1e10^2, tolerance = 1e-10)
  expect_equal(ratios[1, 1], 0.5 * pi / 1e5 * (1 + 1e-10), tolerance = 1e-3)
  # A beyond the range of a double is handled on the log scale
  huge <- .loglogistic_ratio(c(0.5, 2), 800)
  expect_true(all(is.finite(huge)))
  # r_2 is below the smallest double
  expect_gt(huge[1], 0)
  expect_gte(huge[2], 0)
  tiny <- .loglogistic_ratio(c(0.5, 2), -800)
  expect_equal(as.vector(tiny), c(1 / 1.5, 1 / 3), tolerance = 1e-14)
})

# Reference (1 / w) int_{max(q - w, 0)}^{q} F(t) dt
uniform_reference <- function(q, pwindow, cdf) {
  vapply(q, function(qq) {
    lower <- max(qq - pwindow, 0)
    if (qq <= 0) {
      return(0)
    }
    stats::integrate(
      cdf, lower, qq,
      rel.tol = 1e-13, abs.tol = 0, subdivisions = 2000L
    )$value / pwindow
  }, numeric(1))
}

test_that("the analytic CDF matches a reference integral", {
  for (family in families) {
    cdf <- exptilt_cdf(family)
    obj <- uniform_object(family)
    for (pwindow in c(0.5, 1, 2, 7)) {
      q <- sort(c(
        1e-6, 1e-3, 0.3 * pwindow, pwindow - 1e-3, pwindow, pwindow + 1e-3,
        2, 3, 6, 12, 25, 100, 1e4
      ))
      expected <- uniform_reference(q, pwindow, cdf)
      actual <- pcens_cdf(obj, q, pwindow)
      expect_lt(
        max_rel_diff(actual, expected), 1e-8,
        label = paste(family$label, "pwindow", pwindow)
      )
    }
  }
})

test_that("the analytic CDF agrees with use_numeric = TRUE", {
  q <- c(0.05, 0.5, 1.5, 3, 6, 12, 20)
  for (family in families) {
    obj <- uniform_object(family)
    for (pwindow in c(1, 2, 7)) {
      expect_equal(
        pcens_cdf(obj, q, pwindow),
        pcens_cdf(obj, q, pwindow, use_numeric = TRUE),
        tolerance = 1e-6,
        info = paste(family$label, "pwindow", pwindow)
      )
    }
  }
})

test_that("the analytic CDF handles boundary values and the pmf", {
  obj <- uniform_object(families[[3]])
  expect_identical(pcens_cdf(obj, c(-Inf, -1, 0), 2), c(0, 0, 0))
  expect_identical(pcens_cdf(obj, Inf, 2), 1)
  expect_identical(pcens_cdf(obj, numeric(0), 2), numeric(0))
  expect_true(is.na(pcens_cdf(obj, NA_real_, 2)))
  cdf <- exptilt_cdf(families[[3]])
  expected <- diff(uniform_reference(0:12, 3, cdf))
  expect_equal(pcens_pmf(obj, 0:11, pwindow = 3), expected, tolerance = 1e-9)
})

test_that("the CDF is accurate where the difference cancels", {
  for (case in uniform_conditioning_cases) {
    obj <- new_pcens(
      pdist = pllogis_test, dprimary = dunif, primary_args = list(),
      shape = case[1], scale = case[2]
    )
    ref <- loglogistic_censored_reference(case[4], case[1], case[2], 0, case[3])
    actual <- pcens_cdf(obj, case[4], case[3])
    # The survival is resolved to 1e-15 in absolute terms
    expect_lt(
      abs(actual - ref[["cdf"]]) / max(min(ref), 1e-9), 1e-6,
      label = paste(case, collapse = " ")
    )
  }
})

test_that("the uniform PMF is positive and accurate at large q", {
  for (case in uniform_pmf_cases) {
    obj <- new_pcens(
      pdist = pllogis_test, dprimary = dunif, primary_args = list(),
      shape = case[1], scale = case[2]
    )
    ref <- vapply(
      case[4] + 0:1, loglogistic_censored_reference, numeric(2),
      shape = case[1], scale = case[2], rho = 0, pwindow = case[3]
    )
    expected <- -diff(ref["survival", ])
    actual <- pcens_pmf(obj, case[4], pwindow = case[3])
    expect_gt(actual, 0)
    expect_lt(
      abs(actual - expected) / expected, 1e-6,
      label = paste(case, collapse = " ")
    )
  }
})

test_that("the uniform CDF keeps the series where it is well conditioned", {
  obj <- uniform_object(families[[3]])
  # Shape 2, scale 5 and q / w of 10 is well within the tolerance
  ill <- function(q) {
    .loglogistic_uniform_ill(.loglogistic_uniform_terms(obj, q, 1))
  }
  expect_false(ill(10))
  expect_true(ill(1e6))
})

test_that("a shape below 0.01 uses the numerical method", {
  tiny <- list(pdist = pllogis_test, args = list(shape = 0.005, scale = 2))
  expect_identical(
    pcens_cdf(uniform_object(tiny), c(1, 4), 2),
    pcens_cdf(uniform_object(tiny), c(1, 4), 2, use_numeric = TRUE)
  )
})

test_that("use_numeric = TRUE uses the default method", {
  obj <- uniform_object(families[[2]])
  expect_identical(
    pcens_cdf(obj, c(1, 4), 2, use_numeric = TRUE),
    pcens_cdf.default(obj, c(1, 4), 2)
  )
})

test_that("rate and scale parametrisations and defaults are accepted", {
  base <- new_pcens(pllogis_test, dunif, list(), shape = 2, scale = 4)
  by_rate <- new_pcens(
    add_name_attribute(
      function(q, shape, rate) pllogis_test(q, shape, 1 / rate), "pllogis"
    ),
    dunif, list(),
    shape = 2, rate = 0.25
  )
  expect_equal(
    pcens_cdf(by_rate, c(1, 5, 11), 2), pcens_cdf(base, c(1, 5, 11), 2),
    tolerance = 1e-14
  )
  default_scale <- new_pcens(pllogis_test, dunif, list(), shape = 2)
  expect_equal(
    pcens_cdf(default_scale, c(0.5, 2), 1),
    pcens_cdf(default_scale, c(0.5, 2), 1, use_numeric = TRUE),
    tolerance = 1e-6
  )
  expect_error(
    pcens_cdf(new_pcens(pllogis_test, dunif, list(), scale = 2), 1, 2),
    "shape parameter is required for the log-logistic"
  )
})

test_that("flexsurv and actuar log-logistic functions are recognised", {
  skip_if_not_installed("flexsurv")
  skip_if_not_installed("actuar")
  expected <- pcens_cdf(
    new_pcens(pllogis_test, dunif, list(), shape = 2, scale = 4),
    c(1, 5, 11), 2
  )
  for (pdist in list(flexsurv::pllogis, actuar::pllogis)) {
    obj <- new_pcens(pdist, dunif, list(), shape = 2, scale = 4)
    expect_s3_class(obj, "pcens_pllogis_dunif")
    expect_equal(pcens_cdf(obj, c(1, 5, 11), 2), expected, tolerance = 1e-14)
  }
})

test_that("the analytic CDF matches rprimarycensored samples", {
  withr::local_seed(2)
  n <- 20000
  pwindow <- 2
  for (family in families[c(1, 3, 5)]) {
    samples <- loglogistic_samples(n, family, pwindow, 0)
    qs <- unname(stats::quantile(samples, c(0.05, 0.25, 0.5, 0.75, 0.95)))
    expect_lt(
      max(abs(
        vapply(qs, function(q) mean(samples <= q), numeric(1)) -
          pcens_cdf(uniform_object(family), qs, pwindow)
      )),
      0.015,
      label = family$label
    )
  }
})
