# Targeted coverage of the non-parametric internal helpers in
# R/nonparametric_helpers.R that are reached only on error/edge paths.

test_that(".is_valid_simplex rejects non-numeric and NA-bearing inputs", {
  for (bad in list("a", NA_real_, c(0.5, NA_real_), list(0.5, 0.5))) {
    expect_false(
      primarycensored:::.is_valid_simplex(bad),
      info = paste("rejected:", paste(class(bad), collapse = ","))
    )
  }
  # Negative values and bad sums also fail.
  expect_false(primarycensored:::.is_valid_simplex(c(-0.1, 1.1)))
  expect_false(primarycensored:::.is_valid_simplex(c(0.5, 0.6)))
  # A genuine simplex passes.
  expect_true(primarycensored:::.is_valid_simplex(c(0.5, 0.5)))
})

test_that(".sort_eps_names preserves order when suffixes are non-integer", {
  # Lexicographic suffixes 10 < 2 would scramble; integer suffixes work.
  expect_identical(
    primarycensored:::.sort_eps_names(c("eps_10", "eps_2", "eps_1")),
    c("eps_1", "eps_2", "eps_10")
  )
  # Non-integer suffix: helper returns the input unchanged.
  expect_identical(
    primarycensored:::.sort_eps_names(c("eps_a", "eps_b")),
    c("eps_a", "eps_b")
  )
  # Length-0/1 short-circuit.
  expect_identical(primarycensored:::.sort_eps_names(character(0)), character(0))
  expect_identical(primarycensored:::.sort_eps_names("eps_1"), "eps_1")
})

test_that(".make_hazard_transform handles both models and K = 1 (no eps)", {
  rw <- primarycensored:::.make_hazard_transform("rw")
  re <- primarycensored:::.make_hazard_transform("re")
  # K = 1: no eps_*; final hazard pinned to 1.
  expect_equal(rw(list(alpha = -1, log_sigma = log(0.5))), 1, tolerance = 0)
  expect_equal(re(list(alpha = -1, log_sigma = log(0.5))), 1, tolerance = 0)
  # K = 4: RW cumulates eps, RE doesn't. Need K >= 4 for the
  # difference to surface in h[1..K-1] (h[K] is pinned to 1).
  par <- list(
    alpha = -1, log_sigma = log(0.5),
    eps_1 = 0.5, eps_2 = -0.3, eps_3 = 0.4
  )
  h_rw <- rw(par)
  h_re <- re(par)
  expect_length(h_rw, 4L)
  expect_length(h_re, 4L)
  expect_equivalent_numeric <- function(a, b) {
    expect_equal(unname(a), b, tolerance = 1e-12)
  }
  expect_equivalent_numeric(h_rw[4], 1)
  expect_equivalent_numeric(h_re[4], 1)
  expect_false(isTRUE(all.equal(h_rw, h_re)))
})

test_that(".hazard_fit_penalty returns finite penalty amortised across N", {
  par <- list(
    alpha = -1, log_sigma = log(0.5),
    eps_1 = 0.1, eps_2 = -0.2, eps_3 = 0.3
  )
  nlp_per_obs <- primarycensored:::.hazard_fit_penalty(par, N = 100)
  expect_true(is.finite(nlp_per_obs))
  expect_gt(nlp_per_obs, 0)
  # Doubling N halves the per-observation penalty.
  nlp_double <- primarycensored:::.hazard_fit_penalty(par, N = 200)
  expect_equal(nlp_double, nlp_per_obs / 2, tolerance = 1e-12)
  # User priors flow through and shift the penalty.
  shifted <- primarycensored:::.hazard_fit_penalty(
    par, N = 100,
    prior_settings = list(alpha = list(mean = -1), log_sigma = list(sd = 5))
  )
  expect_true(is.finite(shifted))
  expect_false(isTRUE(all.equal(shifted, nlp_per_obs)))
  # K = 1: no eps_* parameters, only alpha/log_sigma contribute.
  nlp_no_eps <- primarycensored:::.hazard_fit_penalty(
    list(alpha = -1, log_sigma = log(0.5)), N = 100
  )
  expect_true(is.finite(nlp_no_eps))
})
