skip_on_cran()

# Simple step PMF: 3 intervals, boundaries 0:3
pmf_vals <- c(0.2, 0.5, 0.3)
boundaries_vals <- 0:3

# Packed params for dist_id 26: [boundaries (K+1=4), pmf (K=3)],
# total length = 7
packed_params <- c(as.numeric(boundaries_vals), pmf_vals)

test_that("Stan pstep_lcdf matches R step CDF on a grid", {
  cum <- cumsum(pmf_vals)
  r_step_cdf <- function(t) {
    if (t < boundaries_vals[1]) {
      return(0)
    }
    if (t >= boundaries_vals[length(boundaries_vals)]) {
      return(1)
    }
    k <- max(which(boundaries_vals[seq_along(pmf_vals)] <= t))
    cum[k]
  }
  grid <- c(-0.5, 0, 0.5, 1, 1.5, 2, 2.5, 2.999, 3, 3.5)
  for (t in grid) {
    r_val <- r_step_cdf(t)
    r_log <- if (r_val == 0) -Inf else log(r_val)
    stan_val <- pstep_lcdf(t, as.numeric(boundaries_vals), pmf_vals)
    expect_equal(stan_val, r_log,
      tolerance = 1e-9,
      info = sprintf("pstep_lcdf mismatch at t = %s", t)
    )
  }
})

test_that("Stan dist_lcdf with dist_id 26 matches pstep_lcdf directly", {
  for (t in c(0.5, 1.5, 2.5, 2.99)) {
    via_dist <- dist_lcdf(t, packed_params, 26L)
    direct <- pstep_lcdf(t, as.numeric(boundaries_vals), pmf_vals)
    expect_equal(via_dist, direct, tolerance = 1e-9)
  }
})

test_that("Stan hazards_to_pmf matches R reference", {
  hazards <- c(0.2, 0.4, 1.0)
  r_pmf <- numeric(length(hazards))
  surv <- 1.0
  for (i in seq_along(hazards)) {
    r_pmf[i] <- hazards[i] * surv
    surv <- surv * (1 - hazards[i])
  }
  stan_pmf <- hazards_to_pmf(hazards)
  expect_equal(stan_pmf, r_pmf, tolerance = 1e-9)
  expect_equal(sum(stan_pmf), 1, tolerance = 1e-9)
})

test_that("Stan analytic step CDF matches R pprimarycensored", {
  # Cross-check the analytic Stan path against the R reference
  # `pprimarycensored` for both supported primaries. The closed-form
  # algebra is also exercised by the R-side `pdiscretestep` tests.
  d_vals <- c(0.5, 1.0, 1.5, 2.0, 2.5, 2.99)

  # Uniform primary, pwindow = 1
  for (d in d_vals) {
    stan_val <- primarycensored_cdf(
      d, 26L, packed_params, 1.0, 0, Inf, 1L, numeric(0)
    )
    r_ref <- pprimarycensored(
      d, pdiscretestep,
      pwindow = 1.0,
      boundaries = boundaries_vals, pmf = pmf_vals
    )
    expect_equal(stan_val, r_ref,
      tolerance = 1e-6,
      info = sprintf("uniform mismatch at d = %s", d)
    )
  }

  # Exponential-growth primary, pwindow = 3, r = 0.3
  for (d in d_vals) {
    stan_val <- primarycensored_cdf(
      d, 26L, packed_params, 3.0, 0, Inf, 2L, c(0.3)
    )
    r_ref <- pprimarycensored(
      d, pdiscretestep,
      dprimary = dexpgrowth, primary_args = list(r = 0.3),
      pwindow = 3.0,
      boundaries = boundaries_vals, pmf = pmf_vals
    )
    expect_equal(stan_val, r_ref,
      tolerance = 1e-6,
      info = sprintf("expgrowth mismatch at d = %s", d)
    )
  }
})
