skip_on_cran()

# Simple step PMF: 3 intervals, boundaries 0:3
pmf_vals <- c(0.2, 0.5, 0.3)
boundaries_vals <- 0:3

# Packed params for dist_id 26: [boundaries (K+1=4), pmf (K=3)],
# total length = 7
packed_params <- c(as.numeric(boundaries_vals), pmf_vals)

test_that("Stan dist_lcdf with dist_id 26 matches pstep_lcdf directly", {
  for (t in c(0.5, 1.5, 2.5, 2.99)) {
    via_dist <- dist_lcdf(t, packed_params, 26L)
    direct <- pstep_lcdf(t, as.numeric(boundaries_vals), pmf_vals)
    expect_equal(via_dist, direct, tolerance = 1e-9)
  }
})

test_that("Stan pstep_lcdf matches R pdiscretestep at the bin edges", {
  # Right-continuous CDF: F jumps up at each right edge, so
  # F(boundaries[i + 1]) = cumsum(pmf)[i]. The previous Stan
  # implementation advanced past the jump and returned cumsum[i + 1].
  edges <- as.numeric(boundaries_vals[-1])
  cum_pmf <- cumsum(pmf_vals)
  for (i in seq_along(edges)) {
    stan_val <- exp(pstep_lcdf(
      edges[i], as.numeric(boundaries_vals), pmf_vals
    ))
    r_val <- pdiscretestep(
      edges[i], boundaries_vals, pmf_vals
    )
    expect_equal(stan_val, cum_pmf[i], tolerance = 1e-9)
    expect_equal(stan_val, r_val, tolerance = 1e-9)
  }
})

test_that("Stan analytic step CDF matches R pprimarycensored", {
  # Cross-check the analytic Stan path against the R reference
  # `pprimarycensored` for both supported primaries.
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
      d, 26L, packed_params, 3.0, 0, Inf, 2L, 0.3
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

test_that(
  "primarycensored_lpmf returns finite mass at the d=1 boundary case",
  {
    # Regression test for the `log(0) = -Inf` autodiff propagation in
    # `discretestep_lcdf`: the lower-edge interval [0, 1] has zero
    # cumulative-before mass, so `log_diff_exp(lcdf(d_upper), lcdf(d))`
    # evaluates against a literal -inf rather than a `log(0)` autodiff
    # value. The vignette's `delay > 1` workaround should be unnecessary.
    d_upper <- 2L
    D <- 3.0
    stan_val <- primarycensored_lpmf(
      1L, 26L, packed_params, 1.0, d_upper, 0.0, D, 1L, numeric(0)
    )
    expect_true(is.finite(stan_val))

    # Cross-validate against R's PMF construction.
    r_ref <- log(
      dprimarycensored(
        1, pdiscretestep,
        pwindow = 1.0, swindow = 1.0, D = D,
        boundaries = boundaries_vals, pmf = pmf_vals
      )
    )
    expect_equal(stan_val, r_ref, tolerance = 1e-9)
  }
)

test_that("discretestep_lcdf handles d=K and d_upper=D boundary cases", {
  # d at the upper support edge (K = 3): convolution still includes the
  # bulk of bin K under the primary smear, so the result must be finite
  # and agree with `pprimarycensored`.
  for (d in c(3, 3 - 1e-6)) {
    stan_val <- primarycensored_cdf(
      d, 26L, packed_params, 1.0, 0, Inf, 1L, numeric(0)
    )
    r_ref <- pprimarycensored(
      d, pdiscretestep,
      pwindow = 1.0, boundaries = boundaries_vals, pmf = pmf_vals
    )
    expect_equal(stan_val, r_ref, tolerance = 1e-6)
  }

  # d_upper = D: `primarycensored_lpmf` must reuse the cached cdf at the
  # upper truncation point and return finite mass.
  stan_val <- primarycensored_lpmf(
    2L, 26L, packed_params, 1.0, 3.0, 0.0, 3.0, 1L, numeric(0)
  )
  expect_true(is.finite(stan_val))
})
