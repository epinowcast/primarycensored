skip_on_cran()

# Hazards over 3 intervals with boundaries 0:3.
# Last hazard must equal 1 so that the PMF sums to 1.
hazard_vals <- c(0.2, 0.4, 1.0)
boundaries_vals <- 0:3

# Packed params for dist_id 27:
# [boundaries (K+1=4), hazards (K=3)], total length = 7
packed_params_27 <- c(as.numeric(boundaries_vals), hazard_vals)

# PMF implied by the hazards (R-side reference)
pmf_vals <- hazards_to_pmf(hazard_vals)
packed_params_26 <- c(as.numeric(boundaries_vals), pmf_vals)

test_that(
  "Stan dist_lcdf with dist_id 27 matches dist_id 26 after hazard conversion",
  {
    for (t in c(0.5, 1.5, 2.5, 2.99)) {
      val_27 <- dist_lcdf(t, packed_params_27, 27L)
      val_26 <- dist_lcdf(t, packed_params_26, 26L)
      expect_equal(val_27, val_26,
        tolerance = 1e-9,
        info = sprintf("mismatch at t = %s", t)
      )
    }
  }
)

test_that(
  "Stan analytic hazard CDF matches R pprimarycensored (uniform primary)",
  {
    for (d in c(0.5, 1.0, 1.5, 2.0, 2.5, 2.99)) {
      stan_val <- primarycensored_cdf(
        d, 27L, packed_params_27, 1.0, 0, Inf, 1L, numeric(0)
      )
      r_ref <- pprimarycensored(
        d, pdiscretehazard,
        pwindow = 1.0,
        boundaries = boundaries_vals, hazards = hazard_vals
      )
      expect_equal(stan_val, r_ref,
        tolerance = 1e-6,
        info = sprintf("uniform mismatch at d = %s", d)
      )
    }
  }
)

test_that(
  "Stan dist_id 27 and 26 give same primary-censored CDF after conversion",
  {
    # The hazard route (27) must agree with the explicit-PMF route (26) for
    # any primary that 26 supports. This implicitly validates the analytic
    # path for non-uniform primaries via the 26 cross-check.
    for (primary in list(
      list(id = 1L, params = numeric(0), w = 1.0),
      list(id = 2L, params = c(0.3), w = 3.0)
    )) {
      for (d in c(0.5, 1.0, 2.0, 2.99)) {
        cdf_27 <- primarycensored_cdf(
          d, 27L, packed_params_27, primary$w, 0, Inf,
          primary$id, primary$params
        )
        cdf_26 <- primarycensored_cdf(
          d, 26L, packed_params_26, primary$w, 0, Inf,
          primary$id, primary$params
        )
        expect_equal(cdf_27, cdf_26,
          tolerance = 1e-9,
          info = sprintf(
            "id=27 vs id=26 mismatch at d = %s, primary = %d",
            d, primary$id
          )
        )
      }
    }
  }
)
