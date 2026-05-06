skip_on_cran()

# Hazards over 3 intervals with boundaries 0:3.
# Last hazard must equal 1 so that the PMF sums to 1.
hazard_vals <- c(0.2, 0.4, 1.0)
boundaries_vals <- 0:3

# Packed params for dist_id 27:
# [boundaries (K+1=4), hazards (K=3)], total length = 7
packed_params_27 <- c(as.numeric(boundaries_vals), hazard_vals)

# R reference: convert hazards to PMF
r_hazards_to_pmf <- function(hazards) {
  K <- length(hazards)
  pmf <- numeric(K)
  surv <- 1.0
  for (i in seq_len(K)) {
    pmf[i] <- hazards[i] * surv
    surv <- surv * (1 - hazards[i])
  }
  pmf
}

pmf_vals <- r_hazards_to_pmf(hazard_vals)

# Packed params for dist_id 26 (step-PMF equivalent)
packed_params_26 <- c(as.numeric(boundaries_vals), pmf_vals)

test_that(
  "Stan dist_lcdf with dist_id 27 matches dist_id 26 after hazard conversion",
  {
    grid <- c(0.5, 1.0, 1.5, 2.0, 2.5, 2.99)
    for (t in grid) {
      val_27 <- dist_lcdf(t, packed_params_27, 27L)
      val_26 <- dist_lcdf(t, packed_params_26, 26L)
      expect_equal(
        val_27, val_26,
        tolerance = 1e-9,
        info = sprintf("dist_lcdf id=27 vs id=26 mismatch at t = %s", t)
      )
    }
  }
)

test_that(
  "Stan analytical hazard+uniform matches R reference (pwindow=1)",
  {
    pwindow <- 1.0
    # R reference: same as step formula but with PMF derived from hazards.
    # Uses (hi - lo) / pwindow for the uniform primary CDF difference.
    r_analytical_uniform <- function(d, boundaries, pmf, pwindow) {
      K <- length(pmf)
      q_lo <- max(d - pwindow, 0)
      cumulative <- 0
      integral <- 0
      for (k in seq_len(K)) {
        bk <- boundaries[k]
        bk1 <- boundaries[k + 1]
        if (bk1 <= q_lo) {
          cumulative <- cumulative + pmf[k]
          next
        }
        if (bk >= d) break
        lo <- max(q_lo, bk)
        hi <- min(d, bk1)
        if (hi > lo) integral <- integral + cumulative * (hi - lo)
        cumulative <- cumulative + pmf[k]
      }
      tail_start <- max(boundaries[K + 1], q_lo)
      if (tail_start < d) integral <- integral + (d - tail_start)
      integral / pwindow
    }

    d_vals <- c(0.5, 1.0, 1.5, 2.0, 2.5, 2.99)
    for (d in d_vals) {
      stan_val <- primarycensored_cdf(
        d, 27L, packed_params_27, pwindow, 0, Inf, 1L, numeric(0)
      )
      r_ref <- r_analytical_uniform(
        d, boundaries_vals, pmf_vals, pwindow
      )
      if (r_ref > 0) {
        expect_equal(
          stan_val, r_ref,
          tolerance = 1e-5,
          info = sprintf(
            "hazard+uniform mismatch at d = %s", d
          )
        )
      }
    }
  }
)

test_that(
  "Stan analytical hazard+expgrowth matches R partition reference",
  {
    pwindow <- 3.0
    r_val <- 0.3

    expgrowth_cdf_r <- function(p, r, pwindow) {
      if (p <= 0) {
        return(0)
      }
      if (p >= pwindow) {
        return(1)
      }
      (exp(r * p) - 1) / (exp(r * pwindow) - 1)
    }

    r_analytical_expgrowth <- function(d, boundaries, pmf, r, pwindow) {
      K <- length(pmf)
      q_lo <- max(d - pwindow, 0)
      cumulative <- 0
      integral <- 0
      for (k in seq_len(K)) {
        bk <- boundaries[k]
        bk1 <- boundaries[k + 1]
        if (bk1 <= q_lo) {
          cumulative <- cumulative + pmf[k]
          next
        }
        if (bk >= d) break
        lo <- max(q_lo, bk)
        hi <- min(d, bk1)
        if (hi > lo) {
          fp_lo <- expgrowth_cdf_r(d - lo, r, pwindow)
          fp_hi <- expgrowth_cdf_r(d - hi, r, pwindow)
          integral <- integral + cumulative * (fp_lo - fp_hi)
        }
        cumulative <- cumulative + pmf[k]
      }
      tail_start <- max(boundaries[K + 1], q_lo)
      if (tail_start < d) {
        fp_tail <- expgrowth_cdf_r(d - tail_start, r, pwindow)
        fp_end <- expgrowth_cdf_r(0, r, pwindow)
        integral <- integral + (fp_tail - fp_end)
      }
      integral
    }

    primary_params_r <- c(r_val)
    d_vals <- c(0.5, 1.0, 2.0, 2.5, 2.99)
    for (d in d_vals) {
      stan_val <- primarycensored_cdf(
        d, 27L, packed_params_27,
        pwindow, 0, Inf,
        2L, primary_params_r
      )
      r_ref <- r_analytical_expgrowth(
        d, boundaries_vals, pmf_vals, r_val, pwindow
      )
      if (r_ref > 0) {
        expect_equal(
          stan_val, r_ref,
          tolerance = 1e-5,
          info = sprintf(
            "hazard+expgrowth mismatch at d = %s", d
          )
        )
      }
    }
  }
)

test_that(
  "Stan dist_id 27 and 26 give same result when hazard PMF equals PMF",
  {
    # Verify that packing the same PMF via hazards route (dist_id=27) gives
    # identical CDF results to the direct PMF route (dist_id=26).
    pwindow <- 1.0
    d_vals <- c(0.5, 1.0, 2.5)
    for (d in d_vals) {
      cdf_27 <- primarycensored_cdf(
        d, 27L, packed_params_27, pwindow, 0, Inf, 1L, numeric(0)
      )
      cdf_26 <- primarycensored_cdf(
        d, 26L, packed_params_26, pwindow, 0, Inf, 1L, numeric(0)
      )
      expect_equal(
        cdf_27, cdf_26,
        tolerance = 1e-9,
        info = sprintf(
          "id=27 vs id=26 primary CDF mismatch at d = %s", d
        )
      )
    }
  }
)
