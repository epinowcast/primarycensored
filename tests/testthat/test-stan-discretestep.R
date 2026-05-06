skip_on_cran()

# Simple step PMF: 3 intervals, boundaries 0:3
pmf_vals <- c(0.2, 0.5, 0.3)
boundaries_vals <- 0:3

# Packed params for dist_id 26: [boundaries (K+1=4), pmf (K=3)],
# total length = 7
packed_params <- c(as.numeric(boundaries_vals), pmf_vals)

# R reference: step CDF at a point t.
# Matches pstep_lcdf() in Stan: on interval [bk, bk+1) the CDF equals
# sum(pmf[1:k-1]) (cumulative before adding bin k's mass).
r_step_cdf <- function(t, boundaries, pmf) {
  K <- length(pmf)
  if (t < boundaries[1]) {
    return(0)
  }
  if (t >= boundaries[K + 1]) {
    return(1)
  }
  cum <- 0
  for (k in seq_len(K)) {
    if (t >= boundaries[k + 1]) {
      cum <- cum + pmf[k]
    } else if (t >= boundaries[k]) {
      cum <- cum + pmf[k]
      break
    }
  }
  cum
}

test_that("Stan pstep_lcdf matches R step CDF on a grid", {
  grid <- c(
    -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 2.999, 3, 3.5
  )
  for (t in grid) {
    r_val <- r_step_cdf(t, boundaries_vals, pmf_vals)
    r_log <- if (r_val == 0) -Inf else log(r_val)
    stan_val <- pstep_lcdf(t, as.numeric(boundaries_vals), pmf_vals)
    expect_equal(
      stan_val, r_log,
      tolerance = 1e-9,
      info = sprintf("pstep_lcdf mismatch at t = %s", t)
    )
  }
})

test_that(
  "Stan dist_lcdf with dist_id 26 matches pstep_lcdf directly",
  {
    grid <- c(0.5, 1.0, 1.5, 2.0, 2.5, 2.99)
    for (t in grid) {
      via_dist <- dist_lcdf(t, packed_params, 26L)
      direct <- pstep_lcdf(
        t, as.numeric(boundaries_vals), pmf_vals
      )
      expect_equal(
        via_dist, direct,
        tolerance = 1e-9,
        info = sprintf("dist_lcdf id=26 mismatch at t = %s", t)
      )
    }
  }
)

test_that("Stan hazards_to_pmf matches R reference", {
  # hazards: last must be 1 for the PMF to sum to 1
  hazards <- c(0.2, 0.4, 1.0)
  # R reference
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

test_that(
  "Stan analytical step+uniform matches R analytical for pwindow=1",
  {
    pwindow <- 1.0
    # Closed-form R reference matching primarycensored_discretestep_lcdf
    # in nonparametric.stan (uniform case).
    # F_obs(d) = (1/pwindow) * integral_{q}^{d} F_step(u) du,
    # q = max(d - pwindow, 0).
    # Uniform reduction: F_primary(d-lo) - F_primary(d-hi) = (hi-lo)/pwindow.
    r_analytical <- function(d, boundaries, pmf, pwindow) {
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
      analytical <- primarycensored_cdf(
        d, 26L, packed_params, pwindow, 0, Inf, 1L, numeric(0)
      )
      r_ref <- r_analytical(
        d, boundaries_vals, pmf_vals, pwindow
      )
      if (r_ref > 0) {
        expect_equal(
          analytical, r_ref,
          tolerance = 1e-5,
          info = sprintf(
            "analytical vs R reference mismatch at d = %s", d
          )
        )
      }
    }
  }
)

test_that(
  "Stan analytical step+expgrowth matches R partition reference",
  {
    pwindow <- 3.0
    r_val <- 0.3 # positive growth rate
    # R reference: F_obs(d) = integral_q^d F_step(u) f_expgrowth(d-u) du
    # using the change-of-variable result with F_primary.
    # expgrowth CDF on [0, pwindow]: F(p) = (exp(r*p) - 1)/(exp(r*pwindow) - 1)
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
        fp_end <- expgrowth_cdf_r(d - d, r, pwindow)
        integral <- integral + (fp_tail - fp_end)
      }
      integral
    }

    # primary_params for expgrowth: [r]
    primary_params_r <- c(r_val)

    d_vals <- c(0.5, 1.0, 2.0, 2.5, 2.99)
    for (d in d_vals) {
      stan_val <- primarycensored_cdf(
        d, 26L, packed_params,
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
            "step+expgrowth mismatch at d = %s", d
          )
        )
      }
    }
  }
)
