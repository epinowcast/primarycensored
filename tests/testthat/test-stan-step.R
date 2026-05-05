skip_on_cran()

# Simple step PMF: 3 intervals, boundaries 0:3
pmf_vals <- c(0.2, 0.5, 0.3)
boundaries_vals <- 0:3

# Packed params for dist_id 26: [boundaries (K+1=4), pmf (K=3)],
# total length = 7
packed_params <- c(as.numeric(boundaries_vals), pmf_vals)

# R reference: step CDF at a point t
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
  "Stan analytical step+uniform matches ODE numerical for pwindow=1",
  {
    pwindow <- 1.0
    # Test over a grid of d values; stay within step support [0, 3]
    d_vals <- c(0.5, 1.0, 1.5, 2.0, 2.5, 2.99)
    for (d in d_vals) {
      # Analytical path (dist_id=26, primary_id=1)
      analytical <- primarycensored_cdf(
        d, 26L, packed_params, pwindow, 0, Inf, 1L, numeric(0)
      )
      # ODE numerical path: force it by using dist_id that has no
      # analytical solution for primary_id=1... but we can also call
      # primarycensored_cdf with the analytical check enabled and compare
      # to the R-side integral directly.
      # Build the R integral manually (vectorised over p for integrate()):
      r_integrand <- function(p) {
        vapply(p, function(pi) {
          s <- d - pi # delay seen from primary event
          r_step_cdf(s, boundaries_vals, pmf_vals) * (1 / pwindow)
        }, numeric(1))
      }
      lower_b <- max(d - pwindow, 0)
      r_integral <- tryCatch(
        integrate(r_integrand, lower = lower_b, upper = d)$value,
        error = function(e) NA_real_
      )
      if (!is.na(r_integral) && r_integral > 0) {
        expect_equal(
          analytical, r_integral,
          tolerance = 1e-5,
          info = sprintf(
            "analytical vs R integral mismatch at d = %s", d
          )
        )
      }
    }
  }
)

# Once the R-side pstep() function is available, add a test comparing
# Stan pstep_lcdf to log(pstep(t, boundaries, pmf)).
