/**
  * Non-parametric step CDF and hazard conversion utilities.
  *
  * The step CDF is defined by K intervals and a PMF over those intervals.
  * `boundaries` is a vector of length K+1 giving interval endpoints
  * [boundaries[1], boundaries[2]), ..., [boundaries[K], boundaries[K+1]).
  * `pmf` is a simplex of length K giving the probability mass in each
  * interval. The hazard parameterisation replaces `pmf` with discrete
  * hazards in [0, 1] whose final entry is 1 so the implied PMF sums to 1.
  */

/**
  * Log CDF of a piecewise-constant (step) distribution
  *
  * Vectorised via `cumulative_sum`; the only loop locates the bin
  * containing `t` and is bounded by the (small) number of bins K.
  *
  * @param t Evaluation point
  * @param boundaries Vector of K+1 interval endpoints (strictly increasing)
  * @param pmf Simplex of K probabilities, one per interval (must sum to 1)
  *
  * @return log(F_step(t)):
  *   negative_infinity() if t < boundaries[1],
  *   0 if t >= boundaries[K + 1],
  *   log of the cumulative mass through the bin containing t otherwise.
  */
real pstep_lcdf(real t, vector boundaries, vector pmf) {
  int K = num_elements(pmf);
  if (t < boundaries[1]) return negative_infinity();
  if (t >= boundaries[K + 1]) return 0;
  vector[K] cum = cumulative_sum(pmf);
  // Locate the bin containing t: largest k with boundaries[k] <= t.
  int k = 1;
  while (k < K && boundaries[k + 1] <= t) k += 1;
  return log(cum[k]);
}

/**
  * Convert discrete hazards to a PMF
  *
  * Each entry satisfies pmf[i] = hazards[i] * prod_{j < i} (1 - hazards[j]).
  * The last hazard must equal 1 so the PMF sums to 1 (caller's
  * responsibility).
  *
  * @param hazards Vector of K hazards in [0, 1], with hazards[K] = 1
  *
  * @return PMF vector of length K
  */
vector hazards_to_pmf(vector hazards) {
  int K = num_elements(hazards);
  // Survival before bin k: S[k] = prod_{j < k} (1 - hazards[j]).
  // S[1] = 1; S[k] = exp(cumulative_sum(log1m_hazards)[k - 1]) for k >= 2.
  vector[K] log1m_h = log1m(hazards);
  vector[K] log_surv;
  log_surv[1] = 0;
  if (K > 1) {
    log_surv[2:K] = cumulative_sum(log1m_h[1:(K - 1)]);
  }
  return hazards .* exp(log_surv);
}

/**
  * Log CDF of a discrete-hazard distribution
  *
  * Sibling of `pstep_lcdf` for the hazard parameterisation; converts
  * hazards to the implied PMF then dispatches to `pstep_lcdf`.
  *
  * @param t Evaluation point
  * @param boundaries Vector of K+1 boundaries (strictly increasing)
  * @param hazards Vector of K hazards in [0, 1] with hazards[K] = 1
  *
  * @return log(F_step(t)) under the implied PMF.
  */
real phazard_lcdf(real t, vector boundaries, vector hazards) {
  return pstep_lcdf(t | boundaries, hazards_to_pmf(hazards));
}

/**
  * Primary event censored log CDF for a step delay (vectorised analytic)
  *
  * Computes log(F_obs(d)) where
  *   F_obs(d) = integral_q^d F_step(u) f_primary(d - u) du,
  *   q = max(d - pwindow, 0).
  * Using f_primary(d - u) du = -d F_primary(d - u), the integral on each
  * sub-interval [lo, hi] where F_step is constant becomes
  *   cumulative * (F_primary(d - lo) - F_primary(d - hi)).
  * This routine builds the lo/hi/cumulative vectors in one pass and reduces
  * via `dot_product`. The tail [boundaries[K+1], d] (where F_step = 1) is
  * added in a single closed-form term.
  *
  * @param d Delay (observation point)
  * @param boundaries Vector of K+1 step boundaries
  * @param pmf Step PMF of length K
  * @param primary_id Primary distribution identifier
  * @param primary_params Primary distribution parameters
  * @param pwindow Primary event window width
  *
  * @return log(F_obs(d)) under convolution with the given primary
  */
real discretestep_lcdf(
  data real d, vector boundaries, vector pmf,
  int primary_id, array[] real primary_params, data real pwindow
) {
  int K = num_elements(pmf);
  real u_min = fmax(d - pwindow, 0);
  real u_max = d;

  // Sub-interval endpoints in u-space, clipped to [u_min, u_max]
  vector[K] lo = fmax(u_min, head(boundaries, K));
  vector[K] hi = fmin(u_max, tail(boundaries, K));

  // F_step is right-continuous and on [b_k, b_{k+1}) takes the value
  // sum_{j < k} pmf[j] (mass before bin k). cumulative_sum(pmf) gives the
  // mass through and including bin k, so we shift right by one.
  vector[K] cum_after = cumulative_sum(pmf);
  vector[K] cum_before;
  cum_before[1] = 0;
  if (K > 1) cum_before[2:K] = cum_after[1:(K - 1)];

  // F_primary differences per sub-interval. Bins where hi <= lo contribute
  // zero, so we only call primary_lcdf when hi > lo.
  vector[K] f_diff = rep_vector(0, K);
  for (k in 1:K) {
    if (hi[k] > lo[k]) {
      real fp_lo = exp(primary_lcdf(d - lo[k] | primary_id, primary_params,
                                    pwindow));
      real fp_hi = exp(primary_lcdf(d - hi[k] | primary_id, primary_params,
                                    pwindow));
      f_diff[k] = fp_lo - fp_hi;
    }
  }

  real integral = dot_product(cum_before, f_diff);

  // Tail region [boundaries[K+1], u_max]: F_step = 1, contributing
  // F_primary(d - tail_start) - F_primary(d - u_max).
  real tail_start = fmax(boundaries[K + 1], u_min);
  if (tail_start < u_max) {
    real fp_tail = exp(primary_lcdf(d - tail_start | primary_id,
                                    primary_params, pwindow));
    real fp_end = exp(primary_lcdf(d - u_max | primary_id, primary_params,
                                   pwindow));
    integral += fp_tail - fp_end;
  }

  return log(integral);
}

/**
  * Primary event censored log CDF for a discrete-hazard delay
  *
  * Wrapper that converts hazards to a PMF and delegates to
  * `discretestep_lcdf`. Provides the analytic CDF for `dist_id == 27`.
  *
  * @param d Delay (observation point)
  * @param boundaries Vector of K+1 step boundaries
  * @param hazards Vector of K hazards in [0, 1] with hazards[K] = 1
  * @param primary_id Primary distribution identifier
  * @param primary_params Primary distribution parameters
  * @param pwindow Primary event window width
  *
  * @return log(F_obs(d)) under convolution with the given primary
  */
real discretehazard_lcdf(
  data real d, vector boundaries, vector hazards,
  int primary_id, array[] real primary_params, data real pwindow
) {
  return discretestep_lcdf(
    d | boundaries, hazards_to_pmf(hazards), primary_id, primary_params,
    pwindow
  );
}
