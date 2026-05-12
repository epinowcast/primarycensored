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
  * Vectorised PMF reduction via `cumulative_sum`. The bin-search index
  * runs on `data`-level inputs (`t`, `boundaries`) and never appears on
  * the autodiff tape, so a small data-only loop is kept here.
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
  if (t < boundaries[2]) return negative_infinity();
  if (t >= boundaries[K + 1]) return 0;
  // Right-continuous CDF with jumps at the right edges
  // boundaries[2], ..., boundaries[K + 1]. F(t) = cum_pmf[k] for
  // t in [boundaries[k + 1], boundaries[k + 2]); equivalently the
  // largest k with boundaries[k + 1] <= t. Boundary-on-jump cases
  // (t == boundaries[k + 1]) advance k, matching R's
  // `findInterval(left.open = FALSE)`.
  int k = 1;
  while (k < K && boundaries[k + 2] <= t) k += 1;
  return log(cumulative_sum(pmf)[k]);
}

/**
  * Convert discrete hazards to a PMF
  *
  * Each entry satisfies pmf[i] = hazards[i] * prod_{j < i} (1 - hazards[j]).
  * The last hazard must equal 1 so the PMF sums to 1 (caller's
  * responsibility). One `log1m`, one `cumulative_sum`, one `exp` -- no
  * per-bin loop on the autodiff tape.
  *
  * @param hazards Vector of K hazards in [0, 1], with hazards[K] = 1
  *
  * @return PMF vector of length K
  */
vector hazards_to_pmf(vector hazards) {
  int K = num_elements(hazards);
  vector[K] log_surv;
  log_surv[1] = 0;
  if (K > 1) {
    log_surv[2:K] = cumulative_sum(log1m(hazards[1:(K - 1)]));
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
  * Vectorised primary log CDF
  *
  * Element-wise wrapper around `primary_lcdf`. Lets the analytic step
  * convolution pull `f_lo` and `f_hi` out of one pair of vector calls
  * without per-bin branching at the reduction site.
  *
  * @param p Vector of primary event times in [0, pwindow]
  * @param primary_id Primary distribution identifier
  * @param primary_params Distribution parameters
  * @param pwindow Primary event window width
  *
  * @return Vector of `log(F_primary(p))`
  */
vector primary_lcdf_vec(vector p, int primary_id,
                        array[] real primary_params, data real pwindow) {
  int N = num_elements(p);
  vector[N] out;
  for (i in 1:N) {
    out[i] = primary_lcdf(p[i] | primary_id, primary_params, pwindow);
  }
  return out;
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
  * The lo/hi/cumulative vectors are built in one pass; `f_lo` and `f_hi`
  * come from two vectorised `primary_lcdf_vec` calls; an `active` mask
  * zeros out empty sub-intervals; the per-bin contributions are reduced
  * via `dot_product`. The tail [boundaries[K+1], d] (where F_step = 1)
  * is added in a single closed-form term, no loop.
  *
  * Boundary case: when the integration support lies entirely at or below
  * `boundaries[2]` the only sub-interval that overlaps is bin 1 (with
  * `cum_before = 0`) and the tail is empty, so the integral is
  * structurally zero. We return `negative_infinity()` directly to keep
  * `log(0)` off the autodiff tape; the value has no parameter dependence
  * in this regime so the gradient is zero, and downstream
  * `log_diff_exp(a, -inf)` evaluates cleanly to `a`.
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

  // Structural-zero short-circuit. Below `boundaries[2]` F_step is zero
  // and the bin-1 contribution carries `cum_before = 0`, so the integral
  // collapses to 0. Returning `negative_infinity()` directly keeps
  // `log(0)` off the autodiff tape so downstream `log_diff_exp(a, -inf)`
  // evaluates cleanly with a zero gradient w.r.t. `pmf`.
  if (u_max <= boundaries[2]) return negative_infinity();

  // Sub-interval endpoints in u-space, clipped to [u_min, u_max].
  vector[K] lo = fmax(u_min, head(boundaries, K));
  vector[K] hi = fmin(u_max, tail(boundaries, K));

  // F_step is right-continuous and on [b_k, b_{k+1}) takes the value
  // sum_{j < k} pmf[j] (mass before bin k). cumulative_sum(pmf) gives
  // the mass through and including bin k, so we shift right by one.
  vector[K] cum_before;
  cum_before[1] = 0;
  if (K > 1) cum_before[2:K] = head(cumulative_sum(pmf), K - 1);

  // 0/1 mask drops bins with `hi <= lo` from the reduction without a
  // branch in the inner expression. Built on `data`-level inputs.
  vector[K] active;
  for (k in 1:K) active[k] = hi[k] > lo[k] ? 1 : 0;

  // F_primary at lo/hi via two vectorised calls; one masked subtraction
  // gives the per-bin difference for the dot product.
  vector[K] f_lo = primary_lcdf_vec(d - lo, primary_id, primary_params,
                                    pwindow);
  vector[K] f_hi = primary_lcdf_vec(d - hi, primary_id, primary_params,
                                    pwindow);
  vector[K] f_diff = (exp(f_lo) - exp(f_hi)) .* active;

  real integral = dot_product(cum_before, f_diff);

  // Tail region [boundaries[K+1], u_max]: F_step = 1, contributing
  // F_primary(d - tail_start) - F_primary(d - u_max).
  real tail_start = fmax(boundaries[K + 1], u_min);
  if (tail_start < u_max) {
    real fp_tail = exp(primary_lcdf(d - tail_start | primary_id,
                                    primary_params, pwindow));
    real fp_end = exp(primary_lcdf(d - u_max | primary_id,
                                   primary_params, pwindow));
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
