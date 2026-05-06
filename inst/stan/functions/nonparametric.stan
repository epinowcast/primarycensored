/**
  * Non-parametric step CDF and hazard conversion utilities.
  *
  * The step CDF is defined by K intervals and a PMF over those intervals.
  * Boundaries is a vector of length K+1 giving interval endpoints
  * [boundaries[1], boundaries[2]), ..., [boundaries[K], boundaries[K+1]).
  * pmf is a simplex of length K giving the probability mass in each interval.
  */

/**
  * Compute the log CDF of a piecewise-constant (step) distribution
  *
  * @param t Evaluation point
  * @param boundaries Vector of K+1 interval endpoints (strictly increasing)
  * @param pmf Simplex of K probabilities, one per interval (must sum to 1)
  *
  * @return log(F_step(t)):
  *   negative_infinity() if t < boundaries[1],
  *   0 if t >= boundaries[K+1],
  *   log of cumulative mass up to the interval containing t otherwise.
  */
real pstep_lcdf(real t, vector boundaries, vector pmf) {
  int K = num_elements(pmf);
  // t is below the support
  if (t < boundaries[1]) return negative_infinity();
  // t is at or above the upper boundary
  if (t >= boundaries[K + 1]) return 0;

  // Accumulate mass for intervals that are fully below t
  real cumulative = 0;
  for (k in 1:K) {
    if (t >= boundaries[k + 1]) {
      cumulative += pmf[k];
    } else if (t >= boundaries[k]) {
      // t is inside interval k; the full mass of this interval is included
      // because the step CDF is right-continuous and we have F(t) = sum up to
      // and including the interval containing t
      cumulative += pmf[k];
      break;
    }
  }
  return log(cumulative);
}

/**
  * Convert discrete hazards to a PMF
  *
  * Each entry satisfies: pmf[i] = hazards[i] * prod_{j < i}(1 - hazards[j]).
  * The last hazard must equal 1 so that the PMF sums to 1 (caller's
  * responsibility).
  *
  * @param hazards Vector of K hazards in [0, 1], with hazards[K] = 1
  *
  * @return PMF vector of length K
  */
vector hazards_to_pmf(vector hazards) {
  int K = num_elements(hazards);
  vector[K] pmf;
  real survival = 1.0;
  for (k in 1:K) {
    pmf[k] = hazards[k] * survival;
    survival *= (1.0 - hazards[k]);
  }
  return pmf;
}

/**
  * Numerically stable conversion from log-hazards to a PMF
  *
  * @param log_hazards Vector of K log-hazards, log(h_k)
  * @param log1m_hazards Vector of K log(1 - h_k) values
  *
  * @return PMF in primal scale, length K
  */
vector log_hazards_to_pmf(vector log_hazards, vector log1m_hazards) {
  int K = num_elements(log_hazards);
  vector[K] pmf;
  real log_survival = 0.0; // log(1) = 0
  for (k in 1:K) {
    pmf[k] = exp(log_hazards[k] + log_survival);
    log_survival += log1m_hazards[k];
  }
  return pmf;
}

/**
  * Compute the primary event censored log CDF analytically for a step delay
  * with any supported primary distribution.
  *
  * The formula integrates F_step over the primary event window using:
  *   F_obs(d) = integral_0^{pwindow} F_step(d - p) f_primary(p) dp.
  * Substituting u = d - p (dp = -du), with clamping to [0, pwindow]:
  *   u_min = max(d - pwindow, 0) = q, u_max = d.
  *   F_obs(d) = integral_{q}^{d} F_step(u) f_primary(d - u) du.
  * Because f_primary(d - u) du = -dF_primary(d - u) we can write:
  *   F_obs(d) = integral_{q}^{d} F_step(u) dF_primary(d - u).
  * On each sub-interval [lo, hi] where F_step equals the constant
  * `cumulative`, the contribution is:
  *   cumulative * [F_primary(d - lo) - F_primary(d - hi)].
  * (The sign reversal comes from the change-of-variable direction.)
  *
  * Reduction to the uniform case: for uniform primary,
  * F_primary(p) = p / pwindow, so
  *   F_primary(d - lo) - F_primary(d - hi) = (hi - lo) / pwindow.
  * Summing gives integral / pwindow, recovering the original formula.
  *
  * @param d Delay time (observation point)
  * @param q Lower bound of integration: max(d - pwindow, 0)
  * @param boundaries Vector of K+1 step boundaries
  * @param pmf Step PMF of length K
  * @param primary_id Primary distribution identifier (1=uniform, 2=expgrowth)
  * @param primary_params Primary distribution parameters
  * @param pwindow Primary event window width
  *
  * @return log(F_obs(d)) under convolution with the given primary distribution
  */
real primarycensored_discretestep_lcdf(
  data real d, real q,
  vector boundaries, vector pmf,
  int primary_id, array[] real primary_params,
  data real pwindow
) {
  int K = num_elements(pmf);
  // F_step is right-continuous: on [bk, bk1) the value is
  // sum(pmf[1:k-1]) (the jump at bk1 has not yet occurred).

  real u_min = q; // lower limit of integration
  real u_max = d; // upper limit of integration
  real integral = 0;
  real cumulative = 0; // running sum of pmf

  for (k in 1:K) {
    // Step interval k: [bk, bk1); F_step = cumulative (before adding pmf[k])
    real bk  = boundaries[k];
    real bk1 = boundaries[k + 1];

    // Skip intervals entirely below u_min
    if (bk1 <= u_min) {
      cumulative += pmf[k];
      continue;
    }
    // Stop once we are entirely above u_max
    if (bk >= u_max) break;

    // Active sub-interval [lo, hi] in u-space
    real lo = fmax(u_min, bk);
    real hi = fmin(u_max, bk1);

    if (hi > lo) {
      // Contribution: cumulative * [F_primary(d-lo) - F_primary(d-hi)]
      real fp_lo = primary_F(d - lo, primary_id, primary_params, pwindow);
      real fp_hi = primary_F(d - hi, primary_id, primary_params, pwindow);
      integral += cumulative * (fp_lo - fp_hi);
    }

    cumulative += pmf[k]; // F_step = cumulative on [bk1, next_bk1)
    if (hi >= u_max) break;
  }

  // Tail region u in [boundaries[K+1], u_max]: F_step = 1
  real tail_start = fmax(boundaries[K + 1], u_min);
  if (tail_start < u_max) {
    real fp_tail = primary_F(
      d - tail_start, primary_id, primary_params, pwindow
    );
    real fp_end = primary_F(
      d - u_max, primary_id, primary_params, pwindow
    );
    integral += (fp_tail - fp_end);
  }

  return log(integral);
}
