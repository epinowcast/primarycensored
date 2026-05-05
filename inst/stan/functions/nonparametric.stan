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
  * with a uniform primary distribution.
  *
  * The analytic formula integrates F_step over the primary event window
  * [max(d - pwindow, 0), d] with a uniform primary density 1/pwindow.
  * Because F_step is piecewise-constant, the integral reduces to a sum over
  * sub-intervals determined by the step boundaries.
  *
  * @param d Delay time (observation point)
  * @param q Lower bound of integration: max(d - pwindow, 0)
  * @param boundaries Vector of K+1 step boundaries
  * @param pmf Step PMF of length K
  * @param pwindow Primary event window width
  *
  * @return log(F_obs(d)) under convolution with uniform primary censoring
  */
real primarycensored_step_uniform_lcdf(
  data real d, real q,
  vector boundaries, vector pmf, data real pwindow
) {
  int K = num_elements(pmf);
  // Integration is over p in [q, d]; shift variable s = d - p is the delay,
  // so we integrate F_step(s) * (1/pwindow) ds for s from 0 to d-q.
  // Equivalently, integrate over p in [q, d] using the uniform density.
  //
  // The integral (1/pwindow) * integral_{q}^{d} F_step(d - p) dp
  // where d - p ranges from 0 (at p=d) to d - q (at p=q).
  // Setting u = d - p: (1/pwindow) * integral_{0}^{d-q} F_step(u) du
  //
  // F_step(u) is constant within each boundary interval, so partition [0, d-q]
  // at the step boundaries and evaluate each piece.

  real u_max = d - q; // upper limit of integration in u
  real integral = 0;
  real prev = 0; // lower end of current sub-interval in u

  for (k in 1:K) {
    // boundaries for u: step interval k spans [boundaries[k], boundaries[k+1])
    real bk = boundaries[k];
    real bk1 = boundaries[k + 1];

    // The sub-interval in u is [max(prev, bk), min(u_max, bk1)]
    real lo = fmax(prev, bk);
    real hi = fmin(u_max, bk1);

    if (hi <= lo) {
      if (bk >= u_max) break;
      continue;
    }

    // Cumulative PMF at step k (F_step equals sum of pmf[1..k] on [bk, bk1))
    real cum_k = sum(pmf[1:k]);
    integral += cum_k * (hi - lo);
    prev = hi;
    if (hi >= u_max) break;
  }

  // Also account for u >= boundaries[K+1]: F_step = 1
  if (prev < u_max) {
    integral += (u_max - prev);
  }

  integral /= pwindow;
  return log(integral);
}
