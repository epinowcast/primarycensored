/**
  * Compute the log of the leading term of the incomplete gamma series
  * @ingroup delay_log_cdfs
  *
  * Returns log(x^a exp(-x) / Gamma(a + 1)) for x = exp(log_x).
  * For `a >= 100` the Stirling series of lgamma(a + 1) cancels the terms of
  * size a analytically, so rounding error does not grow with a.
  *
  * @param log_x Log of the argument, log(x) with x > 0
  * @param a Shape parameter of the Gamma distribution (a > 0)
  *
  * @return log(x^a exp(-x) / Gamma(a + 1))
  */
real gamma_log_lead_logx(real log_x, real a) {
  real x = exp(log_x);
  if (a < 100) {
    return a * log_x - x - lgamma(a + 1);
  }
  real inv_a = inv(a);
  real inv_a2 = square(inv_a);
  real delta = inv_a * (1.0 / 12
                        - inv_a2 * (1.0 / 360
                                    - inv_a2 * (1.0 / 1260
                                                - inv_a2 / 1680)));
  real r = x / a;
  // r - 1 is exact near 1, so log1p avoids cancellation
  real log_r = r > 0.5 && r < 2 ? log1p(r - 1) : log_x - log(a);
  return a * (log_r - (r - 1))
         - 0.5 * (log(a) + log(2 * pi())) - delta;
}

/**
  * Compute the log CDF and log CCDF of a unit rate Gamma distribution from
  * the log of x
  * @ingroup delay_log_cdfs
  *
  * Returns log P(a, x) and log Q(a, x) = log(1 - P(a, x)), the logs of the
  * regularised incomplete gamma functions, for x = exp(log_x).
  * For x < a + 1, P comes from its series.
  * Otherwise Q comes from its continued fraction (modified Lentz).
  * The other tail is `log1m_exp()` of the first, which is accurate because
  * the first is not close to 1 in its region.
  * Both stay finite, with accurate gradients, when x underflows or the
  * other tail is below the smallest double.
  *
  * @param log_x Log of the argument, log(x) with x > 0
  * @param a Shape parameter of the Gamma distribution (a > 0)
  *
  * @return Vector [log P(a, exp(log_x)), log Q(a, exp(log_x))]. Rejects
  * a <= 0, `nan` and `inf`.
  */
vector gamma_log_pq_logx(real log_x, real a) {
  if (!(a > 0) || is_inf(a)) {
    reject("gamma_log_pq_logx: shape must be positive finite, found a = ", a);
  }
  real x = exp(log_x);
  if (log_x == negative_infinity()) {
    return [negative_infinity(), 0]';
  }
  if (is_inf(x)) {
    return [0, negative_infinity()]';
  }
  real log_lead = gamma_log_lead_logx(log_x, a);
  // Both need about 9 sqrt(a) terms at worst, near x = a + 1
  real max_terms = 20 * sqrt(a) + 200;
  if (x < a + 1) {
    real term = 1;
    real total = 1;
    int n = 0;
    while (term >= 1e-17 * total) {
      n += 1;
      if (n > max_terms) {
        reject("gamma_log_pq_logx: series did not converge for log_x = ",
               log_x, " and a = ", a);
      }
      term *= x / (a + n);
      total += term;
    }
    real log_p = log_lead + log(total);
    return [log_p, log1m_exp(log_p)]';
  }
  real b = x + 1 - a;
  real c = 1e300;
  real d = 1 / b;
  real h = d;
  real delta = 1;
  // err tracks |delta - 1|, which each step multiplies by
  // |an * d * d_prev / delta_prev|. Flooring |i - a| in an at 1 keeps it
  // from reaching 0 at an integer shape, where the fraction terminates but
  // its shape gradient does not.
  real err = 1;
  int i = 0;
  while (err >= 1e-15) {
    i += 1;
    if (i > max_terms) {
      reject("gamma_log_pq_logx: continued fraction did not converge for ",
             "log_x = ", log_x, " and a = ", a);
    }
    real an = -i * (i - a);
    real d_prev = d;
    b += 2;
    d = an * d + b;
    if (abs(d) < 1e-300) d = 1e-300;
    c = b + an / c;
    if (abs(c) < 1e-300) c = 1e-300;
    d = 1 / d;
    err *= i * fmax(abs(i - a), 1) * abs(d * d_prev / delta);
    delta = d * c;
    h *= delta;
  }
  real log_q = log_lead + log(a) + log(h);
  return [log1m_exp(log_q), log_q]';
}

/**
  * Compute the log CDF of a unit rate Gamma distribution from the log of x
  * @ingroup delay_log_cdfs
  *
  * @param log_x Log of the argument, log(x) with x > 0
  * @param a Shape parameter of the Gamma distribution (a > 0)
  *
  * @return log P(a, exp(log_x)), the first element of
  * `gamma_log_pq_logx()`
  */
real gamma_lcdf_logx(real log_x, real a) {
  return gamma_log_pq_logx(log_x, a)[1];
}

/**
  * Compute the log CDFs of Gamma(a) and Gamma(a + 1) from the log of x
  * @ingroup delay_log_cdfs
  *
  * Returns [log P(a, x), log P(a + 1, x)] with unit rate, as needed by the
  * uniform primary solution for a Gamma delay.
  * Uses one `gamma_log_pq_logx()` evaluation and the recursion
  * P(a + 1, x) = P(a, x) - x^a exp(-x) / Gamma(a + 1).
  * Below a + 1 it starts from P(a + 1, x) and above from Q(a, x), so that
  * the recursion adds terms rather than subtracting them.
  *
  * @param log_x Log of the argument, log(x) with x > 0
  * @param a Shape parameter of the Gamma distribution (a > 0)
  *
  * @return Vector [log P(a, exp(log_x)), log P(a + 1, exp(log_x))].
  * Rejects a <= 0, `nan` and `inf`.
  */
vector gamma_lcdf_logx_pair(real log_x, real a) {
  if (is_inf(log_x)) {
    real log_p = gamma_lcdf_logx(log_x, a);
    return [log_p, log_p]';
  }
  real log_lead = gamma_log_lead_logx(log_x, a);
  if (exp(log_x) < a + 1) {
    real log_p_ap1 = gamma_lcdf_logx(log_x, a + 1);
    return [log_sum_exp(log_p_ap1, log_lead), log_p_ap1]';
  }
  real log_q = gamma_log_pq_logx(log_x, a)[2];
  return [log1m_exp(log_q), log1m_exp(log_sum_exp(log_q, log_lead))]';
}

/**
  * Compute the log of a Gamma rate
  * @ingroup delay_log_cdfs
  *
  * @param rate Rate parameter of the Gamma distribution
  *
  * @return log(rate). Rejects a rate that is not positive and finite.
  */
real gamma_log_rate(real rate) {
  if (!(rate > 0) || is_inf(rate)) {
    reject("Gamma rate must be positive finite, found ", rate);
  }
  return log(rate);
}

/**
  * Compute the log CDF of the generalised gamma distribution
  * @ingroup delay_log_cdfs
  *
  * Uses the Stacy parameterisation of `flexsurv::pgengamma.orig()` in R.
  * The CDF is the regularised lower incomplete gamma function
  * P(k, (y / scale)^shape), so the Gamma (shape = 1) and Weibull (k = 1)
  * distributions are special cases.
  *
  * @param y Value at which to evaluate the log CDF (y >= 0)
  * @param shape Shape (power) parameter
  * @param scale Scale parameter
  * @param k Shape parameter of the underlying Gamma distribution
  *
  * @return Log CDF of the generalised gamma distribution
  */
real gengamma_lcdf(real y, real shape, real scale, real k) {
  if (y < 0) {
    reject("gengamma_lcdf: y must be non-negative, found y = ", y);
  }
  if (y == 0) {
    return negative_infinity();
  }
  return gamma_lcdf_logx(shape * (log(y) - log(scale)), k);
}

/**
  * Test whether a delay distribution has support only on the non-negative reals
  * @ingroup delay_log_cdfs
  *
  * Used internally to decide whether to short-circuit `dist_lcdf` at
  * `delay <= 0` and whether the ODE / nested CDF calls need to integrate over
  * negative arguments. Returns 1 for distributions with strictly non-negative
  * support, 0 otherwise. IDs match `pcd_distributions$stan_id` in R.
  *
  * @param dist_id Distribution identifier
  * @return 1 if the delay distribution has non-negative support, 0 otherwise.
  */
int dist_has_positive_support(data int dist_id) {
  if (dist_id == 1) return 1;   // Lognormal
  if (dist_id == 2) return 1;   // Gamma
  if (dist_id == 3) return 1;   // Weibull
  if (dist_id == 4) return 1;   // Exponential
  if (dist_id == 5) return 1;   // Generalised gamma
  if (dist_id == 9) return 1;   // Beta (support on [0, 1])
  if (dist_id == 13) return 1;  // Chi-square
  if (dist_id == 16) return 1;  // Inverse Gamma
  if (dist_id == 19) return 1;  // Inverse Chi-square
  if (dist_id == 21) return 1;  // Pareto
  if (dist_id == 22) return 1;  // Scaled inverse Chi-square
  if (dist_id == 29) return 1;  // Log-logistic
  return 0;
}

/**
  * Test whether `lognormal_lcdf` underflows to `-inf` at these arguments
  * @ingroup delay_log_cdfs
  *
  * Underflow makes the autodiff partial `0 / 0`, and Stan's reverse pass
  * chains that `NaN` into `mu` and `sigma` whatever weight the term is later
  * given. Callers must therefore test this before calling `lognormal_lcdf`,
  * rather than checking its result.
  *
  * The threshold is -38 on the standardised scale `(log(y) - mu) / sigma`,
  * inside the region where the CDF is still representable: `log F(y)` is
  * below -726 there, so a term dropped on this test cannot change a result
  * at double precision.
  *
  * @param y Value at which the log CDF would be evaluated
  * @param mu Location parameter on the log scale
  * @param sigma Scale parameter on the log scale
  *
  * @return 1 if `lognormal_lcdf` would underflow or `y` is non-positive,
  *   0 otherwise
  */
int lognormal_lcdf_underflows(real y, real mu, real sigma) {
  if (y <= 0) {
    return 1;
  }
  return (log(y) - mu) / sigma < -38 ? 1 : 0;
}

/**
  * Compute the log CDF of the delay distribution
  * @ingroup delay_log_cdfs
  *
  * @param delay Time delay
  * @param params Distribution parameters
  * @param dist_id Distribution identifier matching pcd_distributions in R:
  *   1: Lognormal, 2: Gamma, 3: Weibull, 4: Exponential,
  *   5: Generalised gamma, 9: Beta, 12: Cauchy, 13: Chi-square,
  *   15: Gumbel, 16: Inverse Gamma, 17: Logistic,
  *   18: Normal, 19: Inverse Chi-square,
  *   20: Double Exponential, 21: Pareto,
  *   22: Scaled Inverse Chi-square, 23: Student's t,
  *   24: Uniform, 25: von Mises,
  *   26: Non-parametric step (params = [boundaries (K+1), pmf (K)],
  *       length 2*K + 1),
  *   27/28: Non-parametric discrete hazard (params = [boundaries (K+1),
  *       hazards (K)], length 2*K + 1; hazards[K] must equal 1). 27 and
  *       28 share this likelihood and only differ in the prior on the
  *       hazards (random walk for 27, IID random effect for 28),
  *   29: Log-logistic (params = [scale, shape], the order of the Stan
  *       built-in `loglogistic_cdf`),
  *
  * @return Log CDF of the delay distribution
  *
  * @code
  * // Example: Lognormal distribution
  * real delay = 5.0;
  * array[2] real params = {0.0, 1.0}; // mean and standard deviation on log scale
  * int dist_id = 1; // Lognormal
  * real log_cdf = dist_lcdf(delay, params, dist_id);
  * @endcode
  */
real dist_lcdf(real delay, array[] real params, int dist_id) {
  if (dist_has_positive_support(dist_id) && delay <= 0) {
    return negative_infinity();
  }

  // IDs match pcd_distributions$stan_id in R
  // Guarded so a lower-tail underflow cannot put a NaN partial on the tape.
  // The downstream `exp(-inf)` differentiates to 0.
  if (dist_id == 1) {
    return lognormal_lcdf_underflows(delay, params[1], params[2])
           ? negative_infinity()
           : lognormal_lcdf(delay | params[1], params[2]);
  }
  else if (dist_id == 2) {
    return gamma_lcdf_logx(log(delay) + gamma_log_rate(params[2]), params[1]);
  }
  else if (dist_id == 3) return weibull_lcdf(delay | params[1], params[2]);
  else if (dist_id == 4) return exponential_lcdf(delay | params[1]);
  else if (dist_id == 5) return gengamma_lcdf(delay | params[1], params[2], params[3]);
  else if (dist_id == 9) return beta_lcdf(delay | params[1], params[2]);
  else if (dist_id == 12) return cauchy_lcdf(delay | params[1], params[2]);
  else if (dist_id == 13) return chi_square_lcdf(delay | params[1]);
  else if (dist_id == 15) return gumbel_lcdf(delay | params[1], params[2]);
  else if (dist_id == 16) return inv_gamma_lcdf(delay | params[1], params[2]);
  else if (dist_id == 17) return logistic_lcdf(delay | params[1], params[2]);
  else if (dist_id == 18) return normal_lcdf(delay | params[1], params[2]);
  else if (dist_id == 19) return inv_chi_square_lcdf(delay | params[1]);
  else if (dist_id == 20) return double_exponential_lcdf(delay | params[1], params[2]);
  else if (dist_id == 21) return pareto_lcdf(delay | params[1], params[2]);
  else if (dist_id == 22) return scaled_inv_chi_square_lcdf(delay | params[1], params[2]);
  else if (dist_id == 23) return student_t_lcdf(delay | params[1], params[2], params[3]);
  else if (dist_id == 24) return uniform_lcdf(delay | params[1], params[2]);
  else if (dist_id == 25) return von_mises_lcdf(delay | params[1], params[2]);
  else if (dist_id == 29) {
    return primarycensored_loglogistic_lcdf(delay | params[1], params[2]);
  }
  else if (dist_id == 26) {
    // Non-parametric step: params = [boundaries (K+1), pmf (K)].
    int K = (size(params) - 1) %/% 2;
    return pstep_lcdf(
      delay | to_vector(segment(params, 1, K + 1)),
              to_vector(segment(params, K + 2, K))
    );
  }
  else if (dist_id == 27 || dist_id == 28) {
    // Non-parametric discrete hazard: params = [boundaries (K+1),
    // hazards (K)] with hazards[K] = 1. RW (27) and RE (28) share the
    // same likelihood; they only differ in the prior.
    int K = (size(params) - 1) %/% 2;
    return phazard_lcdf(
      delay | to_vector(segment(params, 1, K + 1)),
              to_vector(segment(params, K + 2, K))
    );
  }
  else reject("Invalid distribution identifier: ", dist_id);
}

/**
  * Log CDF of the primary distribution on [0, pwindow]
  * @ingroup primary_distribution_log_cdfs
  *
  * Returns log F_primary(p) for the primary event time p in [0, pwindow].
  * Only primary_id values supported by `check_for_analytical` should be
  * passed here. The Stan `_lcdf` convention requires the `|` syntax at
  * call sites.
  *
  * @param p Primary event time in [0, pwindow]
  * @param primary_id Primary distribution identifier (1=uniform, 2=expgrowth,
  *   3=truncated logistic)
  * @param primary_params Distribution parameters (empty for uniform;
  *   [r] for expgrowth; [location, scale] for the truncated logistic)
  * @param pwindow Primary event window width
  *
  * @return log(F_primary(p))
  */
real primary_lcdf(real p, int primary_id, array[] real primary_params,
                  data real pwindow) {
  if (primary_id == 1) {
    // Uniform on [0, pwindow]: built-in uniform_lcdf matches the package
    // primary semantics over [0, pwindow].
    if (p <= 0) return negative_infinity();
    if (p >= pwindow) return 0;
    return uniform_lcdf(p | 0, pwindow);
  } else if (primary_id == 2) {
    return expgrowth_lcdf(p | 0, pwindow, primary_params[1]);
  } else if (primary_id == 3) {
    return tlogis_lcdf(p | 0, pwindow, primary_params[1], primary_params[2]);
  }
  reject("primary_lcdf: unsupported primary_id ", primary_id);
}

/**
  * Compute the log PDF of the primary distribution
  * @ingroup primary_distribution_log_pdfs
  *
  * @param x Value
  * @param primary_id Primary distribution identifier
  * @param params Distribution parameters
  * @param xmin Minimum value
  * @param xmax Maximum value
  *
  * @return Log PDF of the primary distribution
  *
  * @code
  * // Example: Uniform distribution
  * real x = 0.5;
  * int primary_id = 1; // Uniform
  * array[0] real params = {}; // No additional parameters for uniform
  * real xmin = 0;
  * real xmax = 1;
  * real log_pdf = primary_lpdf(x | primary_id, params, xmin, xmax);
  * @endcode
  */
real primary_lpdf(real x, int primary_id, array[] real params, real xmin, real xmax) {
  // Implement switch for different primary distributions
  if (primary_id == 1) return uniform_lpdf(x | xmin, xmax);
  if (primary_id == 2) return expgrowth_lpdf(x | xmin, xmax, params[1]);
  if (primary_id == 3) {
    return tlogis_lpdf(x | xmin, xmax, params[1], params[2]);
  }
  // Add more primary distributions as needed
  reject("Invalid primary distribution identifier");
}

/**
  * ODE system for the primary censored distribution
  * @ingroup ode
  *
  * @param t Time
  * @param y State variables
  * @param theta Parameters
  * @param x_r Real data
  * @param x_i Integer data
  *
  * @return Derivatives of the state variables
  */
vector primarycensored_ode(real t, vector y, array[] real theta,
                            array[] real x_r, array[] int x_i) {
  real d = x_r[1];
  int dist_id = x_i[1];
  int primary_id = x_i[2];
  real pwindow = x_r[2];
  int dist_params_len = x_i[3];
  int primary_params_len = x_i[4];

  // Extract distribution parameters
  array[dist_params_len] real params;
  if (dist_params_len) {
    params = theta[1:dist_params_len];
  }
  array[primary_params_len] real primary_params;
  if (primary_params_len) {
    int primary_loc = num_elements(theta);
    primary_params = theta[primary_loc - primary_params_len + 1:primary_loc];
  }

  real log_cdf = dist_lcdf(t | params, dist_id);
  real log_primary_pdf = primary_lpdf(d - t | primary_id, primary_params, 0, pwindow);

  return rep_vector(exp(log_cdf + log_primary_pdf), 1);
}
