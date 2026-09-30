/**
  * Compute the log of the sum in the incomplete gamma series, from log x
  * @ingroup delay_log_cdfs
  *
  * The series
  *   P(a, x) = x^a exp(-x) / Gamma(a + 1) * S,
  *   S = 1 + x / (a + 1) + x^2 / ((a + 1) (a + 2)) + ...
  * has only positive terms, so there is no cancellation. Successive terms
  * shrink by a factor of at most x / (a + n) and the sum converges for any
  * x. It takes under 60 terms for x < (a + 1) / 2, under 400 for
  * x < 0.9 (a + 1) and about 9 sqrt(a) terms near x = a + 1. Autodiff
  * differentiates the sum directly, so gradients are as accurate as the
  * value.
  *
  * @param log_x Log of the argument, log(x) with x > 0
  * @param a Shape parameter of the Gamma distribution (a > 0)
  *
  * @return log S, or `nan` if the series has not converged after 100000
  * terms (a above about 1e8 near x = a + 1)
  */
real gamma_lseries_sum_logx(real log_x, real a) {
  real x = exp(log_x);
  real term = 1;
  real total = 1;
  int n = 0;
  while (n < 100000) {
    n += 1;
    term *= x / (a + n);
    total += term;
    if (term < 1e-17 * total) {
      return log(total);
    }
  }
  return not_a_number();
}

/**
  * Compute log Q(a, x) by the incomplete gamma continued fraction
  * @ingroup delay_log_cdfs
  *
  * Q(a, x) = 1 - P(a, x) is the regularised upper incomplete gamma
  * function. The continued fraction
  *   Gamma(a, x) = exp(-x) x^a / (x + 1 - a - 1 (1 - a) / (x + 3 - a -
  *   2 (2 - a) / (x + 5 - a - ...)))
  * is evaluated with the modified Lentz algorithm. It converges for
  * x >= a + 1, in about 10 to 20 steps for x > 2 a and about
  * 1.3 sqrt(a) + 10 steps near x = a + 1 for a >= 10. Autodiff
  * differentiates the algorithm, so gradients are as accurate as the value.
  *
  * The numerator of step i is i (i - a), which is zero at i = a when `a` is
  * an integer. The value of the fraction is then final, but its derivative
  * with respect to `a` is not, because the derivative of the numerator is
  * not zero. Stopping on the change in the value therefore leaves the
  * derivative wrong when the fraction needs more than a steps. That is the
  * case for a of about 5 or less near x = a + 1, where the derivative is
  * wrong by 1e-7 to 1e-2 at integer `a`, so this is only called for
  * a >= 10 by `gamma_lcdf_logx()`. For a >= 10 the derivative at integer
  * `a` is accurate to 1e-8.
  *
  * @param log_x Log of the argument, log(x) with x > 0
  * @param a Shape parameter of the Gamma distribution (a > 0)
  *
  * @return log Q(a, exp(log_x)), or `nan` if the fraction has not
  * converged after 100000 steps (a above about 1e8 near x = a + 1)
  */
real gamma_lccdf_cf_logx(real log_x, real a) {
  real x = exp(log_x);
  real tiny = 1e-300;
  real b = x + 1 - a;
  real c = 1 / tiny;
  real d = 1 / b;
  real h = d;
  int i = 0;
  while (i < 100000) {
    i += 1;
    real numerator = -i * (i - a);
    b += 2;
    d = numerator * d + b;
    if (abs(d) < tiny) d = tiny;
    c = b + numerator / c;
    if (abs(c) < tiny) c = tiny;
    d = 1 / d;
    real delta = d * c;
    h *= delta;
    if (abs(delta - 1) < 1e-15) {
      return a * log_x - x - lgamma(a) + log(h);
    }
  }
  return not_a_number();
}

/**
  * Compute the log CDF of a unit rate Gamma distribution from the log of x
  * @ingroup delay_log_cdfs
  *
  * Returns log P(a, x), the log of the regularised lower incomplete gamma
  * function, for x = exp(log_x). Stan's `gamma_lcdf` is accurate in value
  * but not in its partial derivative with respect to `a`, which fails in
  * three regions (see #363 and #381).
  *
  * - Lower tail. `gamma_lcdf` underflows to `-inf` deep in the lower tail,
  *   where the true value is finite. The gradient is inaccurate before that,
  *   by 1e-2 to 0.9 once log P is below about -10 and x < 0.9 (a + 1), and
  *   can be `nan` or 0 for a of about 1000 or more.
  * - Large shape. The gradient is `nan` or throws "n (internal counter)
  *   exceeded 100000 iterations" for a of about 700 or more, anywhere in
  *   the body or the tail.
  * - Upper tail. The gradient is off by about 1e-14, and is `nan` once
  *   Q(a, x) underflows, for a of about 10 or more.
  *
  * The rules are
  * 1. `a >= 10`. Use `gamma_lseries_sum_logx()` for x < a + 1 and
  *    `log1m_exp(gamma_lccdf_cf_logx())` for x >= a + 1. Both are exact.
  * 2. `a < 10` and x < 0.9 (a + 1) with the leading term
  *    x^a exp(-x) / Gamma(a + 1) below exp(-10). Use
  *    `gamma_lseries_sum_logx()`.
  * 3. Otherwise `gamma_lcdf`, which is accurate in the rest of the region.
  *
  * The value changes by under 1e-10 relative to `gamma_lcdf` where the
  * rule changes, because every branch is accurate. If the series or the
  * fraction does not converge the result is from `gamma_lcdf`, which either
  * has the right value or throws.
  *
  * Taking `log_x` rather than x keeps the result finite when x itself would
  * underflow, as it does for a generalised gamma with a large `shape`.
  *
  * Measured against `pgamma(log.p = TRUE)` for a from 0.01 to 1e6, and x
  * from 0.01 a to 30 a. The log CDF has a relative error of 1e-9 or below,
  * except for a < 10 where x > 5 a puts log P below 1e-12 in size, and
  * `gamma_lcdf` returns an absolute error of 1e-16 in it. The derivative
  * with respect to `a` has a relative error of 5e-9 or below for a up to
  * 3e4, 1.4e-6 at 1e5 and 2e-8 at 1e6, ignoring derivatives below 1e-13 in
  * size. The error comes from the leading term, which is formed from terms
  * of size a log(x).
  *
  * @param log_x Log of the argument, log(x) with x > 0
  * @param a Shape parameter of the Gamma distribution (a > 0)
  *
  * @return log P(a, exp(log_x)), `-inf` when `log_x` is `-inf` and 0 when
  * it is `inf`
  */
real gamma_lcdf_logx(real log_x, real a) {
  if (log_x == negative_infinity()) {
    return negative_infinity();
  }
  if (log_x == positive_infinity()) {
    return 0;
  }
  real x = exp(log_x);
  real result = not_a_number();
  if (a >= 10) {
    result = x < a + 1
             ? a * log_x - x - lgamma(a + 1)
               + gamma_lseries_sum_logx(log_x, a)
             : log1m_exp(gamma_lccdf_cf_logx(log_x, a));
  } else if (x < 0.9 * (a + 1)
             && a * log_x - x - lgamma(a + 1) < -10) {
    result = a * log_x - x - lgamma(a + 1)
             + gamma_lseries_sum_logx(log_x, a);
  } else {
    return gamma_lcdf(x | a, 1);
  }
  if (is_nan(result)) {
    return gamma_lcdf(x | a, 1);
  }
  return result;
}

/**
  * Compute the log CDFs of Gamma(a) and Gamma(a + 1) from the log of x
  * @ingroup delay_log_cdfs
  *
  * Returns [log P(a, x), log P(a + 1, x)] with unit rate, as needed by the
  * uniform primary solution for a Gamma delay. The two are related by
  *   P(a + 1, x) = P(a, x) - x^a exp(-x) / Gamma(a + 1),
  * which subtracts two equal terms when x is small compared to a + 1. It is
  * applied only where that is safe, and otherwise both values come from the
  * same series or fraction as `gamma_lcdf_logx()`, so that the pair costs
  * one evaluation and not two.
  *
  * - Series region, `a >= 10` and x < a + 1, or `a < 10` and x < 0.5 (a + 1)
  *   or x < 0.9 (a + 1) with the leading term below exp(-10). The sum S for
  *   `a + 1` is evaluated, and the sum for `a` is 1 + x S / (a + 1) with
  *   no subtraction.
  * - Fraction region, `a >= 10` and x >= a + 1. Q(a + 1, x) is
  *   Q(a, x) + x^a exp(-x) / Gamma(a + 1), a sum of positive terms.
  * - Otherwise `a < 10` and x is at least about half of a + 1, so that the
  *   subtraction loses at most about one bit. P(a, x) is from `gamma_lcdf`.
  *
  * If the series or the fraction does not converge, both values are from
  * `gamma_lcdf_logx()`, which falls back to `gamma_lcdf`.
  *
  * @param log_x Log of the argument, log(x) with x > 0
  * @param a Shape parameter of the Gamma distribution (a > 0)
  *
  * @return Vector [log P(a, exp(log_x)), log P(a + 1, exp(log_x))]
  */
vector gamma_lcdf_logx_pair(real log_x, real a) {
  if (log_x == negative_infinity()) {
    return rep_vector(negative_infinity(), 2);
  }
  if (log_x == positive_infinity()) {
    return rep_vector(0, 2);
  }
  real x = exp(log_x);
  real log_lead = a * log_x - x - lgamma(a + 1);
  vector[2] result = rep_vector(not_a_number(), 2);
  if (x < a + 1
      && (a >= 10 || x < 0.5 * (a + 1)
          || (x < 0.9 * (a + 1) && log_lead < -10))) {
    // S for a + 1, from which S for a follows without subtraction
    real log_sum_kp1 = gamma_lseries_sum_logx(log_x, a + 1);
    result[1] = log_lead
                + log1p_exp(log_x - log(a + 1) + log_sum_kp1);
    result[2] = log_lead + log_x - log(a + 1) + log_sum_kp1;
  } else if (a >= 10) {
    real log_q = gamma_lccdf_cf_logx(log_x, a);
    result[1] = log1m_exp(log_q);
    result[2] = log1m_exp(log_sum_exp(log_q, log_lead));
  } else {
    result[1] = gamma_lcdf(x | a, 1);
    result[2] = log_diff_exp(result[1], log_lead);
  }
  if (is_nan(result[1]) || is_nan(result[2])) {
    return [gamma_lcdf_logx(log_x, a), gamma_lcdf_logx(log_x, a + 1)]';
  }
  return result;
}

/**
  * Compute the log CDF of the generalised gamma distribution
  * @ingroup delay_log_cdfs
  *
  * Uses the Stacy parameterisation of `flexsurv::pgengamma.orig()` in R.
  * The CDF is the regularised lower incomplete gamma function
  * P(k, (y / scale)^shape), so the Gamma (shape = 1) and Weibull (k = 1)
  * distributions are special cases. Uses `gamma_lcdf_logx()` for
  * lower-tail accuracy.
  *
  * @param y Value at which to evaluate the log CDF (y >= 0). Negative y
  * is rejected.
  * @param shape Shape (power) parameter
  * @param scale Scale parameter
  * @param k Shape parameter of the underlying Gamma distribution
  *
  * @return Log CDF of the generalised gamma distribution, `-inf` for y = 0
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
  *       hazards (random walk for 27, IID random effect for 28).
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
    return gamma_lcdf_logx(log(delay) + log(params[2]), params[1]);
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
  * @param primary_id Primary distribution identifier (1=uniform, 2=expgrowth)
  * @param primary_params Distribution parameters (empty for uniform;
  *   [r] for expgrowth)
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
