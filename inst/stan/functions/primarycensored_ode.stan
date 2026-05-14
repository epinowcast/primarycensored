/**
  * Compute the log CDF of the delay distribution
  * @ingroup delay_log_cdfs
  *
  * @param delay Time delay
  * @param params Distribution parameters
  * @param dist_id Distribution identifier matching pcd_distributions in R:
  *   1: Lognormal, 2: Gamma, 3: Weibull, 4: Exponential,
  *   9: Beta, 12: Cauchy, 13: Chi-square,
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
  if (dist_id == 9) return 1;   // Beta (support on [0, 1])
  if (dist_id == 13) return 1;  // Chi-square
  if (dist_id == 16) return 1;  // Inverse Gamma
  if (dist_id == 19) return 1;  // Inverse Chi-square
  if (dist_id == 21) return 1;  // Pareto
  if (dist_id == 22) return 1;  // Scaled inverse Chi-square
  return 0;
}

real dist_lcdf(real delay, array[] real params, int dist_id) {
  if (dist_has_positive_support(dist_id) && delay <= 0) {
    return negative_infinity();
  }

  // IDs match pcd_distributions$stan_id in R
  if (dist_id == 1) return lognormal_lcdf(delay | params[1], params[2]);
  else if (dist_id == 2) return gamma_lcdf(delay | params[1], params[2]);
  else if (dist_id == 3) return weibull_lcdf(delay | params[1], params[2]);
  else if (dist_id == 4) return exponential_lcdf(delay | params[1]);
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
