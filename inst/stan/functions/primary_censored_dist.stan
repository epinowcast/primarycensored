/**
 * Primary event censored distribution functions
 */

/**
  * Compute the log CDF of the delay distribution
  *
  * @param delay Time delay
  * @param params Distribution parameters
  * @param dist_id Distribution identifier
  *   1: Lognormal, 2: Gamma, 3: Normal, 4: Exponential, 5: Weibull,
  *   6: Beta, 7: Cauchy, 8: Chi-square, 9: Inverse Chi-square,
  *   10: Double Exponential, 11: Inverse Gamma, 12: Logistic,
  *   13: Pareto, 14: Scaled Inverse Chi-square, 15: Student's t,
  *   16: Uniform, 17: von Mises
  *
  * @return Log CDF of the delay distribution
  *
  * @code
  * // Example: Lognormal distribution
  * real delay = 5.0;
  * array[2] real params = {0.0, 1.0}; // mean and standard deviation on log scale
  * int dist_id = 1; // Lognormal
  * real log_cdf = dist_lcdf(delay, params, dist_id);
  */
  real dist_lcdf(real delay, array[] real params, int dist_id) {
    if (delay <= 0) return negative_infinity();

    // Use if-else statements to handle different distribution types
    if (dist_id == 1) return lognormal_lcdf(delay | params[1], params[2]);
    else if (dist_id == 2) return gamma_lcdf(delay | params[1], params[2]);
    else if (dist_id == 3) return normal_lcdf(delay | params[1], params[2]);
    else if (dist_id == 4) return exponential_lcdf(delay | params[1]);
    else if (dist_id == 5) return weibull_lcdf(delay | params[1], params[2]);
    else if (dist_id == 6) return beta_lcdf(delay | params[1], params[2]);
    else if (dist_id == 7) return cauchy_lcdf(delay | params[1], params[2]);
    else if (dist_id == 8) return chi_square_lcdf(delay | params[1]);
    else if (dist_id == 9) return inv_chi_square_lcdf(delay | params[1]);
    else if (dist_id == 10) return double_exponential_lcdf(delay | params[1], params[2]);
    else if (dist_id == 11) return inv_gamma_lcdf(delay | params[1], params[2]);
    else if (dist_id == 12) return logistic_lcdf(delay | params[1], params[2]);
    else if (dist_id == 13) return pareto_lcdf(delay | params[1], params[2]);
    else if (dist_id == 14) return scaled_inv_chi_square_lcdf(delay | params[1], params[2]);
    else if (dist_id == 15) return student_t_lcdf(delay | params[1], params[2], params[3]);
    else if (dist_id == 16) return uniform_lcdf(delay | params[1], params[2]);
    else if (dist_id == 17) return von_mises_lcdf(delay | params[1], params[2]);
    else reject("Invalid distribution identifier");
  }

/**
  * Compute the log PDF of the primary distribution
  *
  * @param x Value
  * @param primary_dist_id Primary distribution identifier
  * @param params Distribution parameters
  * @param min Minimum value
  * @param max Maximum value
  *
  * @return Log PDF of the primary distribution
  *
  * @code
  * // Example: Uniform distribution
  * real x = 0.5;
  * int primary_dist_id = 1; // Uniform
  * array[0] real params = {}; // No additional parameters for uniform
  * real min = 0;
  * real max = 1;
  * real log_pdf = primary_dist_lpdf(x, primary_dist_id, params, min, max);
  */
real primary_dist_lpdf(real x, int primary_dist_id, array[] real params, real min, real max) {
  // Implement switch for different primary distributions
  if (primary_dist_id == 1) return uniform_lpdf(x | min, max);
  if (primary_dist_id == 2) return expgrowth_lpdf(x | min, max, params[1]);
  // Add more primary distributions as needed
  reject("Invalid primary distribution identifier");
}

/**
  * Compute the integrand for the primary censored distribution
  *
  * @param p Integration variable
  * @param xc Delay
  * @param theta Distribution parameters
  * @param x_r Real data (contains pwindow)
  * @param x_i Integer data (contains dist_id and primary_dist_id)
  *
  * @return Value of the integrand
  *
  * @code
  * // Example: Lognormal delay distribution with uniform primary distribution
  * real p = 0.5;
  * real xc = 2.0;
  * array[2] real theta = {0.0, 1.0}; // mean and standard deviation on log scale
  * array[1] real x_r = {1.0}; // pwindow
  * array[2] int x_i = {1, 1}; // dist_id = 1 (Lognormal), primary_dist_id = 1 (Uniform)
  * real integrand_value = primary_censored_integrand(p, xc, theta, x_r, x_i);
  */
real primary_censored_integrand(real p, real xc, array[] real theta,
                                data array[] real x_r, data array[] int x_i) {
  real d = xc;
  real pwindow = x_r[1];
  int dist_id = x_i[1];
  int primary_dist_id = x_i[2];
  real d_adj = d - p;

  real log_cdf = dist_lcdf(d_adj | theta, dist_id);
  real log_primary_pdf = primary_dist_lpdf(
    p | primary_dist_id, theta, 0, pwindow
  );

  return exp(log_cdf + log_primary_pdf);
}

/**
  * Compute the integrand for the truncated primary censored distribution
  *
  * @param p Integration variable
  * @param xc Delay
  * @param theta Distribution parameters (including D)
  * @param x_r Real data (contains pwindow)
  * @param x_i Integer data (contains dist_id and primary_dist_id)
  *
  * @return Value of the integrand
  *
  * @code
  * // Example: Truncated Gamma delay distribution with uniform primary distribution
  * real p = 0.5;
  * real xc = 2.0;
  * array[3] real theta = {2.0, 1.0, 5.0}; // shape, scale, and D (truncation point)
  * array[1] real x_r = {1.0}; // pwindow
  * array[2] int x_i = {2, 1}; // dist_id = 2 (Gamma), primary_dist_id = 1 (Uniform)
  * real integrand_value = primary_censored_integrand_truncated(p, xc, theta, x_r, x_i);
  */
real primary_censored_integrand_truncated(real p, real xc, array[] real theta,
                                          data array[] real x_r, data array[] int x_i) {
  real d = xc;
  real pwindow = x_r[1];
  int dist_id = x_i[1];
  int primary_dist_id = x_i[2];
  real d_adj = d - p;
  real D = theta[size(theta)];
  real D_adj = D - p;

  real log_cdf = dist_lcdf(d_adj | theta[1:(size(theta)-1)], dist_id);
  real log_cdf_D = dist_lcdf(D_adj | theta[1:(size(theta)-1)], dist_id);
  real log_primary_pdf = primary_dist_lpdf(
    p | primary_dist_id, theta[1:(size(theta)-1)], 0, pwindow
  );

  return exp(log_cdf - log_cdf_D + log_primary_pdf);
}

/**
  * Compute the primary event censored CDF for a single delay
  *
  * @param d Delay
  * @param dist_id Distribution identifier
  * @param params Array of distribution parameters
  * @param pwindow Primary event window
  * @param D Maximum delay (truncation point)
  * @param primary_dist_id Primary distribution identifier
  * @param primary_params Primary distribution parameters
  *
  * @return Primary event censored CDF, normalized by D if finite (truncation adjustment)
  *
  * @code
  * // Example: Weibull delay distribution with uniform primary distribution
  * real d = 3.0;
  * int dist_id = 5; // Weibull
  * array[2] real params = {2.0, 1.5}; // shape and scale
  * real pwindow = 1.0;
  * real D = positive_infinity();
  * int primary_dist_id = 1; // Uniform
  * array[0] real primary_params = {};
  * real cdf = primary_censored_dist_cdf(d, dist_id, params, pwindow, D, primary_dist_id, primary_params);
  */
real primary_censored_dist_cdf(real d, int dist_id, array[] real params,
                               real pwindow, real D, int primary_dist_id,
                               array[] real primary_params) {
  real result;
  if (d <= 0 || d >= D) {
    return 0;
  }

  if (is_inf(D)) {
    result = integrate_1d(
      primary_censored_integrand,
      0, pwindow, {d}, params, {pwindow}, {dist_id, primary_dist_id}
    );
  } else {
    result = integrate_1d(
      primary_censored_integrand_truncated,
      0, pwindow, {d}, append_array(params, {D}), {pwindow}, {dist_id, primary_dist_id}
    );
  }

  return result;
}

/**
  * Compute the primary event censored log CDF for a single delay
  *
  * @param d Delay
  * @param dist_id Distribution identifier
  * @param params Array of distribution parameters
  * @param pwindow Primary event window
  * @param D Maximum delay (truncation point)
  * @param primary_dist_id Primary distribution identifier
  * @param primary_params Primary distribution parameters
  *
  * @return Primary event censored log CDF, normalized by D if finite (truncation adjustment)
  *
  * @code
  * // Example: Exponential delay distribution with uniform primary distribution
  * real d = 2.0;
  * int dist_id = 4; // Exponential
  * array[1] real params = {0.5}; // rate
  * real pwindow = 1.0;
  * real D = 10.0;
  * int primary_dist_id = 1; // Uniform
  * array[0] real primary_params = {};
  * real log_cdf = primary_censored_dist_lcdf(d, dist_id, params, pwindow, D, primary_dist_id, primary_params);
  */
real primary_censored_dist_lcdf(real d, int dist_id, array[] real params,
                                real pwindow, real D, int primary_dist_id,
                                array[] real primary_params) {
  if (d <= 0 || d >= D) {
    return negative_infinity();
  }
  return log(primary_censored_dist_cdf(d, dist_id, params, pwindow, D, primary_dist_id, primary_params));
}

/**
  * Compute the primary event censored log PMF for a single delay
  *
  * @param d Delay
  * @param dist_id Distribution identifier
  * @param params Array of distribution parameters
  * @param pwindow Primary event window
  * @param swindow Secondary event window
  * @param D Maximum delay (truncation point)
  * @param primary_dist_id Primary distribution identifier
  * @param primary_params Primary distribution parameters
  *
  * @return Primary event censored log PMF, normalized by D if finite (truncation adjustment)
  *
  * @code
  * // Example: Normal delay distribution with uniform primary distribution
  * real d = 1.5;
  * int dist_id = 3; // Normal
  * array[2] real params = {0.0, 1.0}; // mean and standard deviation
  * real pwindow = 1.0;
  * real swindow = 0.1;
  * real D = positive_infinity();
  * int primary_dist_id = 1; // Uniform
  * array[0] real primary_params = {};
  * real log_pmf = primary_censored_dist_lpmf(d, dist_id, params, pwindow, swindow, D, primary_dist_id, primary_params);
  */
real primary_censored_dist_lpmf(real d, int dist_id, array[] real params,
                                real pwindow, real swindow, real D,
                                int primary_dist_id, array[] real primary_params) {
  real log_cdf_upper = primary_censored_dist_lcdf(d + swindow, dist_id, params,
                        pwindow, D, primary_dist_id, primary_params);
  real log_cdf_lower = primary_censored_dist_lcdf(d, dist_id, params,
                        pwindow, D, primary_dist_id, primary_params);
  return log_diff_exp(log_cdf_upper, log_cdf_lower);
}

/**
  * Compute the primary event censored PMF for a single delay
  *
  * @param d Delay
  * @param dist_id Distribution identifier
  * @param params Array of distribution parameters
  * @param pwindow Primary event window
  * @param swindow Secondary event window
  * @param D Maximum delay (truncation point)
  * @param primary_dist_id Primary distribution identifier
  * @param primary_params Primary distribution parameters
  *
  * @return Primary event censored PMF, normalized by D if finite (truncation adjustment)
  *
  * @code
  * // Example: Gamma delay distribution with exponential growth primary distribution
  * real d = 2.5;
  * int dist_id = 2; // Gamma
  * array[2] real params = {2.0, 1.0}; // shape and rate
  * real pwindow = 1.0;
  * real swindow = 0.1;
  * real D = 10.0;
  * int primary_dist_id = 2; // Exponential growth
  * array[1] real primary_params = {0.5}; // growth rate
  * real pmf = primary_censored_dist_pmf(d, dist_id, params, pwindow, swindow, D, primary_dist_id, primary_params);
  */
real primary_censored_dist_pmf(real d, int dist_id, array[] real params,
                               real pwindow, real swindow, real D,
                               int primary_dist_id, array[] real primary_params) {
  return exp(primary_censored_dist_lpmf(d, dist_id, params, pwindow, swindow, D, primary_dist_id, primary_params));
}

/**
  * Compute the primary event censored log PMF for integer delays up to D
  *
  * @param D Maximum delay (truncation point)
  * @param dist_id Distribution identifier
  * @param params Array of distribution parameters
  * @param pwindow Primary event window
  * @param primary_dist_id Primary distribution identifier
  * @param primary_params Primary distribution parameters
  *
  * @return Vector of primary event censored log PMFs for integer delays 0 to D-1
  *
  * This function differs from primary_censored_dist_lpmf in that it:
  * 1. Computes PMFs for all integer delays from 0 to D-1 in one call
  * 2. Assumes integer delays (swindow = 1)
  * 3. Is more computationally efficient for multiple delay calculations
  *
  * @code
  * // Example: Normal delay distribution with uniform primary distribution
  * int D = 10;
  * int dist_id = 3; // Normal
  * array[2] real params = {0.0, 1.0}; // mean and standard deviation
  * real pwindow = 7.0;
  * int primary_dist_id = 1; // Uniform
  * array[0] real primary_params = {};
  * vector[D] log_pmf = primary_censored_sint_lpmf(D, dist_id, params, pwindow, primary_dist_id, primary_params);
  */
vector primary_censored_sint_lpmf(int D, int dist_id, array[] real params,
                                  real pwindow, int primary_dist_id, array[] real primary_params) {
  vector[D] log_cdfs;
  vector[D] log_pmfs;
  real log_normalizer;

  // Compute log CDFs
  for (d in 1:D) {
    log_cdfs[d] = primary_censored_dist_lcdf(d, dist_id, params, pwindow, positive_infinity(), primary_dist_id, primary_params);
  }

  // Compute log normalizer
  log_normalizer = log_cdfs[D];

  // Compute log PMFs
  log_pmfs[1] = log_cdfs[1] - log_normalizer;
  for (d in 2:D) {
    log_pmfs[d] = log_diff_exp(log_cdfs[d], log_cdfs[d-1]) - log_normalizer;
  }

  return log_pmfs;
}

/**
* Compute the primary event censored PMF for integer delays up to D
*
* @param D Maximum delay (truncation point)
* @param dist_id Distribution identifier
* @param params Array of distribution parameters
* @param pwindow Primary event window
* @param primary_dist_id Primary distribution identifier
* @param primary_params Primary distribution parameters
*
* @return Vector of primary event censored PMFs for integer delays 0 to D-1
*
* This function differs from primary_censored_dist_pmf in that it:
* 1. Computes PMFs for all integer delays from 0 to D-1 in one call
* 2. Assumes integer delays (swindow = 1)
* 3. Is more computationally efficient for multiple delay calculations
*
* @code
* // Example: Gamma delay distribution with exponential growth primary distribution
* int D = 10;
* int dist_id = 2; // Gamma
* array[2] real params = {2.0, 1.0}; // shape and rate
* real pwindow = 7.0;
* int primary_dist_id = 2; // Exponential growth
* array[1] real primary_params = {0.5}; // growth rate
* vector[D] pmf = primary_censored_sint_pmf(D, dist_id, params, pwindow, primary_dist_id, primary_params);
*/
vector primary_censored_sint_pmf(int D, int dist_id, array[] real params,
                                 real pwindow, int primary_dist_id, array[] real primary_params) {
  return exp(primary_censored_sint_lpmf(D, dist_id, params, pwindow, primary_dist_id, primary_params));
}
