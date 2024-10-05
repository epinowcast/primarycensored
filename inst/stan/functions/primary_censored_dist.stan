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
  * @endcode
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
  * @endcode
  */
real primary_dist_lpdf(real x, int primary_dist_id, array[] real params, real min, real max) {
  // Implement switch for different primary distributions
  if (primary_dist_id == 1) return uniform_lpdf(x | min, max);
  if (primary_dist_id == 2) return expgrowth_lpdf(x | min, max, params[1]);
  // Add more primary distributions as needed
  reject("Invalid primary distribution identifier");
}

/**
  * ODE system for the primary censored distribution
  *
  * @param t Time
  * @param y State variables
  * @param theta Parameters
  * @param x_r Real data
  * @param x_i Integer data
  *
  * @return Derivatives of the state variables
  */
vector primary_censored_ode(real t, vector y, array[] real theta,
                            array[] real x_r, array[] int x_i) {
  real d = x_r[1];
  int dist_id = x_i[1];
  int primary_dist_id = x_i[2];
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
    int primary_loc = size(theta);
    primary_params = theta[primary_loc - primary_params_len + 1:primary_loc];
  }

  real log_cdf = dist_lcdf(t | params, dist_id);
  real log_primary_pdf = primary_dist_lpdf(d - t | primary_dist_id, primary_params, 0, pwindow);

  return rep_vector(exp(log_cdf + log_primary_pdf), 1);
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
  */
real primary_censored_dist_cdf(data real d, int dist_id, array[] real params,
                               data real pwindow, data real D,
                               int primary_dist_id,
                               array[] real primary_params) {
  real result;
  if (d <= 0) {
    return 0;
  }

  if (d >= D) {
    return 1;
  }

  // Check if an analytical solution exists
  if (check_for_analytical(dist_id, primary_dist_id)) {
    // Use analytical solution
    result = primary_censored_dist_analytical_cdf(
      d | dist_id, params, pwindow, D, primary_dist_id, primary_params
    );
  } else {
    // Use numerical integration for other cases
    real lower_bound = max({d - pwindow, 1e-6});
    array[size(params) + size(primary_params)] real theta = append_array(params, primary_params);
    array[4] int ids = {dist_id, primary_dist_id, size(params), size(primary_params)};

    vector[1] y0 = rep_vector(0.0, 1);
    result = ode_rk45(primary_censored_ode, y0, lower_bound, {d}, theta, {d, pwindow}, ids)[1, 1];

    if (!is_inf(D)) {
      real log_cdf_D = primary_censored_dist_lcdf(
        D | dist_id, params, pwindow, positive_infinity(), primary_dist_id,primary_params
      );
      result = exp(log(result) - log_cdf_D);
    }
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
  * // Example: Weibull delay distribution with uniform primary distribution
  * real d = 3.0;
  * int dist_id = 5; // Weibull
  * array[2] real params = {2.0, 1.5}; // shape and scale
  * real pwindow = 1.0;
  * real D = positive_infinity();
  * int primary_dist_id = 1; // Uniform
  * array[0] real primary_params = {};
  * real log_cdf = primary_censored_dist_lcdf(
  *   d, dist_id, params, pwindow, D, primary_dist_id, primary_params
  * );
  * @endcode
  */
real primary_censored_dist_lcdf(data real d, int dist_id, array[] real params,
                                data real pwindow, data real D,
                                int primary_dist_id,
                                array[] real primary_params) {
  real result;

  if (d <= 0) {
    return negative_infinity();
  }

  if (d >= D) {
    return 0;
  }

  // Check if an analytical solution exists
  if (check_for_analytical(dist_id, primary_dist_id)) {
    result = primary_censored_dist_analytical_lcdf(
      d | dist_id, params, pwindow, positive_infinity(), primary_dist_id, primary_params
    );
  } else {
    // Use numerical integration
    result = log(primary_censored_dist_cdf(
      d | dist_id, params, pwindow, positive_infinity(), primary_dist_id, primary_params
    ));
  }

  // Handle truncation
  if (!is_inf(D)) {
    real log_cdf_D = primary_censored_dist_lcdf(
      D | dist_id, params, pwindow, positive_infinity(), primary_dist_id, primary_params
    );
    result = result - log_cdf_D;
  }

  return result;
}

/**
  * Compute the primary event censored log PMF for a single delay
  *
  * @param d Delay (integer)
  * @param dist_id Distribution identifier
  * @param params Array of distribution parameters
  * @param pwindow Primary event window
  * @param d_upper Upper bound for the delay interval
  * @param D Maximum delay (truncation point)
  * @param primary_dist_id Primary distribution identifier
  * @param primary_params Primary distribution parameters
  *
  * @return Primary event censored log PMF, normalized by D if finite (truncation adjustment)
  *
  * @code
  * // Example: Weibull delay distribution with uniform primary distribution
  * int d = 3;
  * int dist_id = 5; // Weibull
  * array[2] real params = {2.0, 1.5}; // shape and scale
  * real pwindow = 1.0;
  * real d_upper = 4.0;
  * real D = positive_infinity();
  * int primary_dist_id = 1; // Uniform
  * array[0] real primary_params = {};
  * real log_pmf = primary_censored_dist_lpmf(
  *   d, dist_id, params, pwindow, d_upper, D, primary_dist_id, primary_params
  * );
  * @endcode
  */
real primary_censored_dist_lpmf(data int d, int dist_id, array[] real params,
                                data real pwindow, data real d_upper,
                                data real D, int primary_dist_id,
                                array[] real primary_params) {
  if (d_upper > D) {
    reject("Upper truncation point is greater than D. It is ", d_upper,
           " and D is ", D, ". Resolve this by increasing D to be greater or equal to d + swindow or decreasing swindow.");
  }
  if (d_upper <= d) {
    reject("Upper truncation point is less than or equal to d. It is ", d_upper,
           " and d is ", d, ". Resolve this by increasing d to be less than d_upper.");
  }
  real log_cdf_upper = primary_censored_dist_lcdf(
    d_upper | dist_id, params, pwindow, positive_infinity(), primary_dist_id, primary_params
  );
  real log_cdf_lower = primary_censored_dist_lcdf(
    d | dist_id, params, pwindow, positive_infinity(), primary_dist_id, primary_params
  );
  if (!is_inf(D)) {
    real log_cdf_D;

    if (d_upper == D) {
      log_cdf_D = log_cdf_upper;
    } else {
      log_cdf_D = primary_censored_dist_lcdf(
        D | dist_id, params, pwindow, positive_infinity(), primary_dist_id, primary_params
      );
    }
    return log_diff_exp(log_cdf_upper, log_cdf_lower) - log_cdf_D;
  } else {
    return log_diff_exp(log_cdf_upper, log_cdf_lower);
  }
}

/**
  * Compute the primary event censored PMF for a single delay
  *
  * @param d Delay (integer)
  * @param dist_id Distribution identifier
  * @param params Array of distribution parameters
  * @param pwindow Primary event window
  * @param d_upper Upper bound for the delay interval
  * @param D Maximum delay (truncation point)
  * @param primary_dist_id Primary distribution identifier
  * @param primary_params Primary distribution parameters
  *
  * @return Primary event censored PMF, normalized by D if finite (truncation adjustment)
  *
  * @code
  * // Example: Weibull delay distribution with uniform primary distribution
  * int d = 3;
  * real d = 3.0;
  * int dist_id = 5; // Weibull
  * array[2] real params = {2.0, 1.5}; // shape and scale
  * real pwindow = 1.0;
  * real swindow = 0.1;
  * real D = positive_infinity();
  * int primary_dist_id = 1; // Uniform
  * array[0] real primary_params = {};
  * real pmf = primary_censored_dist_pmf(d, dist_id, params, pwindow, swindow, D, primary_dist_id, primary_params);
  * @endcode
  */
real primary_censored_dist_pmf(data int d, int dist_id, array[] real params,
                               data real pwindow, data real d_upper,
                               data real D, int primary_dist_id,
                               array[] real primary_params) {
  return exp(
    primary_censored_dist_lpmf(
      d | dist_id, params, pwindow, d_upper, D, primary_dist_id, primary_params
    )
  );
}

/**
  * Compute the primary event censored log PMF for integer delays up to max_delay
  *
  * @param max_delay Maximum delay to compute PMF for
  * @param D Maximum delay (truncation point), must be at least max_delay + 1
  * @param dist_id Distribution identifier
  * @param params Array of distribution parameters
  * @param pwindow Primary event window
  * @param primary_dist_id Primary distribution identifier
  * @param primary_params Primary distribution parameters
  *
  * @return Vector of primary event censored log PMFs for delays \[0, 1\] to
  * \[max_delay, max_delay + 1\].
  *
  * This function differs from primary_censored_dist_lpmf in that it:
  * 1. Computes PMFs for all integer delays from \[0, 1\] to \[max_delay,
  *    max_delay + 1\] in one call.
  * 2. Assumes integer delays (swindow = 1)
  * 3. Is more computationally efficient for multiple delay calculation as it
  *    reduces the number of integration calls.
  *
  * @code
  * // Example: Weibull delay distribution with uniform primary distribution
  * int max_delay = 10;
  * real D = 15.0;
  * int dist_id = 5; // Weibull
  * array[2] real params = {2.0, 1.5}; // shape and scale
  * real pwindow = 7.0;
  * int primary_dist_id = 1; // Uniform
  * array[0] real primary_params = {};

  * vector[max_delay] log_pmf =
  *   primary_censored_sone_lpmf_vectorized(
  *      max_delay, D, dist_id, params, pwindow, primary_dist_id,
  *      primary_params
  *   );
  * @endcode
  */
vector primary_censored_sone_lpmf_vectorized(
  int max_delay, data real D, int dist_id,
  array[] real params, data real pwindow,
  int primary_dist_id, array[] real primary_params
) {

  int upper_interval = max_delay + 1;
  vector[upper_interval] log_pmfs;
  vector[upper_interval] log_cdfs;
  real log_normalizer;

  // Check if D is at least max_delay + 1
  if (D < upper_interval) {
    reject("D must be at least max_delay + 1");
  }

  // Compute log CDFs
  for (d in 1:upper_interval) {
    log_cdfs[d] = primary_censored_dist_lcdf(
      d | dist_id, params, pwindow, positive_infinity(), primary_dist_id,
      primary_params
    );
  }

  // Compute log normalizer using upper_interval
  if (D > upper_interval) {
    if (is_inf(D)) {
      log_normalizer = 0; // No normalization needed for infinite D
    } else {
      log_normalizer = primary_censored_dist_lcdf(
        D | dist_id, params, pwindow, positive_infinity(),
        primary_dist_id, primary_params
      );
    }
  } else {
    log_normalizer = log_cdfs[upper_interval];
  }

  // Compute log PMFs
  log_pmfs[1] = log_cdfs[1] - log_normalizer;
  for (d in 2:upper_interval) {
    log_pmfs[d] = log_diff_exp(log_cdfs[d], log_cdfs[d-1]) - log_normalizer;
  }

  return log_pmfs;
}

/**
  * Compute the primary event censored PMF for integer delays up to max_delay
  *
  * @param max_delay Maximum delay to compute PMF for
  * @param D Maximum delay (truncation point), must be at least max_delay + 1
  * @param dist_id Distribution identifier
  * @param params Array of distribution parameters
  * @param pwindow Primary event window
  * @param primary_dist_id Primary distribution identifier
  * @param primary_params Primary distribution parameters
  *
  * @return Vector of primary event censored PMFs for integer delays 1 to
  * max_delay
  *
  * This function differs from primary_censored_dist_pmf in that it:
  * 1. Computes PMFs for all integer delays from \[0, 1\] to \[max_delay,
  *    max_delay + 1\] in one call.
  * 2. Assumes integer delays (swindow = 1)
  * 3. Is more computationally efficient for multiple delay calculations
  *
  * @code
  * // Example: Weibull delay distribution with uniform primary distribution
  * int max_delay = 10;
  * real D = 15.0;
  * int dist_id = 5; // Weibull
  * array[2] real params = {2.0, 1.5}; // shape and scale
  * real pwindow = 7.0;
  * int primary_dist_id = 1; // Uniform
  * array[0] real primary_params = {};
  * vector[max_delay] pmf =
  *   primary_censored_sone_pmf_vectorized(
  *      max_delay, D, dist_id, params, pwindow, primary_dist_id, primary_params
  *   );
  * @endcode
  */
vector primary_censored_sone_pmf_vectorized(
  int max_delay, data real D, int dist_id,
  array[] real params, data real pwindow,
  int primary_dist_id,
  array[] real primary_params
) {
  return exp(
    primary_censored_sone_lpmf_vectorized(
      max_delay, D, dist_id, params, pwindow, primary_dist_id, primary_params
    )
  );
}
