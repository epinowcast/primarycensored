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
    * Compute the primary event censored CDF for a single delay
    *
    * @param d Delay
    * @param params Array of distribution parameters. See `dist_lcdf` for details
    * of what parameters are required for each distribution.
    * @param dist_id Distribution identifier. See `dist_lcdf` for details.
    * @param pwindow Primary event window
    * @param D Maximum delay (truncation point). If finite, the distribution is truncated at D.
    *          If set to positive_infinity(), no truncation is applied.
    *
    * @return Primary event censored CDF, normalized by D if finite (truncation adjustment)
    *
    * @code
    * // Example: Gamma distribution with truncation
    * real d = 3.0;
    * array[2] real params = {2.0, 1.5}; // shape and rate
    * int dist_id = 2; // Gamma
    * real pwindow = 7.0;
    * real D = 10.0; // truncation point
    * real cdf = primary_censored_dist_cdf(d, params, dist_id, pwindow, D);
    */
  real primary_censored_dist_cdf(real d, array[] real params, int dist_id, real pwindow, real D) {
    if (D == positive_infinity()) {
      return integrate_1d(
        function(real p, real xc, array[] real theta, array[] real x_r, array[] int x_i) {
          real d_adj = xc - p;
          return exp(dist_lcdf(d_adj, theta, x_i[1]));
        },
        0, pwindow, d, params, {}, {dist_id}
      );
    } else {
      return integrate_1d(
        function(real p, real xc, array[] real theta, array[] real x_r, array[] int x_i) {
          real d_adj = xc - p;
          real D_adj = x_r[1] - p;
          return exp(
            dist_lcdf(d_adj, theta, x_i[1]) - dist_lcdf(D_adj, theta, x_i[1])
          );
        },
        0, pwindow, d, params, {D}, {dist_id}
      );
    }
  }

  /**
    * Compute the primary event censored log CDF for a single delay
    *
    * @param d Delay
    * @param params Array of distribution parameters. See `dist_lcdf` for details
    * of what parameters are required for each distribution.
    * @param dist_id Distribution identifier. See `dist_lcdf` for details.
    * @param pwindow Primary event window
    * @param D Maximum delay (truncation point). If finite, the distribution is
    * truncated at D.  If set to positive_infinity(), no truncation is applied.
    *
    * @return Primary event censored log CDF, normalized by D if finite (truncation adjustment)
    *
    * @code
    * // Example: Normal distribution without truncation
    * real d = 2.5;
    * array[2] real params = {0.0, 1.0}; // mean and standard deviation
    * int dist_id = 3; // Normal
    * real pwindow = 5.0;
    * real D = positive_infinity(); // no truncation
    * real log_cdf = primary_censored_dist_lcdf(d, params, dist_id, pwindow, D);
    */
  real primary_censored_dist_lcdf(real d, array[] real params, int dist_id, real pwindow, real D) {
    return log(primary_censored_dist_cdf(d, params, dist_id, pwindow, D));
  }

  /**
    * Compute the primary event censored log PMF for a single delay
    *
    * @param d Delay
    * @param params Array of distribution parameters. See `dist_lcdf` for details
    * of what parameters are required for each distribution.
    * @param dist_id Distribution identifier. See `dist_lcdf` for details.
    * @param pwindow Primary event window
    * @param swindow Secondary event window
    * @param D Maximum delay (truncation point). If finite, the distribution is
    * truncated at D. If set to positive_infinity(), no truncation is applied.
    *
    * @return Primary event censored log PMF, normalized by D if finite (truncation adjustment)
    *
    * @code
    * // Example: Exponential distribution with truncation
    * real d = 1.5;
    * array[1] real params = {0.5}; // rate
    * int dist_id = 4; // Exponential
    * real pwindow = 4.0;
    * real swindow = 1.0;
    * real D = 8.0; // truncation point
    * real log_pmf = primary_censored_dist_lpmf(d, params, dist_id, pwindow, swindow, D);
    */
  real primary_censored_dist_lpmf(real d, array[] real params, int dist_id, real pwindow, real swindow, real D) {
    if (d <= 0) {
      reject("Delay must be greater than 0");
    } else {
      return log_diff_exp(
        primary_censored_dist_lcdf(d + swindow, params, dist_id, pwindow, D),
        primary_censored_dist_lcdf(d, params, dist_id, pwindow, D)
      );
    }
  }

  /**
    * Compute the primary event censored PMF for a single delay
    *
    * @param d Delay
    * @param params Array of distribution parameters. See `dist_lcdf` for details
    * of what parameters are required for each distribution.
    * @param dist_id Distribution identifier. See `dist_lcdf` for details.
    * @param pwindow Primary event window
    * @param swindow Secondary event window
    * @param D Maximum delay (truncation point). If finite, the distribution is
    * truncated at D. If set to positive_infinity(), no truncation is applied.
    *
    * @return Primary event censored PMF, normalized by D if finite (truncation adjustment)
    *
    * @code
    * // Example: Weibull distribution without truncation
    * real d = 2.0;
    * array[2] real params = {1.5, 2.0}; // shape and scale
    * int dist_id = 5; // Weibull
    * real pwindow = 6.0;
    * real swindow = 0.5;
    * real D = positive_infinity(); // no truncation
    * real pmf = primary_censored_dist_pmf(d, params, dist_id, pwindow, swindow, D);
    */
  real primary_censored_dist_pmf(real d, array[] real params, int dist_id, real pwindow, real swindow, real D) {
    if (d <= 0) {
      reject("Delay must be greater than 0");
    } else {
      return primary_censored_dist_cdf(d + swindow, params, dist_id, pwindow, D) -
              primary_censored_dist_cdf(d, params, dist_id, pwindow, D);
    }
  }
