functions {
  /**
   * Check if an analytical solution exists for the given distribution combination
   *
   * @param dist_id Distribution identifier for the delay distribution
   * @param primary_dist_id Distribution identifier for the primary distribution
   *
   * @return 1 if an analytical solution exists, 0 otherwise
   */
  int check_for_analytical(int dist_id, int primary_dist_id) {
    if (dist_id == 2 && primary_dist_id == 1) return 1; // Gamma delay with Uniform primary
    if (dist_id == 1 && primary_dist_id == 1) return 1; // Lognormal delay with Uniform primary
    return 0; // No analytical solution for other combinations
  }

  /**
   * Compute the primary event censored log CDF analytically for Gamma delay with Uniform primary
   */
  real primary_censored_gamma_uniform_lcdf(real d, real q, array[] real params, real pwindow) {
    real shape = params[1];
    real scale = params[2];
    real log_pgamma_q = gamma_lcdf(q | shape, scale);
    real log_pgamma_d = gamma_lcdf(d | shape, scale);
    real log_pgamma_q_1 = gamma_lcdf(q | shape + 1, scale);
    real log_pgamma_d_1 = gamma_lcdf(d | shape + 1, scale);

    real log_Q_T = log1m_exp(log_pgamma_d);
    real log_Delta_F_T_kp1 = log_diff_exp(log_pgamma_d_1, log_pgamma_q_1);
    real log_Delta_F_T_k = log_diff_exp(log_pgamma_d, log_pgamma_q);

    real log_Q_Splus = log_sum_exp(
      log_Q_T,
      log(shape * scale / pwindow) + log_Delta_F_T_kp1,
      log(q / pwindow) + log_Delta_F_T_k
    );

    return log1m_exp(log_Q_Splus);
  }

  /**
   * Compute the primary event censored log CDF analytically for Lognormal delay with Uniform primary
   */
  real primary_censored_lognormal_uniform_lcdf(real d, real q, array[] real params, real pwindow) {
    real mu = params[1];
    real sigma = params[2];
    real log_plnorm_q = lognormal_lcdf(q | mu, sigma);
    real log_plnorm_d = lognormal_lcdf(d | mu, sigma);
    real log_plnorm_q_sigma2 = lognormal_lcdf(q | mu + sigma^2, sigma);
    real log_plnorm_d_sigma2 = lognormal_lcdf(d | mu + sigma^2, sigma);

    real log_Q_T = log1m_exp(log_plnorm_d);
    real log_Delta_F_T_mu_sigma = log_diff_exp(log_plnorm_d_sigma2, log_plnorm_q_sigma2);
    real log_Delta_F_T = log_diff_exp(log_plnorm_d, log_plnorm_q);

    real log_Q_Splus = log_sum_exp(
      log_Q_T,
      log(exp(mu + 0.5 * sigma^2) / pwindow) + log_Delta_F_T_mu_sigma,
      log(q / pwindow) + log_Delta_F_T
    );

    return log1m_exp(log_Q_Splus);
  }

  /**
   * Compute the primary event censored log CDF analytically for a single delay
   *
   * @param d Delay
   * @param dist_id Distribution identifier
   * @param params Array of distribution parameters
   * @param pwindow Primary event window
   * @param D Maximum delay (truncation point)
   * @param primary_dist_id Primary distribution identifier
   * @param primary_params Primary distribution parameters
   * @param lower_bound Lower bound for the delay interval
   *
   * @return Primary event censored log CDF, normalized by D if finite (truncation adjustment)
   */
  real primary_censored_dist_analytical_lcdf(real d, int dist_id, array[] real params,
                                             real pwindow, real D,
                                             int primary_dist_id,
                                             array[] real primary_params,
                                             real lower_bound) {
    real result;
    real q;
    real log_cdf_D;

    if (d <= 0) return negative_infinity();

    q = lower_bound;

    if (dist_id == 2 && primary_dist_id == 1) {
      // Gamma delay with Uniform primary
      result = primary_censored_gamma_uniform_lcdf(d | q, params, pwindow);
    } else if (dist_id == 1 && primary_dist_id == 1) {
      // Lognormal delay with Uniform primary
      result = primary_censored_lognormal_uniform_lcdf(d | q, params, pwindow);
    } else {
      // No analytical solution available
      return negative_infinity();
    }

    if (!is_inf(D)) {
      log_cdf_D = primary_censored_dist_lcdf(D | dist_id, params, pwindow, positive_infinity(), primary_dist_id, primary_params);
      result = result - log_cdf_D;
    }

    return result;
  }

  /**
   * Compute the primary event censored CDF analytically for a single delay
   *
   * @param d Delay
   * @param dist_id Distribution identifier
   * @param params Array of distribution parameters
   * @param pwindow Primary event window
   * @param D Maximum delay (truncation point)
   * @param primary_dist_id Primary distribution identifier
   * @param primary_params Primary distribution parameters
   * @param lower_bound Lower bound for the delay interval
   *
   * @return Primary event censored CDF, normalized by D if finite (truncation adjustment)
   */
  real primary_censored_dist_analytical_cdf(real d, int dist_id, array[] real params,
                                            real pwindow, real D,
                                            int primary_dist_id,
                                            array[] real primary_params,
                                            real lower_bound) {
    return exp(primary_censored_dist_analytical_lcdf(d | dist_id, params, pwindow, D, primary_dist_id, primary_params, lower_bound));
  }
}
