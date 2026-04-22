/**
  * Check if an analytical solution exists for the given distribution combination
  * @ingroup analytical_solution_helpers
  *
  * @param dist_id Distribution identifier for the delay distribution
  * @param primary_id Distribution identifier for the primary distribution
  *
  * @return 1 if an analytical solution exists, 0 otherwise
  */
int check_for_analytical(int dist_id, int primary_id) {
  if (dist_id == 2 && primary_id == 1) return 1; // Gamma delay with Uniform primary
  if (dist_id == 1 && primary_id == 1) return 1; // Lognormal delay with Uniform primary
  if (dist_id == 3 && primary_id == 1) return 1; // Weibull delay with Uniform primary
  return 0; // No analytical solution for other combinations
}

/**
  * Compute the primary event censored log CDF analytically for Gamma delay with Uniform primary
  * @ingroup primary_event_analytical_distributions
  *
  * @param d Delay time
  * @param q Lower bound of integration (max(d - pwindow, 0))
  * @param params Array of Gamma distribution parameters [shape, rate]
  * @param pwindow Primary event window
  *
  * @return Log of the primary event censored CDF for Gamma delay with Uniform
  * primary
  */
real primarycensored_gamma_uniform_lcdf(data real d, real q, array[] real params, data real pwindow) {
  real shape = params[1];
  real rate = params[2];
  real log_window = log(pwindow);
  // log E where E = k * theta = shape / rate is the mean of the delay
  real log_E = log(shape) - log(rate);

  // F_T(d; k) and the recursion to F_T(d; k+1):
  // P(k+1, y) = P(k, y) - y^k e^{-y} / Gamma(k+1), with y = rate * d
  real log_F_T_d_k = gamma_lcdf(d | shape, rate);
  real gamma_kp1_pdf_log_d
    = shape * log(rate * d) - rate * d - lgamma(shape + 1);
  real log_F_T_d_kp1 = log_diff_exp(log_F_T_d_k, gamma_kp1_pdf_log_d);

  // q-dependent terms. Final algebra is unified; only a guard to avoid
  // log_diff_exp(-inf, -inf) and log(0) when q == 0 (q is data, so autodiff
  // is unaffected by this branch).
  real log_q_F_T_q;    // log(q * F_T(q; k))
  real log_E_tF_T_q;   // log(E * F_T(q; k+1))
  if (q > 0) {
    real log_F_T_q_k = gamma_lcdf(q | shape, rate);
    real gamma_kp1_pdf_log_q
      = shape * log(rate * q) - rate * q - lgamma(shape + 1);
    real log_F_T_q_kp1 = log_diff_exp(log_F_T_q_k, gamma_kp1_pdf_log_q);
    log_q_F_T_q = log(q) + log_F_T_q_k;
    log_E_tF_T_q = log_E + log_F_T_q_kp1;
  } else {
    log_q_F_T_q = negative_infinity();
    log_E_tF_T_q = negative_infinity();
  }

  // Unified form: F_{S+}(d) = (A - B) / w_P with A, B sums of positives:
  //   A = d * F_T(d; k)   + E * F_T(q; k+1)
  //   B = q * F_T(q; k)   + E * F_T(d; k+1)
  // Ordering A >= B is guaranteed by F_{S+}(d) >= 0.
  real log_A = log_sum_exp(log(d) + log_F_T_d_k, log_E_tF_T_q);
  real log_B = log_sum_exp(log_q_F_T_q, log_E + log_F_T_d_kp1);

  return log_diff_exp(log_A, log_B) - log_window;
}

/**
  * Compute the primary event censored log CDF analytically for Lognormal delay with Uniform primary
  * @ingroup primary_event_analytical_distributions
  *
  * @param d Delay time
  * @param q Lower bound of integration (max(d - pwindow, 0))
  * @param params Array of Lognormal distribution parameters [mu, sigma]
  * @param pwindow Primary event window
  *
  * @return Log of the primary event censored CDF for Lognormal delay with
  * Uniform primary
  */
real primarycensored_lognormal_uniform_lcdf(data real d, real q, array[] real params, data real pwindow) {
  real mu = params[1];
  real sigma = params[2];
  real mu_sigma2 = mu + square(sigma);
  real log_window = log(pwindow);
  // log E where E = exp(mu + sigma^2/2) is the mean of the delay
  real log_E = mu + 0.5 * square(sigma);

  real log_F_T_d = lognormal_lcdf(d | mu, sigma);
  real log_tF_T_d = lognormal_lcdf(d | mu_sigma2, sigma);

  // q-dependent terms (guard only to avoid log(0); final algebra is unified).
  real log_q_F_T_q;    // log(q * F_T(q))
  real log_E_tF_T_q;   // log(E * tilde F_T(q))
  if (q > 0) {
    real log_F_T_q = lognormal_lcdf(q | mu, sigma);
    real log_tF_T_q = lognormal_lcdf(q | mu_sigma2, sigma);
    log_q_F_T_q = log(q) + log_F_T_q;
    log_E_tF_T_q = log_E + log_tF_T_q;
  } else {
    log_q_F_T_q = negative_infinity();
    log_E_tF_T_q = negative_infinity();
  }

  // Unified form: F_{S+}(d) = (A - B) / w_P with
  //   A = d * F_T(d) + E * tilde F_T(q)
  //   B = q * F_T(q) + E * tilde F_T(d)
  // Ordering A >= B is guaranteed by F_{S+}(d) >= 0.
  real log_A = log_sum_exp(log(d) + log_F_T_d, log_E_tF_T_q);
  real log_B = log_sum_exp(log_q_F_T_q, log_E + log_tF_T_d);

  return log_diff_exp(log_A, log_B) - log_window;
}

/**
  * Compute the log of the lower incomplete gamma function
  * @ingroup analytical_solution_helpers
  *
  * This function is used in the analytical solution for the primary censored
  * Weibull distribution with uniform primary censoring. It corresponds to the
  * g(t; λ, k) function described in the analytic solutions document.
  *
  * @param t Upper bound of integration
  * @param shape Shape parameter (k) of the Weibull distribution
  * @param scale Scale parameter (λ) of the Weibull distribution
  *
  * @return Log of g(t; λ, k) = γ(1 + 1/k, (t/λ)^k)
  */
real log_weibull_g(real t, real shape, real scale) {
  real x = pow(t * inv(scale), shape);
  real a = 1 + inv(shape);
  return log(gamma_p(a, x)) + lgamma(a);
}

/**
  * Compute the primary event censored log CDF analytically for Weibull delay with Uniform primary
  * @ingroup primary_event_analytical_distributions
  *
  * @param d Delay time
  * @param q Lower bound of integration (max(d - pwindow, 0))
  * @param params Array of Weibull distribution parameters [shape, scale]
  * @param pwindow Primary event window
  *
  * @return Log of the primary event censored CDF for Weibull delay with
  * Uniform primary
  */
real primarycensored_weibull_uniform_lcdf(data real d, real q, array[] real params, data real pwindow) {
  real shape = params[1];
  real scale = params[2];
  real log_window = log(pwindow);
  real log_scale = log(scale);

  // For Weibull: E = scale (lambda) and tilde F_T(t) = g(t; lambda, k), so
  // log(E * tilde F_T(t)) = log(scale) + log_weibull_g(t, shape, scale).
  real log_F_T_d = weibull_lcdf(d | shape, scale);
  real log_E_tF_T_d = log_scale + log_weibull_g(d, shape, scale);

  // q-dependent terms (guard only to avoid log(0); final algebra is unified).
  real log_q_F_T_q;    // log(q * F_T(q))
  real log_E_tF_T_q;   // log(E * tilde F_T(q)) = log(scale * g(q; lambda, k))
  if (q > 0) {
    log_q_F_T_q = log(q) + weibull_lcdf(q | shape, scale);
    log_E_tF_T_q = log_scale + log_weibull_g(q, shape, scale);
  } else {
    log_q_F_T_q = negative_infinity();
    log_E_tF_T_q = negative_infinity();
  }

  // Unified form: F_{S+}(d) = (A - B) / w_P with
  //   A = d * F_T(d)    + scale * g(q; lambda, k)
  //   B = q * F_T(q)    + scale * g(d; lambda, k)
  // Ordering A >= B is guaranteed by F_{S+}(d) >= 0.
  real log_A = log_sum_exp(log(d) + log_F_T_d, log_E_tF_T_q);
  real log_B = log_sum_exp(log_q_F_T_q, log_E_tF_T_d);

  return log_diff_exp(log_A, log_B) - log_window;
}

/**
  * Compute the primary event censored log CDF analytically for a single delay
  * (internal version without truncation)
  * @ingroup primary_event_analytical_distributions
  */
real primarycensored_analytical_lcdf_raw(data real d, int dist_id,
                                         array[] real params,
                                         data real pwindow,
                                         int primary_id) {
  real q = max({d - pwindow, 0});

  if (dist_id == 2 && primary_id == 1) {
    return primarycensored_gamma_uniform_lcdf(d | q, params, pwindow);
  } else if (dist_id == 1 && primary_id == 1) {
    return primarycensored_lognormal_uniform_lcdf(d | q, params, pwindow);
  } else if (dist_id == 3 && primary_id == 1) {
    return primarycensored_weibull_uniform_lcdf(d | q, params, pwindow);
  }
  return negative_infinity();
}

/**
  * Compute the primary event censored log CDF analytically for a single delay
  * @ingroup primary_event_analytical_distributions
  *
  * @param d Delay
  * @param dist_id Distribution identifier
  * @param params Array of distribution parameters
  * @param pwindow Primary event window
  * @param L Minimum delay (lower truncation point)
  * @param D Maximum delay (upper truncation point)
  * @param primary_id Primary distribution identifier
  * @param primary_params Primary distribution parameters
  *
  * @return Primary event censored log CDF, normalized over [L, D] if truncation
  * is applied
  */
real primarycensored_analytical_lcdf(data real d, int dist_id,
                                           array[] real params,
                                           data real pwindow, data real L,
                                           data real D, int primary_id,
                                           array[] real primary_params) {
  if (d <= L) return negative_infinity();
  if (d >= D) return 0;

  real result = primarycensored_analytical_lcdf_raw(
    d, dist_id, params, pwindow, primary_id
  );

  // Apply truncation normalization
  if (!is_inf(D) || L > 0) {
    vector[2] bounds = primarycensored_truncation_bounds(
      L, D, dist_id, params, pwindow, primary_id, primary_params
    );
    real log_cdf_L = bounds[1];
    real log_cdf_D = bounds[2];

    real log_normalizer = primarycensored_log_normalizer(log_cdf_D, log_cdf_L, L);
    result = primarycensored_apply_truncation(result, log_cdf_L, log_normalizer, L);
  }

  return result;
}

/**
  * Compute the primary event censored CDF analytically for a single delay
  * @ingroup primary_event_analytical_distributions
  *
  * @param d Delay
  * @param dist_id Distribution identifier
  * @param params Array of distribution parameters
  * @param pwindow Primary event window
  * @param L Minimum delay (lower truncation point)
  * @param D Maximum delay (upper truncation point)
  * @param primary_id Primary distribution identifier
  * @param primary_params Primary distribution parameters
  *
  * @return Primary event censored CDF, normalized over [L, D] if truncation
  * is applied
  */
real primarycensored_analytical_cdf(data real d, int dist_id,
                                          array[] real params,
                                          data real pwindow, data real L,
                                          data real D, int primary_id,
                                          array[] real primary_params) {
  return exp(primarycensored_analytical_lcdf(d | dist_id, params, pwindow, L, D, primary_id, primary_params));
}
