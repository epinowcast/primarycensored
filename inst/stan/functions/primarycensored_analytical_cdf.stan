/**
  * Check if the analytical solution is built from uniform primary terms
  * @ingroup analytical_solution_helpers
  *
  * These are the delays whose censored CDF with a uniform primary is
  * primarycensored_uniform_lcdf_from_terms() applied to
  * primarycensored_uniform_terms() at d and q.
  *
  * @param dist_id Distribution identifier for the delay distribution
  * @param primary_id Distribution identifier for the primary distribution
  *
  * @return 1 if the solution is built from uniform primary terms, 0
  * otherwise
  */
int check_for_uniform_terms(int dist_id, int primary_id) {
  if (primary_id != 1) return 0;
  return dist_id == 2 || dist_id == 1 || dist_id == 3 || dist_id == 5;
}

/**
  * Check if an analytical solution exists for the given distribution
  * combination
  * @ingroup analytical_solution_helpers
  *
  * The non-parametric step (26) and discrete-hazard (27, 28) delays are
  * analytic for every primary `primary_lcdf` currently supports, the uniform
  * (1) and exponential growth (2). That list is repeated by hand below, so
  * adding a primary to `primary_lcdf` does not extend the analytic path on
  * its own: without a matching update here the new primary silently falls
  * back to numerical integration.
  *
  * @param dist_id Distribution identifier for the delay distribution
  * @param primary_id Distribution identifier for the primary distribution
  *
  * @return 1 if an analytical solution exists, 0 otherwise
  */
int check_for_analytical(int dist_id, int primary_id) {
  // Gamma, Lognormal, Weibull and generalised gamma with a Uniform primary
  if (check_for_uniform_terms(dist_id, primary_id)) return 1;
  // Keep this primary list in sync with `primary_lcdf`; see the note above.
  if (dist_id == 26 || dist_id == 27 || dist_id == 28) {
    return primary_id == 1 || primary_id == 2;
  }
  return 0; // No analytical solution for other combinations
}

/**
  * Check if the analytical solution applies for the given parameters
  * @ingroup analytical_solution_helpers
  *
  * This is check_for_analytical() and, for the exponentially tilted
  * solutions, which check_for_analytical() does not include, that the tilted
  * delay exists, see check_for_tilt_transform(). It chooses the path in
  * primarycensored_cdf() and primarycensored_lcdf().
  *
  * @param dist_id Distribution identifier for the delay distribution
  * @param params Array of delay distribution parameters
  * @param primary_id Distribution identifier for the primary distribution
  * @param primary_params Array of primary distribution parameters
  *
  * @return 1 if the analytical solution applies, 0 if the numerical path is
  * needed
  */
int check_for_analytical_params(int dist_id, array[] real params,
                                int primary_id,
                                array[] real primary_params) {
  if (check_for_exptilt(dist_id, primary_id)) {
    return check_for_tilt_transform(dist_id, -primary_params[1], params);
  }
  return check_for_analytical(dist_id, primary_id);
}

/**
  * Combine the uniform primary terms at d and q into the censored log CDF
  * @ingroup analytical_solution_helpers
  *
  * For a delay T with mean E and a uniform primary over a window of width
  * w_P, the primary event censored CDF at d is
  *   F_{S+}(d) = (A - B) / w_P, with
  *   A = d * F_T(d) + E * tilde F_T(q),
  *   B = q * F_T(q) + E * tilde F_T(d),
  * where q = max(d - w_P, 0) and tilde F_T is the CDF of the partial
  * expectation distribution. Each of A and B is a sum of one term at d and
  * one at q, and those terms depend on d or q alone (see
  * primarycensored_uniform_terms()). Ordering A >= B is guaranteed by
  * F_{S+}(d) >= 0.
  *
  * @param terms_d Terms at d from primarycensored_uniform_terms()
  * @param terms_q Terms at q from primarycensored_uniform_terms()
  * @param pwindow Primary event window
  *
  * @return Log of the primary event censored CDF at d
  */
real primarycensored_uniform_lcdf_from_terms(vector terms_d, vector terms_q,
                                             data real pwindow) {
  real log_A = log_sum_exp(terms_d[1], terms_q[2]);
  real log_B = log_sum_exp(terms_q[1], terms_d[2]);
  // Deep enough into the lower tail every term underflows together. Both
  // are then `-inf` and `log_diff_exp` would give NaN, so return the limit
  // directly.
  if (log_A == negative_infinity() && log_B == negative_infinity()) {
    return negative_infinity();
  }
  return log_diff_exp(log_A, log_B) - log(pwindow);
}

/**
  * Compute the uniform primary terms at t for a Gamma delay
  * @ingroup analytical_solution_helpers
  *
  * @param t Time (d or q)
  * @param params Array of Gamma distribution parameters [shape, rate]
  *
  * @return Vector [log(t * F_T(t; k)), log(E * F_T(t; k + 1))], both
  * `-inf` for t <= 0
  */
vector primarycensored_gamma_uniform_terms(real t,
                                           array[] real params) {
  if (t <= 0) {
    return rep_vector(negative_infinity(), 2);
  }
  real shape = params[1];
  real rate = params[2];
  // log E where E = k * theta = shape / rate is the mean of the delay
  real log_E = log(shape) - log(rate);
  // F_T(t; k) and the recursion to F_T(t; k+1):
  // P(k+1, y) = P(k, y) - y^k e^{-y} / Gamma(k+1), with y = rate * t
  real log_F_T_k = gamma_lcdf(t | shape, rate);
  real gamma_kp1_pdf_log = shape * log(rate * t) - rate * t
                           - lgamma(shape + 1);
  real log_F_T_kp1 = log_diff_exp(log_F_T_k, gamma_kp1_pdf_log);
  return [log(t) + log_F_T_k, log_E + log_F_T_kp1]';
}

/**
  * Compute the uniform primary terms at t for a Lognormal delay
  * @ingroup analytical_solution_helpers
  *
  * Each term is formed whole and dropped whole. Adding a `-inf` log CDF to
  * the parameter-dependent `log(t)` or `log_E` first would leave an edge
  * back to the parameters that `log_sum_exp` differentiates to
  * `exp(-inf - -inf)`. `t <= 0` underflows on the same test.
  *
  * @param t Time (d or q)
  * @param params Array of Lognormal distribution parameters [mu, sigma]
  *
  * @return Vector [log(t * F_T(t)), log(E * tilde F_T(t))]
  */
vector primarycensored_lognormal_uniform_terms(real t,
                                               array[] real params) {
  real mu = params[1];
  real sigma = params[2];
  real mu_sigma2 = mu + square(sigma);
  // log E where E = exp(mu + sigma^2/2) is the mean of the delay
  real log_E = mu + 0.5 * square(sigma);
  real log_t_F_T = lognormal_lcdf_underflows(t, mu, sigma)
                   ? negative_infinity()
                   : log(t) + lognormal_lcdf(t | mu, sigma);
  real log_E_tF_T = lognormal_lcdf_underflows(t, mu_sigma2, sigma)
                    ? negative_infinity()
                    : log_E + lognormal_lcdf(t | mu_sigma2, sigma);
  return [log_t_F_T, log_E_tF_T]';
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
  // gamma_lcdf(x | a, 1) is log(gamma_p(a, x)), but reverse-mode gamma_p()
  // returns zero gradients for x / a > 10
  // (https://github.com/stan-dev/math/issues/2006). gamma_lcdf() has the
  // same value and computes its own gradients without that cutoff.
  return gamma_lcdf(x | a, 1) + lgamma(a);
}

/**
  * Compute the uniform primary terms at t for a Weibull delay
  * @ingroup analytical_solution_helpers
  *
  * For Weibull, E = scale (lambda) and tilde F_T(t) = g(t; lambda, k).
  *
  * @param t Time (d or q)
  * @param params Array of Weibull distribution parameters [shape, scale]
  *
  * @return Vector [log(t * F_T(t)), log(scale * g(t; lambda, k))], both
  * `-inf` for t <= 0
  */
vector primarycensored_weibull_uniform_terms(real t,
                                             array[] real params) {
  if (t <= 0) {
    return rep_vector(negative_infinity(), 2);
  }
  real shape = params[1];
  real scale = params[2];
  return [
    log(t) + weibull_lcdf(t | shape, scale),
    log(scale) + log_weibull_g(t, shape, scale)
  ]';
}

/**
  * Compute the uniform primary terms at t for a generalised gamma delay
  * @ingroup analytical_solution_helpers
  *
  * Uses the Stacy parameterisation of `flexsurv::pgengamma.orig()`, see
  * `gengamma_lcdf`. The mean is E = scale * Gamma(k + 1/shape) / Gamma(k)
  * and the partial expectation distribution is the generalised gamma with k
  * replaced by k + 1/shape, so this generalises the Gamma (shape = 1) and
  * Weibull (k = 1) solutions.
  *
  * @param t Time (d or q)
  * @param params Array of generalised gamma distribution parameters
  * [shape, scale, k]
  *
  * @return Vector [log(t * F_T(t)), log(E * tilde F_T(t))], both `-inf` for
  * t <= 0
  */
vector primarycensored_gengamma_uniform_terms(real t,
                                              array[] real params) {
  if (t <= 0) {
    return rep_vector(negative_infinity(), 2);
  }
  real shape = params[1];
  real scale = params[2];
  real k = params[3];
  real k_shift = k + inv(shape);
  real log_E = log(scale) + lgamma(k_shift) - lgamma(k);
  return [
    log(t) + gengamma_lcdf(t | shape, scale, k),
    log_E + gengamma_lcdf(t | shape, scale, k_shift)
  ]';
}

/**
  * Compute the uniform primary terms at t for a delay distribution
  * @ingroup analytical_solution_helpers
  *
  * @param t Time (d or q)
  * @param dist_id Distribution identifier (1: Lognormal, 2: Gamma,
  *   3: Weibull, 5: Generalised gamma), see check_for_uniform_terms()
  * @param params Array of distribution parameters
  *
  * @return Vector of the two terms at t, see
  * primarycensored_uniform_lcdf_from_terms()
  */
vector primarycensored_uniform_terms(real t, data int dist_id,
                                     array[] real params) {
  if (dist_id == 2) {
    return primarycensored_gamma_uniform_terms(t, params);
  } else if (dist_id == 1) {
    return primarycensored_lognormal_uniform_terms(t, params);
  } else if (dist_id == 3) {
    return primarycensored_weibull_uniform_terms(t, params);
  } else if (dist_id == 5) {
    return primarycensored_gengamma_uniform_terms(t, params);
  }
  reject("Invalid distribution identifier: ", dist_id);
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
real primarycensored_gamma_uniform_lcdf(data real d, real q,
                                        array[] real params,
                                        data real pwindow) {
  return primarycensored_uniform_lcdf_from_terms(
    primarycensored_gamma_uniform_terms(d, params),
    primarycensored_gamma_uniform_terms(q, params), pwindow
  );
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
real primarycensored_lognormal_uniform_lcdf(data real d, real q,
                                            array[] real params,
                                            data real pwindow) {
  return primarycensored_uniform_lcdf_from_terms(
    primarycensored_lognormal_uniform_terms(d, params),
    primarycensored_lognormal_uniform_terms(q, params), pwindow
  );
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
real primarycensored_weibull_uniform_lcdf(data real d, real q,
                                          array[] real params,
                                          data real pwindow) {
  return primarycensored_uniform_lcdf_from_terms(
    primarycensored_weibull_uniform_terms(d, params),
    primarycensored_weibull_uniform_terms(q, params), pwindow
  );
}

/**
  * Compute the primary event censored log CDF analytically for generalised gamma delay with Uniform primary
  * @ingroup primary_event_analytical_distributions
  *
  * @param d Delay time
  * @param q Lower bound of integration (max(d - pwindow, 0))
  * @param params Array of generalised gamma distribution parameters
  * [shape, scale, k]
  * @param pwindow Primary event window
  *
  * @return Log of the primary event censored CDF for generalised gamma delay
  * with Uniform primary
  */
real primarycensored_gengamma_uniform_lcdf(data real d, real q,
                                           array[] real params,
                                           data real pwindow) {
  return primarycensored_uniform_lcdf_from_terms(
    primarycensored_gengamma_uniform_terms(d, params),
    primarycensored_gengamma_uniform_terms(q, params), pwindow
  );
}

/**
  * Compute the primary event censored log CDF analytically for a single delay
  * (internal version without truncation)
  * @ingroup primary_event_analytical_distributions
  */
real primarycensored_analytical_lcdf_raw(data real d, int dist_id,
                                         array[] real params,
                                         data real pwindow,
                                         int primary_id,
                                         array[] real primary_params) {
  real q = max({d - pwindow, 0});

  if (check_for_exptilt(dist_id, primary_id)) {
    if (!check_for_tilt_transform(dist_id, -primary_params[1], params)) {
      reject(
        "The tilted delay distribution does not exist for tilt ",
        primary_params[1], ". Use the numerical path, see ",
        "check_for_analytical_params()."
      );
    }
    return primarycensored_exptilt_lcdf(
      d | dist_id, params, pwindow, primary_params[1]
    );
  }
  if (dist_id == 2 && primary_id == 1) {
    return primarycensored_gamma_uniform_lcdf(d | q, params, pwindow);
  } else if (dist_id == 1 && primary_id == 1) {
    return primarycensored_lognormal_uniform_lcdf(d | q, params, pwindow);
  } else if (dist_id == 3 && primary_id == 1) {
    return primarycensored_weibull_uniform_lcdf(d | q, params, pwindow);
  } else if (dist_id == 5 && primary_id == 1) {
    return primarycensored_gengamma_uniform_lcdf(d | q, params, pwindow);
  } else if (dist_id == 26) {
    // params = [boundaries (K+1), pmf (K)]; length 2*K + 1.
    int K = (size(params) - 1) %/% 2;
    return discretestep_lcdf(
      d | to_vector(segment(params, 1, K + 1)),
          to_vector(segment(params, K + 2, K)),
          primary_id, primary_params, pwindow
    );
  } else if (dist_id == 27 || dist_id == 28) {
    // params = [boundaries (K+1), hazards (K)]; length 2*K + 1. The last
    // hazard must equal 1. RW (27) and RE (28) only differ in their
    // prior so they share this likelihood dispatch.
    int K = (size(params) - 1) %/% 2;
    return discretehazard_lcdf(
      d | to_vector(segment(params, 1, K + 1)),
          to_vector(segment(params, K + 2, K)),
          primary_id, primary_params, pwindow
    );
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
    d, dist_id, params, pwindow, primary_id, primary_params
  );

  // Apply truncation normalization, also for a finite negative L for delays
  // on the reals
  if (!is_inf(D) || L > 0
      || (!is_inf(L) && !dist_has_positive_support(dist_id))) {
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

/**
  * Check if the analytical solution can be vectorised over integer delays
  * @ingroup analytical_solution_helpers
  *
  * The analytical CDF at d combines terms at d and at q = max(d - pwindow, 0).
  * With an integer pwindow q is an integer delay too, so
  * primarycensored_analytical_lcdf_vectorized() can compute the terms once
  * per delay and share them. This needs the analytical solutions built from
  * primarycensored_uniform_terms(), see check_for_uniform_terms(), or the
  * exponentially tilted solutions, see check_for_exptilt(). The
  * non-parametric delays in check_for_analytical() have no such terms.
  *
  * @param dist_id Distribution identifier for the delay distribution
  * @param primary_id Distribution identifier for the primary distribution
  * @param pwindow Primary event window
  *
  * @return 1 if the vectorised analytical solution applies, 0 otherwise
  */
int check_for_analytical_vectorized(int dist_id, int primary_id,
                                    data real pwindow) {
  return (check_for_uniform_terms(dist_id, primary_id)
          || check_for_exptilt(dist_id, primary_id)) &&
    pwindow >= 1 && floor(pwindow) == pwindow;
}

/**
  * Compute the primary event censored log CDF analytically at integer delays
  * @ingroup primary_event_analytical_distributions
  *
  * The log CDF at d combines the terms at d and at q = max(d - pwindow, 0)
  * (see primarycensored_uniform_lcdf_from_terms()). Both are integer delays,
  * so the terms are computed once per delay and used for both, halving the
  * CDF evaluations. The values are the same as from
  * primarycensored_analytical_lcdf() at each delay without truncation.
  * For an exponentially tilted primary the terms are those of
  * primarycensored_exptilt_lcdf(), with the form chosen at each delay as
  * there.
  * Only for cases where check_for_analytical_vectorized() and
  * check_for_analytical_params() are 1.
  *
  * @param start First delay to compute
  * @param n Last delay to compute, and the length of the result
  * @param dist_id Distribution identifier
  * @param params Array of distribution parameters
  * @param pwindow Primary event window, a positive integer
  * @param primary_id Primary distribution identifier
  * @param primary_params Primary distribution parameters
  *
  * @return Vector whose element d is the log CDF at d, for d in start:n.
  * Elements before start are not computed.
  */
vector primarycensored_analytical_lcdf_vectorized(
  data int start, data int n, data int dist_id, array[] real params,
  data real pwindow, data int primary_id, array[] real primary_params
) {
  int pw = to_int(pwindow);
  vector[n] log_cdfs;
  if (check_for_exptilt(dist_id, primary_id)) {
    real rho = primary_params[1];
    int positive = dist_has_positive_support(dist_id);
    // Endpoints below 0 share the entry for 0 for non-negative delays
    int first = positive ? max(start - pw, 0) : start - pw;
    if (abs(rho) * pwindow < (positive ? 1e-2 : 1e-3)) {
      // moments[t - first + 1] holds the moments at endpoint t
      array[n - first + 1] vector[3] moments;
      for (t in first:n) {
        moments[t - first + 1] = primarycensored_tilt_moments(
          t, dist_id, params
        );
      }
      for (d in start:n) {
        int q_index = (positive ? max(d - pw, 0) : d - pw) - first + 1;
        log_cdfs[d] = primarycensored_exptilt_small_window_lcdf_from_terms(
          moments[d - first + 1], moments[q_index], rho, pwindow
        );
      }
    } else {
      // terms[t - first + 1] holds the terms at endpoint t
      array[n - first + 1] vector[4] terms;
      for (t in first:n) {
        terms[t - first + 1] = append_row(
          log_tilt_transform_pair(t, dist_id, 0, params),
          log_tilt_transform_pair(t, dist_id, -rho, params)
        );
      }
      for (d in start:n) {
        if (positive && d < pwindow && abs(rho) * d < 1e-2) {
          log_cdfs[d] = primarycensored_exptilt_small_delay_lcdf_from_terms(
            primarycensored_tilt_moments(d, dist_id, params), rho, pwindow
          );
        } else {
          int q_index = (positive ? max(d - pw, 0) : d - pw) - first + 1;
          log_cdfs[d] = primarycensored_exptilt_lcdf_from_terms(
            terms[d - first + 1], terms[q_index], d, rho, pwindow
          );
        }
      }
    }
    return log_cdfs;
  }
  // terms[t + 1] holds the terms at delay t
  array[n + 1] vector[2] terms;
  for (t in max(start - pw, 0):n) {
    terms[t + 1] = primarycensored_uniform_terms(t, dist_id, params);
  }
  for (d in start:n) {
    log_cdfs[d] = primarycensored_uniform_lcdf_from_terms(
      terms[d + 1], terms[max(d - pw, 0) + 1], pwindow
    );
  }
  return log_cdfs;
}

/*
 * Truncated exponential-moment transforms of delay distributions
 *
 * T_f(xi; tau) = int_{-inf}^{tau} exp(xi u) f(u) du, with lower limit 0 for
 * delays on the non-negative reals. A delay is added by a branch in
 * check_for_tilt_transform(), log_tilt_transform_pair() and
 * primarycensored_tilt_moments(). The R equivalents are the
 * `.pcens_tilt_*()` generics.
 */

/**
  * Log of the difference of two exponentials, zero when it would be negative
  * @ingroup tilt_transforms
  *
  * Unlike log_diff_exp() this is finite for a difference that is zero to
  * rounding, and a zero subtrahend returns the first argument.
  *
  * @param a Log of the larger term
  * @param b Log of the smaller term
  *
  * @return log(exp(a) - exp(b)), or `-inf` if a <= b
  */
real primarycensored_log_diff_exp(real a, real b) {
  if (b == negative_infinity()) return a;
  if (a <= b) return negative_infinity();
  return log_diff_exp(a, b);
}

/**
  * Log of the standard normal CDF with an exact derivative
  * @ingroup tilt_transforms
  *
  * The derivative of std_normal_lcdf() has a relative error of about 1e-5,
  * which small tilt differences amplify. For z < 0 this uses erfc(), and
  * below -37 the asymptotic series
  * Phi(z) = phi(z) / (-z) (1 - 1 / z^2 + 3 / z^4 - ...), so that the value
  * and the derivative are exact.
  *
  * @param z Point
  *
  * @return log(Phi(z))
  */
real primarycensored_log_std_normal_cdf(real z) {
  if (z >= 0) return log(Phi(z));
  if (z > -37) return log(0.5 * erfc(-z / sqrt(2)));
  real inv_z2 = 1 / square(z);
  real series = 1 - inv_z2 * (1 - 3 * inv_z2 * (1 - 5 * inv_z2 * (
    1 - 7 * inv_z2 * (1 - 9 * inv_z2 * (1 - 11 * inv_z2 * (
      1 - 13 * inv_z2)))
  )));
  return std_normal_lpdf(z) - log(-z) + log(series);
}

/**
  * Log of the regularised lower incomplete gamma function from its series
  * @ingroup tilt_transforms
  *
  * Built from elementary operations so that the shape derivative is exact,
  * unlike gamma_lcdf(), whose shape derivative is inaccurate in parts of the
  * bulk and NaN for shapes of about 200 or more. The series
  * P(shape, x) = x^shape exp(-x) / Gamma(shape + 1) *
  *   sum_k x^k / ((shape + 1) ... (shape + k))
  * is for x < shape + 1.
  *
  * @param x Point, positive, below shape + 1
  * @param shape Shape, positive
  *
  * @return log P(shape, x)
  */
real primarycensored_log_gamma_p_series(real x, real shape) {
  real term = 1;
  real total = 1;
  real max_terms = 10 * sqrt(shape) + 150;
  int k = 0;
  while (term >= 1e-17 * total) {
    k += 1;
    if (k > max_terms) {
      reject("The gamma lower series did not converge for x = ", x,
             " and shape = ", shape);
    }
    term *= x / (shape + k);
    total += term;
  }
  return shape * log(x) - x - lgamma(shape + 1) + log(total);
}

/**
  * Log of the regularised upper incomplete gamma function from its
  * continued fraction
  * @ingroup tilt_transforms
  *
  * The Legendre continued fraction, evaluated by the modified Lentz method,
  * for x >= shape + 1. As for primarycensored_log_gamma_p_series() the
  * derivatives are exact, unlike gamma_lccdf().
  *
  * @param x Point, positive, at least shape + 1
  * @param shape Shape, positive
  *
  * @return log Q(shape, x)
  */
real primarycensored_log_gamma_q_fraction(real x, real shape) {
  real b = x + 1 - shape;
  real d = 1 / b;
  real c = 1e300;
  real h = d;
  real max_terms = 10 * sqrt(shape) + 150;
  real del = 0;
  int i = 0;
  while (abs(del - 1) > 1e-15) {
    i += 1;
    if (i > max_terms) {
      reject("The gamma upper fraction did not converge for x = ", x,
             " and shape = ", shape);
    }
    real an = -i * (i - shape);
    b += 2;
    d = an * d + b;
    if (abs(d) < 1e-300) d = 1e-300;
    c = b + an / c;
    if (abs(c) < 1e-300) c = 1e-300;
    d = 1 / d;
    del = d * c;
    h *= del;
  }
  return shape * log(x) - x - lgamma(shape) + log(h);
}

/**
  * Log of the regularised lower and upper incomplete gamma functions
  * @ingroup tilt_transforms
  *
  * The tail that is not close to 1 is evaluated directly and the other is
  * log1m_exp() of it, including for an upper tail below the smallest double.
  *
  * @param x Point, positive
  * @param shape Shape, positive
  *
  * @return Vector [log P(shape, x), log Q(shape, x)]
  */
vector primarycensored_log_gamma_pq(real x, real shape) {
  if (x < shape + 1) {
    real log_lower = primarycensored_log_gamma_p_series(x, shape);
    return [log_lower, log1m_exp(log_lower)]';
  }
  real log_upper = primarycensored_log_gamma_q_fraction(x, shape);
  return [log1m_exp(log_upper), log_upper]';
}

/**
  * Test whether the gamma lower tail underflows at these arguments
  * @ingroup tilt_transforms
  *
  * Terms below 1e-300 of a probability are dropped as `-inf`. For
  * x < shape + 1 this uses the bound
  * P(shape, x) <= x^shape exp(-x) / (Gamma(shape + 1) (1 - x / (shape + 1))).
  *
  * @param x Point, positive
  * @param shape Shape
  *
  * @return 1 if the regularised lower incomplete gamma function may underflow,
  *   0 otherwise
  */
int gamma_lcdf_underflows(real x, real shape) {
  if (x >= shape + 1) return 0;
  return shape * log(x) - lgamma(shape + 1) - x - log1m(x / (shape + 1))
         < -700;
}

/**
  * Test whether the gamma upper tail underflows at these arguments
  * @ingroup tilt_transforms
  *
  * As for gamma_lcdf_underflows(), using
  * Q(shape, x) = x^(shape - 1) exp(-x) / Gamma(shape) x / (x - shape + 1)
  * for x beyond shape + 1.
  *
  * @param x Point, positive
  * @param shape Shape
  *
  * @return 1 if the regularised upper incomplete gamma function may underflow,
  *   0 otherwise
  */
int gamma_lccdf_underflows(real x, real shape) {
  if (x <= shape + 1) return 0;
  return (shape - 1) * log(x) - x - lgamma(shape) + log(x / (x - shape + 1))
         < -700;
}

/**
  * Check if the tilt transform is closed form for a delay and tilt
  * @ingroup tilt_transforms
  *
  * The exponential (4) and gamma (2) forms need the delay with the rate
  * lowered by xi, which exists if rate - xi > 0. The normal (18) form has no
  * restriction. Callers use the numerical path when this is 0.
  *
  * @param dist_id Distribution identifier for the delay distribution
  * @param xi Tilt. The exponentially tilted window with tilt rho needs
  *   xi = -rho
  * @param params Array of distribution parameters, as for dist_lcdf()
  *
  * @return 1 if the transform is closed form and the tilted delay exists, 0
  * otherwise
  */
int check_for_tilt_transform(int dist_id, real xi, array[] real params) {
  if (dist_id == 4) return params[1] - xi > 0;
  if (dist_id == 2) return params[2] - xi > 0;
  if (dist_id == 18) return 1;
  return 0;
}

/**
  * Log of the tilt transform over the lower and the upper part of the support
  * @ingroup tilt_transforms
  *
  * The lower transform is T_f(xi; t), the delay CDF for xi = 0, and `-inf`
  * for t <= 0 for delays on the non-negative reals. The upper transform is
  * T_f(xi; Inf) - T_f(xi; t). Only defined where check_for_tilt_transform()
  * is 1.
  *
  * The gamma is the gamma CDF with the rate lowered by xi times the total
  * (rate / (rate - xi))^shape, with both tails from
  * primarycensored_log_gamma_pq(). The log total is -shape log1m(xi / rate)
  * rather than a difference of logs, so the rate derivative does not lose a
  * tail far below 1 to cancellation.
  *
  * @param t Point
  * @param dist_id Distribution identifier: 2 (Gamma), 4 (Exponential) or 18
  *   (Normal), see check_for_tilt_transform()
  * @param xi Tilt
  * @param params Array of distribution parameters, as for dist_lcdf()
  *
  * @return Vector [log T_f(xi; t), log(T_f(xi; Inf) - T_f(xi; t))]
  */
vector log_tilt_transform_pair(real t, int dist_id, real xi,
                               array[] real params) {
  if (dist_id == 2) {
    real shape = params[1];
    real rate = params[2];
    real tilted_rate = rate - xi;
    real log_total = -shape * log1m(xi / rate);
    if (t <= 0) return [negative_infinity(), log_total]';
    real x = t * tilted_rate;
    if (gamma_lcdf_underflows(x, shape)) {
      return [negative_infinity(), log_total]';
    }
    if (gamma_lccdf_underflows(x, shape)) {
      return [log_total, negative_infinity()]';
    }
    vector[2] log_tails = primarycensored_log_gamma_pq(x, shape);
    return [log_total + log_tails[1], log_total + log_tails[2]]';
  } else if (dist_id == 4) {
    real rate = params[1];
    real tilted_rate = rate - xi;
    real log_total = -log1m(xi / rate);
    if (t <= 0) return [negative_infinity(), log_total]';
    return [
      log_total + log1m_exp(-tilted_rate * t), log_total - tilted_rate * t
    ]';
  } else if (dist_id == 18) {
    // The upper tail is the lower tail of the reflected normal
    real mu = params[1];
    real sigma = params[2];
    real z = (t - mu - xi * square(sigma)) / sigma;
    real log_total = xi * mu + 0.5 * square(xi * sigma);
    return [
      log_total + primarycensored_log_std_normal_cdf(z),
      log_total + primarycensored_log_std_normal_cdf(-z)
    ]';
  }
  reject("Invalid distribution identifier: ", dist_id);
}

/**
  * Log moments of a gamma delay about a point
  * @ingroup tilt_transforms
  *
  * The log of G_k(t) = int_0^t (t - u)^k f(u) du for k = 1, 2, 3, from gamma
  * CDFs with the shape raised by k.
  *
  * @param t Point, positive
  * @param shape Shape
  * @param rate Rate
  *
  * @return Vector [log G_1(t), log G_2(t), log G_3(t)]
  */
vector primarycensored_gamma_tilt_moments(real t, real shape, real rate) {
  if (gamma_lcdf_underflows(t * rate, shape)
      || gamma_lcdf_underflows(t * rate, shape + 1)
      || gamma_lcdf_underflows(t * rate, shape + 2)
      || gamma_lcdf_underflows(t * rate, shape + 3)) {
    return rep_vector(negative_infinity(), 3);
  }
  // The CDFs are 1, so these are moments of the whole distribution
  if (gamma_lccdf_underflows(t * rate, shape + 3)) {
    real m1 = shape / rate;
    real m2 = shape * (shape + 1) / square(rate);
    real m3 = shape * (shape + 1) * (shape + 2) / pow(rate, 3);
    return [
      log(t - m1),
      log(square(t) - 2 * t * m1 + m2),
      log(pow(t, 3) - 3 * square(t) * m1 + 3 * t * m2 - m3)
    ]';
  }
  real log_t = log(t);
  real x = t * rate;
  real log_m0 = primarycensored_log_gamma_pq(x, shape)[1];
  real log_m1 = log(shape) - log(rate)
                + primarycensored_log_gamma_pq(x, shape + 1)[1];
  real log_m2 = log(shape) + log(shape + 1) - 2 * log(rate)
                + primarycensored_log_gamma_pq(x, shape + 2)[1];
  real log_m3 = log(shape) + log(shape + 1) + log(shape + 2) - 3 * log(rate)
                + primarycensored_log_gamma_pq(x, shape + 3)[1];
  real log_g1 = primarycensored_log_diff_exp(log_t + log_m0, log_m1);
  real log_h = primarycensored_log_diff_exp(log_t + log_m1, log_m2);
  real log_g2 = primarycensored_log_diff_exp(log_t + log_g1, log_h);
  real log_a = primarycensored_log_diff_exp(log_t + log_m2, log_m3);
  real log_b = primarycensored_log_diff_exp(log_t + log_h, log_a);
  real log_g3 = primarycensored_log_diff_exp(log_t + log_g2, log_b);
  return [log_g1, log_g2, log_g3]';
}

/**
  * Log moments of a delay about a point
  * @ingroup tilt_transforms
  *
  * The log of G_k(t) = int (t - u)^k f(u) du up to t for k = 1, 2, 3, used by
  * the small tilt forms. `-inf` for t <= 0 for delays on the non-negative
  * reals.
  *
  * The exponential is the gamma with shape 1. The normal with
  * z = (t - mu) / sigma has G_1 = sigma (phi(z) + z Phi(z)),
  * G_2 = sigma^2 ((z^2 + 1) Phi(z) + z phi(z)) and
  * G_3 = sigma^3 (z (z^2 + 3) Phi(z) + (z^2 + 2) phi(z)).
  *
  * @param t Point
  * @param dist_id Distribution identifier: 2 (Gamma), 4 (Exponential) or 18
  *   (Normal), see check_for_tilt_transform()
  * @param params Array of distribution parameters, as for dist_lcdf()
  *
  * @return Vector [log G_1(t), log G_2(t), log G_3(t)]
  */
vector primarycensored_tilt_moments(real t, int dist_id,
                                    array[] real params) {
  if (dist_id == 2) {
    if (t <= 0) return rep_vector(negative_infinity(), 3);
    return primarycensored_gamma_tilt_moments(t, params[1], params[2]);
  } else if (dist_id == 4) {
    if (t <= 0) return rep_vector(negative_infinity(), 3);
    return primarycensored_gamma_tilt_moments(t, 1, params[1]);
  } else if (dist_id == 18) {
    real mu = params[1];
    real sigma = params[2];
    real z = (t - mu) / sigma;
    real log_phi = std_normal_lpdf(z);
    real log_Phi = primarycensored_log_std_normal_cdf(z);
    real log_g1;
    real log_g2;
    real log_g3;
    if (z > -1) {
      // The direct form has a derivative at z = 0, unlike the log form
      log_g1 = log(exp(log_phi) + z * exp(log_Phi));
      log_g2 = log((square(z) + 1) * exp(log_Phi) + z * exp(log_phi));
      log_g3 = log(
        (square(z) + 2) * exp(log_phi) + z * (square(z) + 3) * exp(log_Phi)
      );
    } else {
      // The terms have opposite signs, so the log form keeps the tail
      log_g1 = primarycensored_log_diff_exp(log_phi, log(-z) + log_Phi);
      log_g2 = primarycensored_log_diff_exp(
        log1p(square(z)) + log_Phi, log(-z) + log_phi
      );
      log_g3 = primarycensored_log_diff_exp(
        log(square(z) + 2) + log_phi, log(-z) + log(square(z) + 3) + log_Phi
      );
    }
    return [
      log(sigma) + log_g1, 2 * log(sigma) + log_g2, 3 * log(sigma) + log_g3
    ]';
  }
  reject("Invalid distribution identifier: ", dist_id);
}

/*
 * Exponentially tilted primary event window
 *
 * For a window of width w with density rho exp(rho z) / (exp(rho w) - 1) on
 * [0, w], with q = d - w, F the delay CDF and J(x) = T_f(-rho; x) the tilt
 * transform above,
 *   F_rho(d) = F(q) + (exp(rho d) (J(d) - J(q)) - (F(d) - F(q))) /
 *              (exp(rho w) - 1).
 * Every term depends on d or q alone, so the terms are computed once per
 * endpoint in primarycensored_analytical_lcdf_vectorized().
 *
 * The direct form cancels as rho goes to zero, see ?pcens_cdf_exptilt. For
 * |rho| * pwindow below 1e-2 (1e-3 for delays on the reals) the small window
 * form replaces it. For delays on the non-negative reals with d < pwindow and
 * |rho| * d below 1e-2 the small delay form replaces it.
 */

/**
  * Check if the analytical solution is the exponentially tilted solution
  * @ingroup exponential_tilt_solutions
  *
  * The tilted primary is primary_id 2 with primary_params = [rho]. Whether
  * the solution applies for given parameters is check_for_analytical_params().
  *
  * @param dist_id Distribution identifier for the delay distribution
  * @param primary_id Distribution identifier for the primary distribution
  *
  * @return 1 if the delay has an exponentially tilted solution and the
  * primary is exponentially tilted, 0 otherwise
  */
int check_for_exptilt(int dist_id, int primary_id) {
  return primary_id == 2 && (dist_id == 2 || dist_id == 4 || dist_id == 18);
}

/**
  * Log of a difference between two points from either tail
  * @ingroup exponential_tilt_solutions
  *
  * With L + U constant, L(d) - L(q) equals U(q) - U(d). This uses the form
  * with the smaller ratio of the terms.
  *
  * @param lower_d Log lower tail quantity at d
  * @param lower_q Log lower tail quantity at q
  * @param upper_d Log upper tail quantity at d
  * @param upper_q Log upper tail quantity at q
  *
  * @return Log of the difference, `-inf` if it is zero to rounding
  */
real primarycensored_tail_diff(real lower_d, real lower_q, real upper_d,
                               real upper_q) {
  // NaN from terms that both underflow is a zero difference
  if (is_nan(lower_q - lower_d) || lower_q - lower_d <= upper_d - upper_q) {
    return primarycensored_log_diff_exp(lower_d, lower_q);
  }
  return primarycensored_log_diff_exp(upper_q, upper_d);
}

/**
  * Combine the exponentially tilted terms at d and q into the log CDF
  * @ingroup exponential_tilt_solutions
  *
  * The direct form, for where the small window and small delay forms do not
  * apply. The difference of each pair of terms is taken with
  * primarycensored_tail_diff().
  *
  * @param terms_d Terms at d, the lower and upper terms of
  *   log_tilt_transform_pair() for xi = 0 and xi = -rho
  * @param terms_q Terms at q = d - pwindow
  * @param d Delay
  * @param rho Tilt, not zero
  * @param pwindow Primary event window
  *
  * @return Log of the primary event censored CDF at d
  */
real primarycensored_exptilt_lcdf_from_terms(vector terms_d, vector terms_q,
                                             data real d, real rho,
                                             data real pwindow) {
  real log_diff_f = primarycensored_tail_diff(
    terms_d[1], terms_q[1], terms_d[2], terms_q[2]
  );
  real log_diff_j = primarycensored_tail_diff(
    terms_d[3], terms_q[3], terms_d[4], terms_q[4]
  );
  real log_num;
  real log_den;
  if (rho > 0) {
    log_num = primarycensored_log_diff_exp(rho * d + log_diff_j, log_diff_f);
    log_den = log_diff_exp(rho * pwindow, 0);
  } else {
    log_num = primarycensored_log_diff_exp(log_diff_f, rho * d + log_diff_j);
    log_den = log1m_exp(rho * pwindow);
  }
  // log_sum_exp of two `-inf` has a NaN derivative
  if (terms_q[1] == negative_infinity() && log_num == negative_infinity()) {
    return negative_infinity();
  }
  // Rounding can put the log CDF above 0
  return fmin(log_sum_exp(terms_q[1], log_num - log_den), 0);
}

/**
  * Combine the moments at d and q into the small tilt log CDF
  * @ingroup exponential_tilt_solutions
  *
  * The uniform window limit with its corrections to second order in the tilt,
  * for the small window form. With G_k(t) = int (t - u)^k f(u) du,
  *   F_rho(d) = (G_1(d) - G_1(q)) / w
  *     + rho (G_2(d) - w G_1(d) - G_2(q) - w G_1(q)) / (2 w)
  *     + rho^2 ((G_3(d) - G_3(q)) / 6 - w (G_2(d) + G_2(q)) / 4
  *              + w^2 (G_1(d) - G_1(q)) / 12) / w.
  *
  * @param moments_d Moments [log G_1, log G_2, log G_3] at d from
  *   primarycensored_tilt_moments()
  * @param moments_q Moments at q = d - pwindow
  * @param rho Tilt
  * @param pwindow Primary event window
  *
  * @return Log of the primary event censored CDF at d
  */
real primarycensored_exptilt_small_window_lcdf_from_terms(
  vector moments_d, vector moments_q, real rho, data real pwindow
) {
  real scale = moments_d[1];
  if (scale == negative_infinity()) return negative_infinity();
  real g1_q = exp(moments_q[1] - scale);
  real g2_d = exp(moments_d[2] - scale);
  real g2_q = exp(moments_q[2] - scale);
  real g3_d = exp(moments_d[3] - scale);
  real g3_q = exp(moments_q[3] - scale);
  real relative = (1 - g1_q)
                  + 0.5 * rho * (g2_d - pwindow - g2_q - pwindow * g1_q)
                  + square(rho) * ((g3_d - g3_q) / 6
                                   - pwindow * (g2_d + g2_q) / 4
                                   + square(pwindow) * (1 - g1_q) / 12);
  if (relative <= 0) return negative_infinity();
  return fmin(scale + log(relative) - log(pwindow), 0);
}

/**
  * Compute the small delay log CDF from the moments at d
  * @ingroup exponential_tilt_solutions
  *
  * For the small delay form the terms at d - pwindow are zero and
  * F_rho(d) = rho (G_1(d) + rho G_2(d) / 2 + rho^2 G_3(d) / 6) /
  *   (exp(rho w) - 1).
  *
  * @param moments_d Moments [log G_1, log G_2, log G_3] at d from
  *   primarycensored_tilt_moments()
  * @param rho Tilt, not zero
  * @param pwindow Primary event window
  *
  * @return Log of the primary event censored CDF at d
  */
real primarycensored_exptilt_small_delay_lcdf_from_terms(
  vector moments_d, real rho, data real pwindow
) {
  if (moments_d[1] == negative_infinity()) return negative_infinity();
  real log_den = rho > 0 ? log_diff_exp(rho * pwindow, 0)
                         : log1m_exp(rho * pwindow);
  return fmin(
    moments_d[1] + log(abs(rho)) - log_den
    + log1p(0.5 * rho * exp(moments_d[2] - moments_d[1])
            + square(rho) * exp(moments_d[3] - moments_d[1]) / 6),
    0
  );
}

/**
  * Compute the primary event censored log CDF for an exponentially tilted
  * primary
  * @ingroup exponential_tilt_solutions
  *
  * Chooses the direct form or the small tilt forms. A zero width window
  * gives the delay log CDF. Only for check_for_exptilt() is 1 and
  * check_for_tilt_transform() is 1 for -rho.
  *
  * @param d Delay
  * @param dist_id Distribution identifier: 2 (Gamma), 4 (Exponential) or 18
  *   (Normal)
  * @param params Array of distribution parameters
  * @param pwindow Primary event window
  * @param rho Tilt, the exponential growth rate of the primary
  *
  * @return Log of the primary event censored CDF at d
  */
real primarycensored_exptilt_lcdf(data real d, int dist_id,
                                  array[] real params, data real pwindow,
                                  real rho) {
  int positive = dist_has_positive_support(dist_id);
  if (positive && d <= 0) {
    return negative_infinity();
  }
  if (pwindow == 0) return dist_lcdf(d | params, dist_id);
  real q = d - pwindow;
  if (abs(rho) * pwindow < (positive ? 1e-2 : 1e-3)) {
    return primarycensored_exptilt_small_window_lcdf_from_terms(
      primarycensored_tilt_moments(d, dist_id, params),
      primarycensored_tilt_moments(q, dist_id, params), rho, pwindow
    );
  }
  if (positive && d < pwindow && abs(rho) * d < 1e-2) {
    return primarycensored_exptilt_small_delay_lcdf_from_terms(
      primarycensored_tilt_moments(d, dist_id, params), rho, pwindow
    );
  }
  return primarycensored_exptilt_lcdf_from_terms(
    append_row(
      log_tilt_transform_pair(d, dist_id, 0, params),
      log_tilt_transform_pair(d, dist_id, -rho, params)
    ),
    append_row(
      log_tilt_transform_pair(q, dist_id, 0, params),
      log_tilt_transform_pair(q, dist_id, -rho, params)
    ),
    d, rho, pwindow
  );
}
