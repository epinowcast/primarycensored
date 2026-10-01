/**
  * Check if the analytical solution is built from uniform primary terms
  * @ingroup analytical_solution_helpers
  *
  * These are the delays whose censored CDF with a uniform primary is
  * primarycensored_uniform_lcdf_from_terms() applied to
  * primarycensored_uniform_terms() at d and q. Some only hold for part of
  * the parameter space, see check_uniform_terms_params().
  *
  * @param dist_id Distribution identifier for the delay distribution
  * @param primary_id Distribution identifier for the primary distribution
  *
  * @return 1 if the solution is built from uniform primary terms, 0
  * otherwise
  */
int check_for_uniform_terms(int dist_id, int primary_id) {
  if (primary_id != 1) return 0;
  return dist_id == 1 || dist_id == 2 || dist_id == 3 || dist_id == 4 ||
    dist_id == 5 || dist_id == 9 || dist_id == 13 || dist_id == 16 ||
    dist_id == 18 || dist_id == 19 || dist_id == 21 || dist_id == 22;
}

/**
  * Check the uniform primary terms are valid for the delay parameters
  * @ingroup analytical_solution_helpers
  *
  * The Inverse gamma terms need a finite mean, `alpha > 1`, and the Inverse
  * chi-square (19) and Scaled inverse chi-square (22) need `nu > 2`.
  * Otherwise the numerical path is used.
  *
  * @param dist_id Distribution identifier for the delay distribution
  * @param params Array of distribution parameters
  *
  * @return 1 if the uniform primary terms are valid for the parameters, 0
  * otherwise
  */
int check_uniform_terms_params(int dist_id, array[] real params) {
  if (dist_id == 16) return params[1] > 1;
  if (dist_id == 19 || dist_id == 22) return params[1] > 2;
  return 1;
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
  // Delays with uniform primary terms and a Uniform primary, see
  // check_for_uniform_terms()
  if (check_for_uniform_terms(dist_id, primary_id)) return 1;
  // Keep this primary list in sync with `primary_lcdf`; see the note above.
  if (dist_id == 26 || dist_id == 27 || dist_id == 28) {
    return primary_id == 1 || primary_id == 2;
  }
  return 0; // No analytical solution for other combinations
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
  * Delays with support on the reals use q = d - w_P. Delays with one
  * positive antiderivative G of F_T have terms [log G(t), -inf].
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
  real log_cdf = log_diff_exp(log_A, log_B) - log(pwindow);
  // Rounding can put the upper tail a few ulp above 1. Not fmin, which
  // would turn a NaN into 0.
  return log_cdf > 0 ? 0 : log_cdf;
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
  * Compute the uniform primary terms at t for an Exponential delay
  * @ingroup analytical_solution_helpers
  *
  * The antiderivative of F_T is G(t) = t - 1 / rate + exp(-rate * t) / rate,
  * so the terms are [log G(t), -inf]. A series is used below `rate * t = 0.1`
  * where the difference cancels.
  *
  * @param t Time (d or q)
  * @param params Array of Exponential distribution parameters [rate]
  *
  * @return Vector [log G(t), -inf], both `-inf` for t <= 0
  */
vector primarycensored_exponential_uniform_terms(real t,
                                                 array[] real params) {
  if (t <= 0) {
    return rep_vector(negative_infinity(), 2);
  }
  real rate = params[1];
  real x = rate * t;
  real log_x_shortfall;
  if (x < 0.1) {
    // x - 1 + exp(-x) = x^2 / 2 * sum_k 2 (-x)^k / (k + 2)!
    real term = 1;
    real series = 1;
    for (k in 1:10) {
      term *= -x / (k + 2);
      series += term;
    }
    log_x_shortfall = 2 * log(x) - log2() + log(series);
  } else {
    log_x_shortfall = log_diff_exp(log(x), log1m_exp(-x));
  }
  return [log_x_shortfall - log(rate), negative_infinity()]';
}

/**
  * Compute the uniform primary terms at t for a Normal delay
  * @ingroup analytical_solution_helpers
  *
  * The positive antiderivative of F_T is
  * G(t) = sigma * (z * Phi(z) + phi(z)) with z = (t - mu) / sigma.
  *
  * @param t Time (d or q)
  * @param params Array of Normal distribution parameters [mu, sigma]
  *
  * @return Vector [log G(t), -inf]
  */
vector primarycensored_normal_uniform_terms(real t, array[] real params) {
  real mu = params[1];
  real sigma = params[2];
  real z = (t - mu) / sigma;
  real log_g;
  if (z < -10) {
    // z * Phi(z) and phi(z) cancel, so use the asymptotic series
    real y = inv(square(z));
    real term = 1;
    real series = 1;
    for (k in 1:20) {
      term *= -(2 * k + 1) * y;
      series += term;
    }
    log_g = std_normal_lpdf(z) + log(y) + log(series);
  } else {
    // erfc has an exact autodiff derivative, unlike std_normal_lcdf
    log_g = log(exp(std_normal_lpdf(z)) + z * 0.5 * erfc(-z / sqrt2()));
  }
  return [log(sigma) + log_g, negative_infinity()]';
}

/**
  * Compute the uniform primary terms at t for a Beta delay
  * @ingroup analytical_solution_helpers
  *
  * The mean is E = a / (a + b) and the partial expectation distribution is
  * the Beta with shape1 `a + 1`. For t >= 1 both CDFs are 1.
  *
  * @param t Time (d or q)
  * @param params Array of Beta distribution parameters [a, b]
  *
  * @return Vector [log(t * F_T(t)), log(E * tilde F_T(t))], both `-inf` for
  * t <= 0
  */
vector primarycensored_beta_uniform_terms(real t, array[] real params) {
  if (t <= 0) {
    return rep_vector(negative_infinity(), 2);
  }
  real a = params[1];
  real b = params[2];
  real log_E = log(a) - log(a + b);
  if (t >= 1) {
    return [log(t), log_E]';
  }
  return [
    log(t) + beta_lcdf(t | a, b),
    log_E + beta_lcdf(t | a + 1, b)
  ]';
}

/**
  * Compute the log CDF of an Inverse gamma delay
  * @ingroup delay_log_cdfs
  *
  * The CDF is Q(alpha, x) with x = beta / y. It is evaluated on the log scale
  * from the series for P(alpha, x) when x <= alpha + 1 and the continued
  * fraction for Q(alpha, x) (modified Lentz) above. Unlike `inv_gamma_lcdf`
  * it does not underflow in the lower tail. At an integer `alpha` the
  * fraction terminates, so extra iterations keep the shape gradient exact.
  *
  * @param y Positive value at which the log CDF is evaluated
  * @param alpha Shape parameter
  * @param beta Scale parameter
  *
  * @return Log CDF at y
  */
real primarycensored_inv_gamma_lcdf(real y, real alpha, real beta) {
  real x = beta / y;
  if (x > alpha + 1) {
    int near_int = abs(alpha - round(alpha)) < 1e-9;
    real b = x + 1 - alpha;
    real c = 1e300;
    real d = inv(b);
    real h = d;
    for (i in 1:2000) {
      real an = -i * (i - alpha);
      b += 2;
      d = an * d + b;
      d = inv(abs(d) < 1e-300 ? 1e-300 : d);
      c = b + an / c;
      c = abs(c) < 1e-300 ? 1e-300 : c;
      h *= d * c;
      if (abs(d * c - 1) < 1e-15 && (i > alpha + 50 || !near_int)) break;
    }
    return (alpha - 1) * log(x) - x - lgamma(alpha) + log(x) + log(h);
  }
  real term = 1;
  real series = 1;
  for (n in 1:2000) {
    term *= x / (alpha + n);
    series += term;
    if (term < 1e-17 * series) break;
  }
  return log1m_exp(alpha * log(x) - x - lgamma(alpha + 1) + log(series));
}

/**
  * Compute the uniform primary terms at t for an Inverse gamma delay
  * @ingroup analytical_solution_helpers
  *
  * The mean is E = beta / (alpha - 1) and the partial expectation
  * distribution is the Inverse gamma with shape `alpha - 1`, so `alpha > 1`
  * is required.
  *
  * @param t Time (d or q)
  * @param params Array of Inverse gamma distribution parameters
  * [alpha, beta]
  *
  * @return Vector [log(t * F_T(t)), log(E * tilde F_T(t))], both `-inf` for
  * t <= 0
  */
vector primarycensored_invgamma_uniform_terms(real t, array[] real params) {
  real alpha = params[1];
  real beta = params[2];
  if (alpha <= 1) {
    reject("Inverse gamma uniform primary terms need shape > 1, shape is ",
           alpha);
  }
  if (t <= 0) {
    return rep_vector(negative_infinity(), 2);
  }
  real log_E = log(beta) - log(alpha - 1);
  return [
    log(t) + primarycensored_inv_gamma_lcdf(t | alpha, beta),
    log_E + primarycensored_inv_gamma_lcdf(t | alpha - 1, beta)
  ]';
}

/**
  * Compute the uniform primary terms at t for a Pareto delay
  * @ingroup analytical_solution_helpers
  *
  * The antiderivative of F_T from `y_min` is
  * G(t) = y_min * (delta + expm1(-(alpha - 1) * L) / (alpha - 1)), with
  * delta = t / y_min - 1 and L = log(t / y_min). It holds for every
  * `alpha > 0`. A series in L is used when `L * max(1, |alpha - 1|) < 0.1`
  * where the terms cancel.
  *
  * @param t Time (d or q)
  * @param params Array of Pareto distribution parameters [y_min, alpha]
  *
  * @return Vector [log G(t), -inf], both `-inf` for t <= y_min
  */
vector primarycensored_pareto_uniform_terms(real t, array[] real params) {
  real y_min = params[1];
  real alpha = params[2];
  if (t <= y_min) {
    return rep_vector(negative_infinity(), 2);
  }
  real L = log(t) - log(y_min);
  real e = alpha - 1;
  real log_G_over_y_min;
  if (L * fmax(1, abs(e)) < 0.1) {
    // sum_{n >= 2} L^n / n! * (1 + (-1)^n e^(n - 1))
    real power = L;
    real e_power = 1;
    real g = 0;
    for (n in 2:14) {
      power *= L / n;
      e_power *= e;
      g += power * (1 + (n % 2 == 0 ? e_power : -e_power));
    }
    log_G_over_y_min = log(g);
  } else {
    real u = -e * L;
    // expm1(u) / u, with a series near 0 so the derivative is finite
    real expm1_over_u = abs(u) < 1e-4
      ? 1 + u * (0.5 + u * (inv(6) + u * inv(24))) : expm1(u) / u;
    log_G_over_y_min = log(expm1(L) - L * expm1_over_u);
  }
  return [log(y_min) + log_G_over_y_min, negative_infinity()]';
}

/**
  * Compute the uniform primary terms at t for a delay distribution
  * @ingroup analytical_solution_helpers
  *
  * @param t Time (d or q)
  * @param dist_id Distribution identifier (1: Lognormal, 2: Gamma,
  *   3: Weibull, 4: Exponential, 5: Generalised gamma, 9: Beta,
  *   13: Chi-square, 16: Inverse gamma, 18: Normal, 19: Inverse chi-square,
  *   21: Pareto, 22: Scaled inverse chi-square), see
  *   check_for_uniform_terms()
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
  } else if (dist_id == 4) {
    return primarycensored_exponential_uniform_terms(t, params);
  } else if (dist_id == 18) {
    return primarycensored_normal_uniform_terms(t, params);
  } else if (dist_id == 16) {
    return primarycensored_invgamma_uniform_terms(t, params);
  } else if (dist_id == 13) {
    // Gamma with shape nu / 2 and rate 1 / 2
    return primarycensored_gamma_uniform_terms(t, {params[1] / 2, 0.5});
  } else if (dist_id == 19) {
    // Inverse gamma with shape nu / 2 and scale 1 / 2
    return primarycensored_invgamma_uniform_terms(
      t, {params[1] / 2, 0.5}
    );
  } else if (dist_id == 22) {
    // Inverse gamma with shape nu / 2 and scale nu * s^2 / 2
    return primarycensored_invgamma_uniform_terms(
      t, {params[1] / 2, params[1] * square(params[2]) / 2}
    );
  } else if (dist_id == 21) {
    return primarycensored_pareto_uniform_terms(t, params);
  } else if (dist_id == 9) {
    return primarycensored_beta_uniform_terms(t, params);
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
  // Delays on the reals are not clipped at 0
  real q = dist_has_positive_support(dist_id) ? max({d - pwindow, 0.0})
                                              : d - pwindow;

  if (dist_id == 2 && primary_id == 1) {
    return primarycensored_gamma_uniform_lcdf(d | q, params, pwindow);
  } else if (dist_id == 1 && primary_id == 1) {
    return primarycensored_lognormal_uniform_lcdf(d | q, params, pwindow);
  } else if (dist_id == 3 && primary_id == 1) {
    return primarycensored_weibull_uniform_lcdf(d | q, params, pwindow);
  } else if (dist_id == 5 && primary_id == 1) {
    return primarycensored_gengamma_uniform_lcdf(d | q, params, pwindow);
  } else if (check_for_uniform_terms(dist_id, primary_id)) {
    // The remaining delays with uniform primary terms
    if (dist_id == 9 && q >= 1) return 0;
    return primarycensored_uniform_lcdf_from_terms(
      primarycensored_uniform_terms(d, dist_id, params),
      primarycensored_uniform_terms(q, dist_id, params), pwindow
    );
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

  // Apply truncation normalization. Skip when F(L) = 0 makes it a no-op
  // (positive support, L <= 0), as in primarycensored_cdf().
  if (!is_inf(D) || L > 0 ||
      (!is_inf(L) && !dist_has_positive_support(dist_id))) {
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
  * The analytical uniform primary CDF at d combines terms at d and at
  * q = max(d - pwindow, 0) (d - pwindow for delays on the reals). With
  * an integer pwindow q is an integer delay too, so
  * primarycensored_analytical_lcdf_vectorized() can compute the terms once
  * per delay and share them. This needs the analytical solutions built from
  * primarycensored_uniform_terms(), see check_for_uniform_terms(). The
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
  return check_for_uniform_terms(dist_id, primary_id) &&
    pwindow >= 1 && floor(pwindow) == pwindow;
}

/**
  * Compute the primary event censored log CDF analytically at integer delays
  * @ingroup primary_event_analytical_distributions
  *
  * The log CDF at d combines the terms at d and at q (see
  * primarycensored_uniform_lcdf_from_terms()). Both are integer delays, so
  * the terms are computed once per delay and used for both, halving the CDF
  * evaluations. The values are the same as from
  * primarycensored_analytical_lcdf() at each delay without truncation.
  * Only for cases where check_for_analytical_vectorized() and
  * check_uniform_terms_params() are 1.
  *
  * @param start First delay to compute
  * @param n Last delay to compute, and the length of the result
  * @param dist_id Distribution identifier
  * @param params Array of distribution parameters
  * @param pwindow Primary event window, a positive integer
  *
  * @return Vector whose element d is the log CDF at d, for d in start:n.
  * Elements before start are not computed.
  */
vector primarycensored_analytical_lcdf_vectorized(data int start,
                                                  data int n,
                                                  data int dist_id,
                                                  array[] real params,
                                                  data real pwindow) {
  int pw = to_int(pwindow);
  vector[n] log_cdfs;
  int first = dist_has_positive_support(dist_id) ? max(start - pw, 0)
                                                 : start - pw;
  // terms[t - first + 1] holds the terms at delay t
  array[n - first + 1] vector[2] terms;
  for (t in first:n) {
    terms[t - first + 1] = primarycensored_uniform_terms(t, dist_id, params);
  }
  for (d in start:n) {
    int q = dist_has_positive_support(dist_id) ? max(d - pw, 0) : d - pw;
    log_cdfs[d] = (dist_id == 9 && q >= 1) ? 0
      : primarycensored_uniform_lcdf_from_terms(
          terms[d - first + 1], terms[q - first + 1], pwindow
        );
  }
  return log_cdfs;
}
