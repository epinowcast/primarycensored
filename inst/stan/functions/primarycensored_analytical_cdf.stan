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
  return dist_id == 2 || dist_id == 1 || dist_id == 3 || dist_id == 5
         || dist_id == 29;
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
  // Gamma, Lognormal, Weibull, generalised gamma and log-logistic with a
  // Uniform primary
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
  * delay exists, see check_for_tilt_transform(). The log-logistic solutions
  * need a shape of at least 0.01. It chooses the path in
  * primarycensored_cdf() and primarycensored_lcdf().
  *
  * @param dist_id Distribution identifier for the delay distribution
  * @param params Array of delay distribution parameters
  * @param primary_id Distribution identifier for the primary distribution
  * @param primary_params Array of primary distribution parameters
  * @param pwindow Primary event window
  *
  * @return 1 if the analytical solution applies, 0 if the numerical path is
  * needed
  */
int check_for_analytical_params(int dist_id, array[] real params,
                                int primary_id, array[] real primary_params,
                                data real pwindow) {
  if (check_for_exptilt(dist_id, primary_id)) {
    return check_for_tilt_transform(dist_id, -primary_params[1], params);
  }
  // The log-logistic partial expectation needs a shape of at least 0.01
  if (dist_id == 29 && primary_id == 1) {
    return params[2] >= 0.01;
  }
  return check_for_analytical(dist_id, primary_id);
}

/**
  * Check if the analytical solution is the path to use at a delay
  * @ingroup analytical_solution_helpers
  *
  * This is check_for_analytical_params() and, for the exponentially tilted
  * solutions, that the tilt transform is the path to use at d, see
  * check_for_tilt_transform_at(). It chooses the path in
  * primarycensored_cdf(), primarycensored_lcdf() and
  * primarycensored_lcdf_vectorized().
  *
  * @param dist_id Distribution identifier for the delay distribution
  * @param params Array of delay distribution parameters
  * @param primary_id Distribution identifier for the primary distribution
  * @param primary_params Array of primary distribution parameters
  * @param d Delay, the largest delay if there are several
  * @param pwindow Primary event window
  *
  * @return 1 if the analytical solution is the path to use at d, 0 if the
  * numerical path is needed
  */
int check_for_analytical_delay(int dist_id, array[] real params,
                               int primary_id, array[] real primary_params,
                               data real d, data real pwindow) {
  if (!check_for_analytical_params(dist_id, params, primary_id,
                                   primary_params, pwindow)) {
    return 0;
  }
  if (check_for_exptilt(dist_id, primary_id)) {
    return check_for_tilt_transform_at(
      dist_id, -primary_params[1], params, d, pwindow
    );
  }
  return 1;
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
  * Both terms come from one `gamma_lcdf_logx_pair()` evaluation. The
  * subtraction in `primarycensored_uniform_lcdf_from_terms()` amplifies their
  * rounding error by about d / pwindow. Rejects a rate that is not positive
  * and finite.
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
  real log_rate = gamma_log_rate(rate);
  // log E where E = k * theta = shape / rate is the mean of the delay
  real log_E = log(shape) - log_rate;
  // Both CDFs share one series or fraction evaluation at log(rate * t)
  real log_t = log(t);
  vector[2] log_F_T = gamma_lcdf_logx_pair(log_t + log_rate, shape);
  return [log_t + log_F_T[1], log_E + log_F_T[2]]';
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
  * The subtraction in `primarycensored_uniform_lcdf_from_terms()` loses
  * precision for large `k * shape`, with errors above 1e-6 possible once it
  * exceeds about 12000.
  * The two terms use separate `gamma_lcdf_logx()` calls, since the shapes
  * k and k + 1 / shape do not differ by one.
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
  real log_t = log(t);
  real log_x = shape * (log_t - log(scale));
  return [
    log_t + gamma_lcdf_logx(log_x, k),
    log_E + gamma_lcdf_logx(log_x, k_shift)
  ]';
}

// Log-logistic delay. With scale s and shape k, F(x) = 1 / (1 + (x / s)^(-k))
// and the partial expectation is int_0^t u f(u) du = t F(t) r_{1/k}(A) for
// A = (t / s)^k and r_a(A) = (1 + A) / A^(a + 1) int_0^A y^a / (1 + y)^2 dy.
// The ratio r_a is computed from three series that do not cancel, as the R
// function `.loglogistic_ratio()`. The parameters are [scale, shape], the
// order of Stan's loglogistic_cdf().

/**
  * Compute the log CDF of the log-logistic distribution
  * @ingroup delay_log_cdfs
  *
  * @param y Value at which to evaluate the log CDF (y > 0)
  * @param scale Scale parameter
  * @param shape Shape parameter
  *
  * @return Log CDF of the log-logistic distribution
  */
real primarycensored_loglogistic_lcdf(real y, real scale, real shape) {
  return -log1p_exp(-shape * (log(y) - log(scale)));
}

/**
  * Compute the log survival of the log-logistic distribution
  * @ingroup delay_log_cdfs
  *
  * @param y Value at which to evaluate the log survival (y > 0)
  * @param scale Scale parameter
  * @param shape Shape parameter
  *
  * @return Log survival of the log-logistic distribution
  */
real primarycensored_loglogistic_lccdf(real y, real scale, real shape) {
  return -log1p_exp(shape * (log(y) - log(scale)));
}

/**
  * Partial moment ratio of the log-logistic for A up to 1
  * @ingroup analytical_solution_helpers
  *
  * @param a Power, non-negative
  * @param log_A Log of A, at most 0
  * @param tol Relative tolerance of the series
  *
  * @return r_a(A)
  */
real loglogistic_ratio_small(real a, real log_A, real tol) {
  real A = exp(log_A);
  real B = A / (1 + A);
  real weight = exp(-a * log1p(A));
  real total = weight / (a + 1);
  for (m in 0:5000) {
    weight *= B * (a + m) / (m + 1);
    real term = weight / (a + m + 2);
    total += term;
    if (term <= tol * total) break;
  }
  return total;
}

/**
  * Partial moment M_a(1) of the log-logistic from the digamma function
  * @ingroup analytical_solution_helpers
  *
  * @param a Power, non-negative
  *
  * @return M_a(1)
  */
real loglogistic_partial_at_one(real a) {
  if (a <= 0) return 0.5;
  return 0.5 * (a * (digamma((a + 1) / 2) - digamma(a / 2)) - 1);
}

/**
  * Power series of M_a(A) - M_a(1) in tau = (A - 1) / (A + 1)
  * @ingroup analytical_solution_helpers
  *
  * The recurrence runs on the terms times tau^j so a large power does not
  * overflow.
  *
  * @param a Power, non-negative
  * @param tau Point, at most 1 / 2
  * @param tol Relative tolerance of the series
  *
  * @return M_a(A) - M_a(1)
  */
real loglogistic_tau_integral(real a, real tau, real tol) {
  real shift = 2 * a * tau;
  real tau_sq = square(tau);
  real previous = 0;
  real current = 1;
  real total = 1;
  for (j in 0:5000) {
    real following = (shift * current + (j - 1) * tau_sq * previous)
                     / (j + 1);
    previous = current;
    current = following;
    total += current / (j + 2);
    if (current <= tol * total) break;
  }
  return 0.5 * tau * total;
}

/**
  * Term of the tail series of the log-logistic partial moment ratio
  * @ingroup analytical_solution_helpers
  *
  * A^(-a) int_3^A y^(a - m - 2) dy, finite at a = m + 1.
  *
  * @param a Power, non-negative
  * @param m Term of the series in 1 / y
  * @param log_A Log of the upper limit, above log(3)
  * @param shrink (3 / A)^a, passed so it is computed once
  *
  * @return The term
  */
real loglogistic_tail_term(real a, int m, real log_A, real shrink) {
  real ell = log_A - log(3);
  real x = a - m - 1;
  real z = x * ell;
  if (abs(z) < 1) {
    // Series for (1 - exp(-z)) / z, to keep its derivative in a
    real exprel = abs(z) < 1e-2
                  ? 1 - z / 2 + square(z) / 6 - z * square(z) / 24
                    + square(square(z)) / 120 - z * square(square(z)) / 720
                  : -expm1(-z) / z;
    return exp(-(m + 1) * log_A) * ell * exprel;
  }
  return (exp(-(m + 1) * log_A)
          - shrink * pow(3.0, -(m + 1))) / x;
}

/**
  * Partial moment ratio of the log-logistic for A above the tail start
  * @ingroup analytical_solution_helpers
  *
  * @param a Power, non-negative
  * @param log_A Log of A, above log(3)
  * @param tol Relative tolerance of the series
  *
  * @return r_a(A)
  */
real loglogistic_ratio_large(real a, real log_A, real tol) {
  // M_a(3) / 3^a, the integral up to the start of the tail series
  real at_start = (loglogistic_partial_at_one(a)
                   + loglogistic_tau_integral(a, 0.5, tol))
                  * exp(-a * log(3));
  real shrink = exp(-a * (log_A - log(3)));
  real total = shrink * at_start;
  for (m in 0:200) {
    real term = (m + 1) * loglogistic_tail_term(a, m, log_A, shrink);
    total += m % 2 == 0 ? term : -term;
    if (abs(term) <= tol * abs(total)) break;
  }
  return (1 + exp(-log_A)) * total;
}

/**
  * Partial moment ratio r_a(A) of the log-logistic
  * @ingroup analytical_solution_helpers
  *
  * @param a Power, non-negative
  * @param log_A Log of A
  * @param tol Relative tolerance of the series
  *
  * @return r_a(A), at most 1 / (a + 1)
  */
real loglogistic_moment_ratio(real a, real log_A, real tol) {
  if (log_A <= 0) {
    return loglogistic_ratio_small(a, log_A, tol);
  }
  if (log_A > log(3)) {
    return loglogistic_ratio_large(a, log_A, tol);
  }
  // 1 < A <= 3
  real A = exp(log_A);
  real tau = (A - 1) / (A + 1);
  return (loglogistic_partial_at_one(a) + loglogistic_tau_integral(a, tau, tol))
         * (1 + A) / A * exp(-a * log_A);
}

/**
  * Rounding error scale of a value held on the log scale
  * @ingroup analytical_solution_helpers
  *
  * @param x Value on the log scale
  *
  * @return x + log(1 + |x|), `-inf` for `-inf`
  */
real loglogistic_log_error_scale(real x) {
  if (x == negative_infinity()) return negative_infinity();
  return x + log1p(abs(x));
}

/**
  * Check if a window difference of the partial expectation is ill conditioned
  * @ingroup analytical_solution_helpers
  *
  * The CDF is a difference of t F(t) terms at d and d - pwindow over pwindow.
  * The bound on its error is 4e-16 times the two scaled terms. It is too
  * large where it exceeds 1e-8 of the smaller of the CDF and the survival,
  * not taken below 1e-15 in the upper tail, or 1e-6 of the density of the
  * CDF, not taken below 1e-13 of the smaller tail. The R equivalent is
  * `.loglogistic_uniform_ill()`.
  *
  * @param log_t_cdf_d Log of d F(d), see
  * primarycensored_loglogistic_uniform_terms()
  * @param log_t_cdf_q Log of q F(q), `-inf` for q <= 0
  * @param d Delay
  * @param pwindow Primary event window
  * @param log_cdf Log CDF at d from the solution, NaN where it failed
  *
  * @return 1 if the numerical CDF is needed, 0 otherwise
  */
int loglogistic_window_ill_conditioned(real log_t_cdf_d, real log_t_cdf_q,
                                       real d, data real pwindow,
                                       real log_cdf) {
  if (is_nan(log_cdf)) return 1;
  if (log_t_cdf_d == negative_infinity()) return 0;
  real log_error = log(4e-16) - log(pwindow) + log_sum_exp(
    loglogistic_log_error_scale(log_t_cdf_d),
    loglogistic_log_error_scale(log_t_cdf_q)
  );
  real log_cdf_q = log_t_cdf_q == negative_infinity()
                   ? negative_infinity() : log_t_cdf_q - log(d - pwindow);
  real log_density = primarycensored_log_diff_exp(
    log_t_cdf_d - log(d), log_cdf_q
  ) - log(pwindow);
  real log_capped = fmin(log_cdf, 0);
  real log_tail = fmin(log_capped, log1m_exp(log_capped));
  real log_limit = log(1e-8) + log_tail;
  if (log_capped > -log2()) {
    log_limit = fmax(log_limit, log(1e-15));
  }
  log_limit = fmin(
    log_limit, fmax(log(1e-6) + log_density, log(1e-13) + log_tail)
  );
  return log_error > log_limit;
}

/**
  * ODE system for the numerical log-logistic primary event censored CDF
  * @ingroup ode
  *
  * The integrand is the delay CDF, or the survival, at d - x, integrated
  * from x = 0 so the resolution does not fall as d grows. It is divided by
  * its bound so the absolute tolerance is relative to the integral. The
  * bound is the delay CDF at d for the CDF.
  *
  * @param x Primary event time
  * @param y State, the scaled integral
  * @param theta Parameters [scale, shape, log of the scale of the integrand]
  * @param x_r Real data [d]
  * @param x_i Integer data [1 to integrate the survival]
  *
  * @return Derivative of the state
  */
vector loglogistic_numeric_ode(real x, vector y, array[] real theta,
                               array[] real x_r, array[] int x_i) {
  real u = x_r[1] - x;
  real log_delay;
  if (u <= 0) {
    log_delay = x_i[1] ? 0 : negative_infinity();
  } else if (x_i[1]) {
    log_delay = primarycensored_loglogistic_lccdf(u | theta[1], theta[2]);
  } else {
    log_delay = primarycensored_loglogistic_lcdf(u | theta[1], theta[2]);
  }
  return rep_vector(exp(log_delay - theta[3]), 1);
}

/**
  * Log of the numerical log-logistic primary event censored CDF
  * @ingroup analytical_solution_helpers
  *
  * The CDF used where the uniform solution is not accurate. It integrates
  * the delay CDF over the primary event time with an ODE solver, see
  * loglogistic_numeric_ode(), with a relative tolerance of 1e-13 and an
  * absolute tolerance of 1e-20 for the scaled integrand. The survival is
  * integrated where the delay CDF at the window midpoint is above a half, so
  * the upper tail keeps its relative precision.
  *
  * @param d Delay
  * @param params Array of log-logistic parameters [scale, shape]
  * @param pwindow Primary event window
  *
  * @return Log of the numerical primary event censored CDF with a uniform
  * primary, at most 0, `-inf` for d <= 0
  */
real loglogistic_numeric_lcdf(data real d, array[] real params,
                              data real pwindow) {
  if (d <= 0) return negative_infinity();
  real x_end = fmin(d, pwindow);
  real mid = d - 0.5 * pwindow;
  int use_survival = mid > 0
              && primarycensored_loglogistic_lcdf(mid | params[1], params[2])
                 > -log2();
  real log_scale = use_survival
                   ? 0
                   : primarycensored_loglogistic_lcdf(d | params[1], params[2]);
  real integral = ode_rk45_tol(
    loglogistic_numeric_ode, rep_vector(0.0, 1), 0.0, {x_end}, 1e-13, 1e-20,
    1000000, {params[1], params[2], log_scale}, {d}, {use_survival}
  )[1, 1];
  if (use_survival) {
    real survival = integral / pwindow;
    if (d < pwindow) {
      survival += 1 - d / pwindow;
    }
    return survival >= 1 ? negative_infinity() : log1m(survival);
  }
  return integral > 0 ? fmin(log(integral) + log_scale - log(pwindow), 0)
                      : negative_infinity();
}

/**
  * Compute the uniform primary terms at t for a log-logistic delay
  * @ingroup analytical_solution_helpers
  *
  * The partial expectation is t F(t) r_{1/shape}(A), see
  * loglogistic_moment_ratio().
  *
  * @param t Time (d or q)
  * @param params Array of log-logistic distribution parameters
  * [scale, shape]
  *
  * @return Vector [log(t * F_T(t)), log(E * tilde F_T(t))], both `-inf` for
  * t <= 0
  */
vector primarycensored_loglogistic_uniform_terms(real t,
                                                 array[] real params) {
  if (t <= 0) {
    return rep_vector(negative_infinity(), 2);
  }
  real shape = params[2];
  real log_ratio = shape * (log(t) - log(params[1]));
  real log_t_F_T = log(t) - log1p_exp(-log_ratio);
  return [
    log_t_F_T,
    log_t_F_T + log(loglogistic_moment_ratio(1 / shape, log_ratio, 1e-16))
  ]';
}

/**
  * Compute the uniform primary terms at t for a delay distribution
  * @ingroup analytical_solution_helpers
  *
  * @param t Time (d or q)
  * @param dist_id Distribution identifier (1: Lognormal, 2: Gamma,
  *   3: Weibull, 5: Generalised gamma, 29: Log-logistic), see
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
  } else if (dist_id == 29) {
    return primarycensored_loglogistic_uniform_terms(t, params);
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
  * Compute the primary event censored log CDF analytically for log-logistic
  * delay with Uniform primary
  * @ingroup primary_event_analytical_distributions
  *
  * @param d Delay time
  * @param q Lower bound of integration (max(d - pwindow, 0))
  * @param params Array of log-logistic distribution parameters
  * [scale, shape]
  * @param pwindow Primary event window
  *
  * @return Log of the primary event censored CDF for log-logistic delay with
  * Uniform primary, from loglogistic_numeric_lcdf() where the solution is
  * ill conditioned
  */
real primarycensored_loglogistic_uniform_lcdf(data real d, real q,
                                              array[] real params,
                                              data real pwindow) {
  vector[2] terms_d = primarycensored_loglogistic_uniform_terms(d, params);
  vector[2] terms_q = primarycensored_loglogistic_uniform_terms(q, params);
  real log_cdf = primarycensored_uniform_lcdf_from_terms(
    terms_d, terms_q, pwindow
  );
  if (loglogistic_window_ill_conditioned(terms_d[1], terms_q[1], d, pwindow,
                                         log_cdf)) {
    return loglogistic_numeric_lcdf(d | params, pwindow);
  }
  return fmin(log_cdf, 0);
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
  } else if (dist_id == 29 && primary_id == 1) {
    return primarycensored_loglogistic_uniform_lcdf(d | q, params, pwindow);
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
  * there. Where that uses the ODE, its values are shifted by their difference
  * from the analytical value at the delay before, so that the PMF across the
  * switch is a difference of ODE values.
  * The log-logistic uses the numerical CDF where the solution is ill
  * conditioned, as primarycensored_loglogistic_uniform_lcdf() does.
  * Only for cases where check_for_analytical_vectorized() is 1 and
  * check_for_analytical_delay() is 1 at n.
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
    int n_endpoints = n - first + 1;
    // Each delay uses the form of primarycensored_exptilt_form(), with each
    // endpoint computed once
    array[n] int form;
    array[n_endpoints] int needs_moments = rep_array(0, n_endpoints);
    array[n_endpoints] int needs_terms = rep_array(0, n_endpoints);
    for (d in start:n) {
      int d_index = d - first + 1;
      int q_index = (positive ? max(d - pw, 0) : d - pw) - first + 1;
      form[d] = primarycensored_exptilt_form(d, dist_id, rho, pwindow);
      if (form[d] == 1) {
        needs_moments[d_index] = 1;
        needs_moments[q_index] = 1;
      } else if (form[d] == 2) {
        needs_moments[d_index] = 1;
      } else {
        needs_terms[d_index] = 1;
        needs_terms[q_index] = 1;
      }
    }
    // moments and terms at endpoint t are at index t - first + 1
    array[n_endpoints] vector[3] moments;
    array[n_endpoints] vector[5] terms;
    vector[4] context = rep_vector(0, 4);
    if (sum(needs_terms) > 0) {
      context = log_tilt_transform_context(dist_id, -rho, params);
    }
    for (t in first:n) {
      int t_index = t - first + 1;
      if (needs_moments[t_index]) {
        moments[t_index] = primarycensored_tilt_moments(t, dist_id, params);
      }
      if (needs_terms[t_index]) {
        terms[t_index] = primarycensored_exptilt_terms(
          t, dist_id, rho, params, context
        );
      }
    }
    // Offset of the ODE CDF from the analytical CDF at the delay before
    real shift = 0;
    int is_ode = 1;
    for (d in start:n) {
      int d_index = d - first + 1;
      int q_index = (positive ? max(d - pw, 0) : d - pw) - first + 1;
      if (form[d] == 1) {
        log_cdfs[d] = primarycensored_exptilt_small_window_lcdf_from_terms(
          moments[d_index], moments[q_index], rho, pwindow
        );
        is_ode = 0;
      } else if (form[d] == 2) {
        log_cdfs[d] = primarycensored_exptilt_small_delay_lcdf_from_terms(
          moments[d_index], rho, pwindow
        );
        is_ode = 0;
      } else {
        log_cdfs[d] = primarycensored_exptilt_lcdf_from_terms(
          terms[d_index], terms[q_index], d, rho, pwindow
        );
        if (is_nan(log_cdfs[d])) {
          if (!is_ode) {
            shift = exp(log_cdfs[d - 1]) - primarycensored_numeric_cdf(
              d - 1 | dist_id, params, pwindow, primary_id, primary_params
            );
          }
          log_cdfs[d] = log(shift + primarycensored_numeric_cdf(
            d | dist_id, params, pwindow, primary_id, primary_params
          ));
          is_ode = 1;
        } else {
          is_ode = 0;
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
  // The numerical CDF can be out of order by its tolerance, which would make
  // the PMF NaN
  int previous_numerical = 0;
  for (d in start:n) {
    vector[2] terms_q = terms[max(d - pw, 0) + 1];
    log_cdfs[d] = primarycensored_uniform_lcdf_from_terms(
      terms[d + 1], terms_q, pwindow
    );
    if (dist_id == 29) {
      int numerical = loglogistic_window_ill_conditioned(
        terms[d + 1][1], terms_q[1], d, pwindow, log_cdfs[d]
      );
      if (numerical) {
        log_cdfs[d] = loglogistic_numeric_lcdf(d | params, pwindow);
      }
      log_cdfs[d] = fmin(log_cdfs[d], 0);
      if (d > start && (numerical || previous_numerical)) {
        log_cdfs[d] = fmax(log_cdfs[d], log_cdfs[d - 1]);
      }
      previous_numerical = numerical;
    }
  }
  return log_cdfs;
}

// Tilt transforms T_f(xi; tau) = int exp(xi u) f(u) du up to tau, from 0 for
// delays on the non-negative reals. The R equivalents are `.pcens_tilt_*()`.

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
  * Log of the incomplete gamma series for the tilt transform of a Stacy
  * family delay
  * @ingroup tilt_transforms
  *
  * The Weibull (k = 1) and generalised gamma have CDF P(k, (t / scale)^shape)
  * and transform
  *   T_f(xi; t) = sum_n (xi scale)^n / n! Gamma(k + n / shape) / Gamma(k)
  *                P(k + n / shape, (t / scale)^shape).
  * The terms for negative xi alternate in sign, and the sum of their absolute
  * values is the transform at -xi. Positive and negative terms are summed
  * apart until the remainder, bounded by a geometric series once
  * n + 1 > |xi| t, is below 1e-17 of that sum.
  *
  * @param t Point, positive
  * @param xi Tilt, not zero
  * @param shape Power of the Stacy family
  * @param scale Scale
  * @param k Shape of the underlying gamma
  *
  * @return Vector [log T_f(xi; t), log of the sum of the absolute terms over
  *   the result], or `nan` for both where the sum of the absolute terms is
  *   more than 1e5 times the result, |xi| t is 400 or more, or 400 terms do
  *   not converge
  */
vector primarycensored_stacy_series(real t, real xi, real shape, real scale,
                                    real k) {
  real log_x = shape * (log(t) - log(scale));
  int alternating = xi < 0;
  real log_z = log(abs(xi)) + log(scale);
  real log_w = log(abs(xi)) + log(t);
  real w = exp(log_w);
  vector[2] unreliable = rep_vector(not_a_number(), 2);
  if (w >= 400) return unreliable;
  real log_k = lgamma(k);
  real log_pos = 0;
  real log_neg = negative_infinity();
  real log_first = 0;
  for (n in 0:400) {
    real s = k + n / shape;
    real log_term = n * log_z - lgamma(n + 1) + lgamma(s) - log_k
                    + gamma_lcdf_logx(log_x, s);
    if (n == 0) {
      log_first = log_term;
      log_pos = log_term;
    } else if (alternating && n % 2 == 1) {
      log_neg = log_sum_exp(log_neg, log_term);
    } else {
      log_pos = log_sum_exp(log_pos, log_term);
    }
    real log_abs = log_sum_exp(log_pos, log_neg);
    // The transform is at most the CDF for negative xi
    if (alternating && log_abs - log_first > log(1e5)) return unreliable;
    if (n + 1 > w
        && log_term + log_w - log(n + 1 - w) < log_abs + log(1e-17)) {
      if (log_pos <= log_neg) return unreliable;
      real log_sum = primarycensored_log_diff_exp(log_pos, log_neg);
      if (log_abs - log_sum > log(1e5)) return unreliable;
      return [log_sum, log_abs - log_sum]';
    }
  }
  return unreliable;
}

/**
  * Check if the tilt transform is closed form for a delay and tilt
  * @ingroup tilt_transforms
  *
  * The exponential (4) and gamma (2) forms need the delay with the rate
  * lowered by xi, which exists if rate - xi > 0. The normal (18) form and the
  * Weibull (3) and generalised gamma (5) series have no restriction. The
  * lognormal (1) form needs xi sigma^2 exp(mu) below exp(690) in magnitude.
  * Callers use the numerical path when this is 0.
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
  if (dist_id == 18 || dist_id == 3 || dist_id == 5) return 1;
  if (dist_id == 1) {
    if (params[2] <= 0) return 0;
    if (xi >= 0) return 1;
    return log(-xi) + 2 * log(params[2]) + params[1] < 690;
  }
  return 0;
}

/**
  * Check if the tilt transform is the path to use at a point
  * @ingroup tilt_transforms
  *
  * This is check_for_tilt_transform() and, for the lognormal (1) with
  * xi > 0, that the series is the path to use at t. The series costs about
  * xi t + 9 sqrt(xi t) + 30 terms and is limited to 20000 terms. It is used
  * up to xi t of 60, and beyond that only where xi w is above 2, as the ODE
  * is inaccurate in the lower tail there. The cut-off of 60 is set for the
  * ODE and differs from the one in R, which is set for the R numerical path.
  *
  * @param dist_id Distribution identifier for the delay distribution
  * @param xi Tilt. The exponentially tilted window with tilt rho needs
  *   xi = -rho
  * @param params Array of distribution parameters, as for dist_lcdf()
  * @param t Point, the largest of the points if there are several
  * @param pwindow Primary event window
  *
  * @return 1 if the transform is closed form and is the path to use at t,
  * 0 otherwise
  */
int check_for_tilt_transform_at(int dist_id, real xi, array[] real params,
                                data real t, data real pwindow) {
  if (!check_for_tilt_transform(dist_id, xi, params)) return 0;
  if (dist_id == 1 && xi > 0 && t > 0) {
    if (xi * t + 9 * sqrt(xi * t) + 30 > 20000) return 0;
    return xi * t <= 60 || xi * pwindow > 2;
  }
  return 1;
}

/**
  * Nodes of the 32 point Gauss-Legendre rule on [-1, 1]
  * @ingroup tilt_transforms
  *
  * @return The 16 positive nodes, the others are their negatives
  */
vector primarycensored_gauss_legendre_nodes() {
  return [
    0.048307665687738324, 0.14447196158279649, 0.23928736225213706,
    0.33186860228212767, 0.42135127613063533, 0.50689990893222936,
    0.5877157572407623, 0.66304426693021523, 0.73218211874028971,
    0.79448379596794239, 0.84936761373256997, 0.89632115576605209,
    0.93490607593773967, 0.96476225558750639, 0.98561151154526838,
    0.99726386184948157
  ]';
}

/**
  * Log weights of the 32 point Gauss-Legendre rule on [-1, 1]
  * @ingroup tilt_transforms
  *
  * @return The log weights of the 16 positive nodes, shared by the negatives
  */
vector primarycensored_gauss_legendre_log_weights() {
  return [
    -2.3377969318791947, -2.3471775191742132, -2.3661171972104831,
    -2.3949868407368724, -2.4343797887852916, -2.4851635920714812,
    -2.5485636934665021, -2.626297960160425, -2.7207977655059796,
    -2.8355865687925124, -2.9759677006607701, -3.1503787891021764,
    -3.3733722294882642, -3.6733185431605069, -4.1181622816923937,
    -4.9591900848988502
  ]';
}

/**
  * Log of a Gaussian weighted integral of the tilt over one panel
  * @ingroup tilt_transforms
  *
  * The log of the integral of exp(-rho exp(mu + sigma z)) phi(z) over [a, b]
  * by the 32 point Gauss-Legendre rule, where phi is the standard normal
  * density.
  *
  * @param a Lower limit
  * @param b Upper limit
  * @param mu Mean of the log of the delay
  * @param sigma Standard deviation of the log of the delay
  * @param rho Tilt, at least 0
  *
  * @return Log of the integral, `-inf` if b <= a
  */
real primarycensored_lognormal_tilt_one_panel(real a, real b, real mu,
                                              real sigma, real rho) {
  if (b <= a) return negative_infinity();
  real half = 0.5 * (b - a);
  real mid = 0.5 * (a + b);
  vector[16] nodes = primarycensored_gauss_legendre_nodes();
  vector[16] log_weights = primarycensored_gauss_legendre_log_weights();
  vector[32] terms;
  for (i in 1:16) {
    real z_up = mid + half * nodes[i];
    real z_down = mid - half * nodes[i];
    terms[i] = log_weights[i] - rho * exp(mu + sigma * z_up)
               - 0.5 * square(z_up);
    terms[16 + i] = log_weights[i] - rho * exp(mu + sigma * z_down)
                    - 0.5 * square(z_down);
  }
  return log_sum_exp(terms) + log(half) - 0.5 * log(2 * pi());
}

/**
  * Log of a Gaussian weighted integral of the tilt over a range
  * @ingroup tilt_transforms
  *
  * As primarycensored_lognormal_tilt_one_panel(), with the range split into
  * ceil(sigma / 1.8) equal panels. One panel loses accuracy for sigma above
  * 1.8.
  *
  * @param a Lower limit
  * @param b Upper limit
  * @param mu Mean of the log of the delay
  * @param sigma Standard deviation of the log of the delay
  * @param rho Tilt, at least 0
  *
  * @return Log of the integral, `-inf` if b <= a
  */
real primarycensored_lognormal_tilt_panel(real a, real b, real mu,
                                          real sigma, real rho) {
  if (b <= a) return negative_infinity();
  int n_panels = 1;
  while (n_panels * 1.8 < sigma) {
    n_panels += 1;
  }
  if (n_panels == 1) {
    return primarycensored_lognormal_tilt_one_panel(a, b, mu, sigma, rho);
  }
  real width = (b - a) / n_panels;
  vector[n_panels] panels;
  for (j in 1:n_panels) {
    panels[j] = primarycensored_lognormal_tilt_one_panel(
      a + (j - 1) * width, a + j * width, mu, sigma, rho
    );
  }
  return log_sum_exp(panels);
}

/**
  * Mode, limits and total of the tilted lognormal integrand
  * @ingroup tilt_transforms
  *
  * The log integrand in z = (log u - mu) / sigma is
  * -rho exp(mu + sigma z) - z^2 / 2. It is concave with a single mode
  * z0 = -W0(rho sigma^2 exp(mu)) / sigma. Beyond the limits lo and hi it is
  * below exp(-40) of its peak. None of these depends on the point, so they
  * are computed once per call, see log_tilt_transform_context().
  *
  * @param mu Mean of the log of the delay
  * @param sigma Standard deviation of the log of the delay
  * @param rho Tilt, above 0
  *
  * @return Vector [z0, lo, hi, log total]
  */
vector primarycensored_lognormal_tilt_bump(real mu, real sigma, real rho) {
  real w0 = lambert_w0(exp(log(rho) + 2 * log(sigma) + mu));
  real z0 = -w0 / sigma;
  real u0 = exp(mu - w0);
  real hi = fmin(
    z0 + sqrt(80 / (1 + w0)), fmax((log(u0 + 40 / rho) - mu) / sigma, -z0)
  );
  real lo = fmax(
    z0 - sqrt(80), -sqrt(square(z0) + 2 * w0 / square(sigma) + 80)
  );
  real log_total = log_sum_exp(
    primarycensored_lognormal_tilt_panel(lo, z0, mu, sigma, rho),
    primarycensored_lognormal_tilt_panel(z0, hi, mu, sigma, rho)
  );
  return [z0, lo, hi, log_total]';
}

/**
  * Log of the lognormal tilt transform by quadrature for a negative tilt
  * @ingroup tilt_transforms
  *
  * The lower and upper transform for xi = -rho < 0. The transform on the
  * side of the mode that t is on is one panel that starts or ends at t. The
  * other is the difference from the total, which does not cancel as it is at
  * least the mass on the far side of the mode.
  *
  * @param t Point, positive
  * @param mu Mean of the log of the delay
  * @param sigma Standard deviation of the log of the delay
  * @param rho Tilt, above 0
  * @param bump Output of primarycensored_lognormal_tilt_bump()
  *
  * @return Vector [log T_f(-rho; t), log(T_f(-rho; Inf) - T_f(-rho; t))]
  */
vector primarycensored_lognormal_tilt_quadrature(real t, real mu,
                                                 real sigma, real rho,
                                                 vector bump) {
  real z = (log(t) - mu) / sigma;
  // The slope of the log integrand at t, positive left of the mode
  real slope = -rho * sigma * t - z;
  if (z < bump[1]) {
    real width = 80 / (slope + sqrt(square(slope) + 80));
    real log_lower = primarycensored_lognormal_tilt_panel(
      z - width, z, mu, sigma, rho
    );
    return [
      log_lower, primarycensored_log_diff_exp(bump[4], log_lower)
    ]';
  }
  real curvature = 1 + rho * square(sigma) * t;
  real width = 80 / (abs(slope) + sqrt(square(slope) + 80 * curvature));
  // The integrand is below exp(-40) of its value at t beyond the point where
  // the tilt alone has decayed by 40
  width = fmin(
    width, fmax(abs(z), (log(t + 40 / rho) - mu) / sigma) - z
  );
  real log_upper = primarycensored_lognormal_tilt_panel(
    z, z + width, mu, sigma, rho
  );
  return [primarycensored_log_diff_exp(bump[4], log_upper), log_upper]';
}

/**
  * Log of the lognormal tilt transform by series for a positive tilt
  * @ingroup tilt_transforms
  *
  * The lower transform for xi > 0 as the sum of the positive terms
  * xi^k m_k(t) / k! over the partial moments
  * m_k(t) = exp(k mu + k^2 sigma^2 / 2) Phi((log t - mu) / sigma - k sigma).
  * The sum stops once past k = xi t a term is below exp(-41) of the sum. The
  * total diverges, so there is no upper transform. It rejects past 20000
  * terms, see check_for_tilt_transform_at().
  *
  * @param t Point, positive
  * @param mu Mean of the log of the delay
  * @param sigma Standard deviation of the log of the delay
  * @param xi Tilt, above 0
  *
  * @return log T_f(xi; t)
  */
real primarycensored_lognormal_tilt_series(real t, real mu, real sigma,
                                           real xi) {
  real z = (log(t) - mu) / sigma;
  real log_xi_mu = log(xi) + mu;
  real log_factorial = 0;
  real total = primarycensored_log_std_normal_cdf(z);
  for (k in 1:20000) {
    log_factorial += log(k);
    real log_term = k * log_xi_mu + 0.5 * square(k * sigma) - log_factorial
                    + primarycensored_log_std_normal_cdf(z - k * sigma);
    total = log_sum_exp(total, log_term);
    if (k > xi * t && log_term < total - 41) return total;
  }
  reject(
    "The lognormal tilt transform needs more than 20000 terms for t = ", t,
    " and xi = ", xi, ". Use the numerical path."
  );
}

/**
  * Log of the lognormal tilt transform over the lower and upper support
  * @ingroup tilt_transforms
  *
  * The upper transform is `inf` for xi > 0, as the total diverges.
  *
  * @param t Point
  * @param mu Mean of the log of the delay
  * @param sigma Standard deviation of the log of the delay
  * @param xi Tilt
  * @param bump Output of primarycensored_lognormal_tilt_bump() for -xi,
  *   used for xi < 0
  *
  * @return Vector [log T_f(xi; t), log(T_f(xi; Inf) - T_f(xi; t))]
  */
vector primarycensored_lognormal_tilt_pair(real t, real mu, real sigma,
                                           real xi, vector bump) {
  if (xi == 0) {
    if (t <= 0) return [negative_infinity(), 0]';
    real z = (log(t) - mu) / sigma;
    return [
      primarycensored_log_std_normal_cdf(z),
      primarycensored_log_std_normal_cdf(-z)
    ]';
  }
  if (xi > 0) {
    if (t <= 0) return [negative_infinity(), positive_infinity()]';
    return [
      primarycensored_lognormal_tilt_series(t, mu, sigma, xi),
      positive_infinity()
    ]';
  }
  if (t <= 0) return [negative_infinity(), bump[4]]';
  return primarycensored_lognormal_tilt_quadrature(t, mu, sigma, -xi, bump);
}

/**
  * Quantities shared by the tilt transforms of a delay at every point
  * @ingroup tilt_transforms
  *
  * Computed once per call and passed to log_tilt_transform_pair_shared(). It
  * has length 4 and is 0 for delays whose transform has no shared setup.
  *
  * @param dist_id Distribution identifier, see check_for_tilt_transform()
  * @param xi Tilt
  * @param params Array of distribution parameters, as for dist_lcdf()
  *
  * @return Vector of 4 shared quantities
  */
vector log_tilt_transform_context(int dist_id, real xi, array[] real params) {
  if (dist_id == 1 && xi < 0) {
    return primarycensored_lognormal_tilt_bump(params[1], params[2], -xi);
  }
  return rep_vector(0, 4);
}

/**
  * Log of the tilt transform over the lower and the upper part of the support
  * @ingroup tilt_transforms
  *
  * The lower transform is T_f(xi; t), the delay CDF for xi = 0, and `-inf`
  * for t <= 0 for delays on the non-negative reals. The upper transform is
  * T_f(xi; Inf) - T_f(xi; t). The third element is the log of the factor by
  * which the terms of the lower transform exceed it, 0 unless the transform is
  * a series. Only defined where check_for_tilt_transform() is 1.
  *
  * The gamma is the gamma CDF with the rate lowered by xi times the total
  * (rate / (rate - xi))^shape, with both tails from gamma_log_pq_logx().
  *
  * The lognormal (1) has no closed form, see
  * primarycensored_lognormal_tilt_pair(). Its upper transform is `inf` for
  * xi > 0, where the total diverges.
  *
  * @param t Point
  * @param dist_id Distribution identifier: 1 (Lognormal), 2 (Gamma),
  *   3 (Weibull), 4 (Exponential), 5 (Generalised gamma) or 18 (Normal), see
  *   check_for_tilt_transform()
  * @param xi Tilt
  * @param params Array of distribution parameters, as for dist_lcdf()
  * @param context Output of log_tilt_transform_context() for dist_id, xi
  *   and params
  *
  * @return Vector [log T_f(xi; t), log(T_f(xi; Inf) - T_f(xi; t)), log loss].
  *   The Weibull and generalised gamma have `nan` for the lower transform
  *   where the series is not reliable, and for the upper transform if xi is
  *   not 0.
  */
vector log_tilt_transform_pair_shared(real t, int dist_id, real xi,
                                      array[] real params, vector context) {
  if (dist_id == 2) {
    real shape = params[1];
    real rate = params[2];
    real log_total = -shape * log1m(xi / rate);
    if (t <= 0) return [negative_infinity(), log_total, 0]';
    return append_row(
      log_total + gamma_log_pq_logx(log(t) + log(rate - xi), shape), 0
    );
  } else if (dist_id == 4) {
    real rate = params[1];
    real tilted_rate = rate - xi;
    real log_total = -log1m(xi / rate);
    if (t <= 0) return [negative_infinity(), log_total, 0]';
    return [
      log_total + log1m_exp(-tilted_rate * t), log_total - tilted_rate * t, 0
    ]';
  } else if (dist_id == 18) {
    // The upper tail is the lower tail of the reflected normal
    real mu = params[1];
    real sigma = params[2];
    real z = (t - mu - xi * square(sigma)) / sigma;
    real log_total = xi * mu + 0.5 * square(xi * sigma);
    return [
      log_total + primarycensored_log_std_normal_cdf(z),
      log_total + primarycensored_log_std_normal_cdf(-z), 0
    ]';
  } else if (dist_id == 3 || dist_id == 5) {
    real shape = params[1];
    real scale = params[2];
    real k = dist_id == 3 ? 1.0 : params[3];
    if (xi == 0) {
      if (t <= 0) return [negative_infinity(), 0, 0]';
      return append_row(
        gamma_log_pq_logx(shape * (log(t) - log(scale)), k), 0
      );
    }
    if (t <= 0) return [negative_infinity(), not_a_number(), 0]';
    vector[2] series = primarycensored_stacy_series(t, xi, shape, scale, k);
    return [series[1], not_a_number(), series[2]]';
  } else if (dist_id == 1) {
    vector[2] pair = primarycensored_lognormal_tilt_pair(
      t, params[1], params[2], xi, context
    );
    return append_row(pair, 0);
  }
  reject("Invalid distribution identifier: ", dist_id);
}

/**
  * Log of the tilt transform over the lower and the upper part of the support
  * @ingroup tilt_transforms
  *
  * As log_tilt_transform_pair_shared() with the shared quantities computed
  * for this point alone.
  *
  * @param t Point
  * @param dist_id Distribution identifier, see log_tilt_transform_pair_shared()
  * @param xi Tilt
  * @param params Array of distribution parameters, as for dist_lcdf()
  *
  * @return Vector [log T_f(xi; t), log(T_f(xi; Inf) - T_f(xi; t)), log loss]
  */
vector log_tilt_transform_pair(real t, int dist_id, real xi,
                               array[] real params) {
  return log_tilt_transform_pair_shared(
    t, dist_id, xi, params, log_tilt_transform_context(dist_id, xi, params)
  );
}

/**
  * Log moments about a point from log partial moments
  * @ingroup tilt_transforms
  *
  * The log of G_k(t) = int_0^t (t - u)^k f(u) du for k = 1, 2, 3 from the
  * partial moments m_k(t) = int_0^t u^k f(u) du, with G_1 = t m_0 - m_1,
  * G_2 = t G_1 - (t m_1 - m_2) and
  * G_3 = t G_2 - (t (t m_1 - m_2) - (t m_2 - m_3)), each a difference of
  * positive integrals.
  *
  * @param log_t Log of the point
  * @param log_m Vector [log m_0, log m_1, log m_2, log m_3]
  *
  * @return Vector [log G_1(t), log G_2(t), log G_3(t)]
  */
vector primarycensored_moments_from_partial(real log_t, vector log_m) {
  real log_g1 = primarycensored_log_diff_exp(log_t + log_m[1], log_m[2]);
  real log_h = primarycensored_log_diff_exp(log_t + log_m[2], log_m[3]);
  real log_g2 = primarycensored_log_diff_exp(log_t + log_g1, log_h);
  real log_a = primarycensored_log_diff_exp(log_t + log_m[3], log_m[4]);
  real log_b = primarycensored_log_diff_exp(log_t + log_h, log_a);
  real log_g3 = primarycensored_log_diff_exp(log_t + log_g2, log_b);
  return [log_g1, log_g2, log_g3]';
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
  real log_t = log(t);
  real log_x = log_t + log(rate);
  // log_mk is the log of int_0^t u^k f(u) du
  vector[2] log_p01 = gamma_lcdf_logx_pair(log_x, shape);
  vector[2] log_p23 = gamma_lcdf_logx_pair(log_x, shape + 2);
  real log_m0 = log_p01[1];
  real log_m1 = log(shape) - log(rate) + log_p01[2];
  real log_m2 = log(shape) + log(shape + 1) - 2 * log(rate) + log_p23[1];
  real log_m3 = log(shape) + log(shape + 1) + log(shape + 2) - 3 * log(rate)
                + log_p23[2];
  return primarycensored_moments_from_partial(
    log_t, [log_m0, log_m1, log_m2, log_m3]'
  );
}

/**
  * Log moments of a Stacy family delay about a point
  * @ingroup tilt_transforms
  *
  * As primarycensored_gamma_tilt_moments(), from the partial moments
  * int_0^t u^j f(u) du = scale^j Gamma(k + j / shape) / Gamma(k)
  * P(k + j / shape, (t / scale)^shape).
  *
  * @param t Point, positive
  * @param shape Power of the Stacy family
  * @param scale Scale
  * @param k Shape of the underlying gamma
  *
  * @return Vector [log G_1(t), log G_2(t), log G_3(t)]
  */
vector primarycensored_stacy_tilt_moments(real t, real shape, real scale,
                                          real k) {
  real log_x = shape * (log(t) - log(scale));
  real k1 = k + 1 / shape;
  real k2 = k + 2 / shape;
  real k3 = k + 3 / shape;
  real log_k = lgamma(k);
  real log_t = log(t);
  real log_m0 = gamma_lcdf_logx(log_x, k);
  real log_m1 = log(scale) + lgamma(k1) - log_k
                + gamma_lcdf_logx(log_x, k1);
  real log_m2 = 2 * log(scale) + lgamma(k2) - log_k
                + gamma_lcdf_logx(log_x, k2);
  real log_m3 = 3 * log(scale) + lgamma(k3) - log_k
                + gamma_lcdf_logx(log_x, k3);
  return primarycensored_moments_from_partial(
    log_t, [log_m0, log_m1, log_m2, log_m3]'
  );
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
  * G_3 = sigma^3 (z (z^2 + 3) Phi(z) + (z^2 + 2) phi(z)). The lognormal
  * uses its partial moments
  * m_k(t) = exp(k mu + k^2 sigma^2 / 2) Phi(z - k sigma).
  *
  * @param t Point
  * @param dist_id Distribution identifier: 1 (Lognormal), 2 (Gamma),
  *   3 (Weibull), 4 (Exponential), 5 (Generalised gamma) or 18 (Normal), see
  *   check_for_tilt_transform()
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
  } else if (dist_id == 3 || dist_id == 5) {
    if (t <= 0) return rep_vector(negative_infinity(), 3);
    return primarycensored_stacy_tilt_moments(
      t, params[1], params[2], dist_id == 3 ? 1.0 : params[3]
    );
  } else if (dist_id == 1) {
    if (t <= 0) return rep_vector(negative_infinity(), 3);
    real mu = params[1];
    real sigma = params[2];
    real log_t = log(t);
    real z = (log_t - mu) / sigma;
    vector[4] log_m;
    for (k in 0:3) {
      log_m[k + 1] = k * mu + 0.5 * square(k * sigma)
                     + primarycensored_log_std_normal_cdf(z - k * sigma);
    }
    return primarycensored_moments_from_partial(log_t, log_m);
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
 * The direct form cancels as rho goes to zero, so small tilt forms replace
 * it, see primarycensored_exptilt_form().
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
  return primary_id == 2
         && (dist_id == 1 || dist_id == 2 || dist_id == 3 || dist_id == 4
             || dist_id == 5 || dist_id == 18);
}

/**
  * Compute the exponentially tilted terms at an endpoint
  * @ingroup exponential_tilt_solutions
  *
  * @param t Endpoint, d or d - pwindow
  * @param dist_id Distribution identifier, see check_for_exptilt()
  * @param rho Tilt
  * @param params Array of distribution parameters
  * @param context Output of log_tilt_transform_context() for -rho
  *
  * @return Vector [log F(t), log(1 - F(t)), log J(t), log(J(Inf) - J(t)),
  * log loss of J], see log_tilt_transform_pair_shared(). The upper tail term
  * of J is `inf` where the total diverges, as for a lognormal delay with
  * rho < 0.
  */
vector primarycensored_exptilt_terms(real t, int dist_id, real rho,
                                     array[] real params, vector context) {
  return append_row(
    log_tilt_transform_pair(t, dist_id, 0, params)[1:2],
    log_tilt_transform_pair_shared(t, dist_id, -rho, params, context)
  );
}

/**
  * Choose the form of the exponentially tilted log CDF at a delay
  * @ingroup exponential_tilt_solutions
  *
  * The small window form is used for |rho| * pwindow below 1e-2, or 1e-3 for
  * delays on the reals, and |rho| * (|d| + pwindow) below 0.1, as it
  * subtracts terms of the size of |d| + pwindow. Otherwise, for delays on
  * the non-negative reals with d < pwindow and |rho| * d below 5e-3, the
  * small delay form is used. Both truncate at third order in the tilt, with
  * a relative error below (|rho| pwindow)^3 / 160 and (|rho| d)^3 / 24, so
  * below 1e-8 at these limits. The direct form is used elsewhere.
  *
  * @param d Delay
  * @param dist_id Distribution identifier for the delay distribution
  * @param rho Tilt
  * @param pwindow Primary event window
  *
  * @return 1 for the small window form, 2 for the small delay form and 0 for
  * the direct form
  */
int primarycensored_exptilt_form(data real d, int dist_id, real rho,
                                 data real pwindow) {
  int positive = dist_has_positive_support(dist_id);
  if (abs(rho) * pwindow < (positive ? 1e-2 : 1e-3)
      && abs(rho) * (abs(d) + pwindow) < 0.1) {
    return 1;
  }
  if (positive && d < pwindow && abs(rho) * d < 5e-3) return 2;
  return 0;
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
  * @param upper_d Log upper tail quantity at d, `nan` if there is none
  * @param upper_q Log upper tail quantity at q, `nan` if there is none
  *
  * @return Log of the difference, `-inf` if it is zero to rounding. The
  * lower tail terms are used where the upper tail terms are `inf`.
  */
real primarycensored_tail_diff(real lower_d, real lower_q, real upper_d,
                               real upper_q) {
  if (is_nan(upper_d) || is_nan(upper_q)) {
    return primarycensored_log_diff_exp(lower_d, lower_q);
  }
  // NaN from terms that both underflow is a zero difference
  if (is_nan(lower_q - lower_d) || is_inf(upper_q)
      || lower_q - lower_d <= upper_d - upper_q) {
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
  * The result is `nan` if J at d or q is. Where J has no upper tail, it is
  * also `nan` if the terms of the numerator over the denominator exceed 1e5.
  *
  * @param terms_d Terms at d, from primarycensored_exptilt_terms()
  * @param terms_q Terms at q = d - pwindow
  * @param d Delay
  * @param rho Tilt, not zero
  * @param pwindow Primary event window
  *
  * @return Log of the primary event censored CDF at d, or `nan`
  */
real primarycensored_exptilt_lcdf_from_terms(vector terms_d, vector terms_q,
                                             data real d, real rho,
                                             data real pwindow) {
  if (is_nan(terms_d[3]) || is_nan(terms_q[3])) return not_a_number();
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
  real result = fmin(log_sum_exp(terms_q[1], log_num - log_den), 0);
  int no_tail = is_nan(terms_d[4]) || is_nan(terms_q[4]);
  real log_size = rho * d + terms_d[3] + fmax(terms_d[5], terms_q[5]);
  if (no_tail && result > negative_infinity()
      && log_size - log_den > log(1e5)) {
    return not_a_number();
  }
  return result;
}

/**
  * Combine the moments at d and q into the small tilt log CDF
  * @ingroup exponential_tilt_solutions
  *
  * The uniform window limit with its corrections to second order in the tilt,
  * for the small window form of primarycensored_exptilt_form(). With G_k(t) = int (t - u)^k f(u) du,
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
  * For the small delay form of primarycensored_exptilt_form() the terms at
  * d - pwindow are zero and
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
  * gives the delay log CDF. Where the direct form is `nan` it uses the ODE
  * path of primarycensored_numeric_cdf(). Only for check_for_exptilt() is 1
  * and check_for_tilt_transform() is 1 for -rho.
  *
  * @param d Delay
  * @param dist_id Distribution identifier: 1 (Lognormal), 2 (Gamma),
  *   3 (Weibull), 4 (Exponential), 5 (Generalised gamma) or 18 (Normal)
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
  int form = primarycensored_exptilt_form(d, dist_id, rho, pwindow);
  if (form == 1) {
    return primarycensored_exptilt_small_window_lcdf_from_terms(
      primarycensored_tilt_moments(d, dist_id, params),
      primarycensored_tilt_moments(q, dist_id, params), rho, pwindow
    );
  }
  if (form == 2) {
    return primarycensored_exptilt_small_delay_lcdf_from_terms(
      primarycensored_tilt_moments(d, dist_id, params), rho, pwindow
    );
  }
  vector[4] context = log_tilt_transform_context(dist_id, -rho, params);
  real log_cdf = primarycensored_exptilt_lcdf_from_terms(
    primarycensored_exptilt_terms(d, dist_id, rho, params, context),
    primarycensored_exptilt_terms(q, dist_id, rho, params, context),
    d, rho, pwindow
  );
  if (is_nan(log_cdf)) {
    return log(
      primarycensored_numeric_cdf(d | dist_id, params, pwindow, 2, {rho})
    );
  }
  return log_cdf;
}
