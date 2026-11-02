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
    dist_id == 18 || dist_id == 19 || dist_id == 21 || dist_id == 22 ||
    dist_id == 29;
}

/**
  * Check the uniform primary terms are valid for the delay parameters
  * @ingroup analytical_solution_helpers
  *
  * The Inverse gamma terms need a finite mean, `alpha > 1`, and the Inverse
  * chi-square (19) and Scaled inverse chi-square (22) need `nu > 2`. The
  * Log-logistic (29) partial expectation needs a shape of at least 0.01.
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
  if (dist_id == 29) return params[2] >= 0.01;
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
  * Check if the analytical solution applies for the given parameters
  * @ingroup analytical_solution_helpers
  *
  * This is check_for_analytical() with check_uniform_terms_params() and, for
  * the exponentially tilted solutions, which check_for_analytical() does not
  * include, that the tilted delay exists, see check_for_tilt_transform(). It
  * chooses the path in primarycensored_cdf() and primarycensored_lcdf().
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
  return check_for_analytical(dist_id, primary_id) &&
    check_uniform_terms_params(dist_id, params);
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
  * F_T(t; k) comes from F_T(t; k + 1) and the recursion
  * P(k, y) = P(k + 1, y) + y^k e^{-y} / Gamma(k + 1), with y = rate * t,
  * which adds terms so stays accurate in the lower tail. Rejects a rate that
  * is not positive and finite.
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
  real log_t = log(t);
  real log_y = log_t + log_rate;
  real log_F_T_kp1 = gamma_lcdf_logx(log_y, shape + 1);
  real log_lead = shape * log_y - exp(log_y) - lgamma(shape + 1);
  return [
    log_t + log_sum_exp(log_F_T_kp1, log_lead),
    log_E + log_F_T_kp1
  ]';
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
  * it does not underflow in the lower tail.
  *
  * @param y Positive value at which the log CDF is evaluated
  * @param alpha Shape parameter
  * @param beta Scale parameter
  *
  * @return Log CDF at y. Rejects if the series or continued fraction does
  * not converge.
  */
real primarycensored_inv_gamma_lcdf(real y, real alpha, real beta) {
  real x = beta / y;
  // Both need about 9 sqrt(alpha) terms at worst, near x = alpha + 1
  real max_terms = 20 * sqrt(alpha) + 200;
  if (x > alpha + 1) {
    // Within 1e-9 of an integer the fraction terminates, or nearly so, at
    // i = round(alpha). Its value has then converged but its shape gradient
    // has not, so the loop runs 50 terms past alpha before it stops.
    int near_int = abs(alpha - round(alpha)) < 1e-9;
    if (near_int) max_terms += alpha + 50;
    real b = x + 1 - alpha;
    // Modified Lentz starts c at infinity and floors |c| and |d| at 1e-300
    // to avoid dividing by zero
    real c = 1e300;
    real d = inv(b);
    real h = d;
    int converged = 0;
    int i = 0;
    while (!converged) {
      i += 1;
      if (i > max_terms) {
        reject("primarycensored_inv_gamma_lcdf: continued fraction did not ",
               "converge for x = ", x, " and alpha = ", alpha);
      }
      real an = -i * (i - alpha);
      b += 2;
      d = an * d + b;
      d = inv(abs(d) < 1e-300 ? 1e-300 : d);
      c = b + an / c;
      c = abs(c) < 1e-300 ? 1e-300 : c;
      h *= d * c;
      // An update within 1e-15 of 1 leaves h unchanged in double precision
      converged = abs(d * c - 1) < 1e-15 && (!near_int || i > alpha + 50);
    }
    return (alpha - 1) * log(x) - x - lgamma(alpha) + log(x) + log(h);
  }
  real term = 1;
  real series = 1;
  int n = 0;
  // A term below 1e-17 of the sum leaves it unchanged in double precision
  while (term >= 1e-17 * series) {
    n += 1;
    if (n > max_terms) {
      reject("primarycensored_inv_gamma_lcdf: series did not converge for ",
             "x = ", x, " and alpha = ", alpha);
    }
    term *= x / (alpha + n);
    series += term;
  }
  return log1m_exp(alpha * log(x) - x - lgamma(alpha + 1) + log(series));
}

/**
  * Compute the uniform primary terms at t for an Inverse gamma delay
  * @ingroup analytical_solution_helpers
  *
  * The mean is E = beta / (alpha - 1) and the partial expectation
  * distribution is the Inverse gamma with shape `alpha - 1`, so `alpha > 1`
  * is required. This is not checked here, so callers must first apply
  * check_uniform_terms_params().
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
  *   21: Pareto, 22: Scaled inverse chi-square, 29: Log-logistic), see
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
  // Delays on the reals are not clipped at 0
  real q = dist_has_positive_support(dist_id) ? max({d - pwindow, 0.0})
                                              : d - pwindow;

  if (dist_id == 29 && primary_id == 1) {
    return primarycensored_loglogistic_uniform_lcdf(d | q, params, pwindow);
  }
  if (check_for_uniform_terms(dist_id, primary_id)) {
    // The Beta CDF is exactly 1 once the window is above its support
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
  * The analytical CDF at d combines terms at d and at q = max(d - pwindow, 0)
  * (d - pwindow for delays on the reals). With an integer pwindow q is an
  * integer delay too, so primarycensored_analytical_lcdf_vectorized() can
  * compute the terms once per delay and share them. This needs the
  * analytical solutions built from primarycensored_uniform_terms(), see
  * check_for_uniform_terms(), or the exponentially tilted solutions, see
  * check_for_exptilt(). The
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
  * The log CDF at d combines the terms at d and at q (see
  * primarycensored_uniform_lcdf_from_terms()). Both are integer delays, so
  * the terms are computed once per delay and used for both, halving the CDF
  * evaluations. The values are the same as from
  * primarycensored_analytical_lcdf() at each delay without truncation.
  * For an exponentially tilted primary see
  * primarycensored_exptilt_lcdf_vectorized().
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
    return primarycensored_exptilt_lcdf_vectorized(
      start, n, dist_id, params, pw, primary_params[1]
    );
  }
  int first = dist_has_positive_support(dist_id) ? max(start - pw, 0)
                                                 : start - pw;
  // terms[t - first + 1] holds the terms at delay t
  array[n - first + 1] vector[2] terms;
  for (t in first:n) {
    terms[t - first + 1] = primarycensored_uniform_terms(t, dist_id, params);
  }
  // The numerical CDF can be out of order by its tolerance, which would make
  // the PMF NaN
  int previous_numerical = 0;
  for (d in start:n) {
    int q = dist_has_positive_support(dist_id) ? max(d - pw, 0) : d - pw;
    vector[2] terms_q = terms[q - first + 1];
    log_cdfs[d] = (dist_id == 9 && q >= 1) ? 0
      : primarycensored_uniform_lcdf_from_terms(
          terms[d - first + 1], terms_q, pwindow
        );
    if (dist_id == 29) {
      int numerical = loglogistic_window_ill_conditioned(
        terms[d - first + 1][1], terms_q[1], d, pwindow, log_cdfs[d]
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
