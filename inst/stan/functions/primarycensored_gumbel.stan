/*
 * Truncated Gumbel primary event window
 *
 * For a primary window of width w with the truncated Gumbel density of
 * tgumbel.stan on [0, w], G(z) = exp(-exp(-(z - mu) / beta)) and
 * D = G(w) - G(0), the primary event censored CDF at delay d is
 *   F_G(d) = F(q) + B(d) / D, with q = d - w and
 *   B(d) = (1 - G(0)) (F(d) - F(q))
 *          + sum_{n >= 1} (-c)^n / n! (T_f(n / beta; d) - T_f(n / beta; q)),
 * where c = exp(-(d - mu) / beta), F is the delay CDF and T_f(xi; t) is the
 * tilt transform of the delay, see tilt_transform.stan. F(x) = T_f(xi; x) = 0
 * for x <= 0 for delays on the non-negative reals. The terms depend on d or
 * q alone through T_f, so the transforms at every tilt n / beta are computed
 * once per endpoint by primarycensored_gumbel_terms() and shared between
 * neighbouring delays in primarycensored_gumbel_lcdf_vectorized().
 *
 * The series alternates. Its terms are at most s0^n / n! times the delay
 * mass in the window, where s0 = exp(mu / beta), so the series is used only
 * for mu / beta below log(15), with n_terms from gumbel_n_terms(). The terms
 * are evaluated on the log scale and the odd and even terms are accumulated
 * separately. The relative error is estimated for each delay, and the
 * numerical path, see primarycensored_gumbel_numeric_cdf(), is used where it
 * is above 1e-8, see gumbel_error_tolerance(). The
 * transform is needed at the tilts n / beta, which are positive, so the
 * exponential and gamma delays need a rate above n_terms / beta, see
 * check_for_gumbel_params().
 *
 * A delay distribution plugs in through tilt_transform.stan and is then
 * supported by the window.
 */

/**
  * Largest value of mu / beta for which the Gumbel series is used
  * @ingroup truncated_gumbel_solutions
  *
  * The terms of the series are as large as exp(mu / beta), and the number
  * of terms grows with it, so the series is not used beyond log(15).
  *
  * @return log(15)
  */
real gumbel_max_log_s0() {
  return log(15);
}

/**
  * Largest estimated relative error for which the Gumbel series is used
  * @ingroup truncated_gumbel_solutions
  *
  * The estimate is 10 to 1000 times the actual error, which was at most
  * 5e-10 relative for an estimate of at most 1e-8 in tests against a tight
  * reference. The numerical path, see primarycensored_gumbel_numeric_cdf(),
  * is accurate to about 1e-11 but 3 to 80 times slower, and loses relative
  * precision for a CDF below 1e-10, where the series is relative. So the
  * series is kept up to the estimate where its error is of the order of
  * 1e-9, and replaced by the numerical path above it. The two are within
  * 1e-9 of each other at the switch.
  *
  * @return 1e-8
  */
real gumbel_error_tolerance() {
  return 1e-8;
}

/**
  * Number of terms of the Gumbel series
  * @ingroup truncated_gumbel_solutions
  *
  * The smallest N, at least 2, such that the bound s0^n / n! on the terms
  * after N is below 1e-18 max(s0, 1), where s0 = exp(mu / beta).
  *
  * @param log_s0 mu / beta
  *
  * @return The number of terms
  */
int gumbel_n_terms(real log_s0) {
  int n = 2;
  while (n < 200
         && n * log_s0 - lgamma(n + 1) - fmax(log_s0, 0) >= log(1e-18)) {
    n += 1;
  }
  return n;
}

/**
  * Check if the analytical solution is the truncated Gumbel solution
  * @ingroup truncated_gumbel_solutions
  *
  * The delay needs a tilt transform, see check_for_tilt_transform(), and the
  * truncated Gumbel primary is primary_id 4 with primary_params = [mu, beta].
  *
  * @param dist_id Distribution identifier for the delay distribution
  * @param primary_id Distribution identifier for the primary distribution
  *
  * @return 1 if the delay has a truncated Gumbel solution and the primary is
  * the truncated Gumbel, 0 otherwise. Whether it applies for given
  * parameters is check_for_analytical_params().
  */
int check_for_gumbel(int dist_id, int primary_id) {
  return primary_id == 4 && (dist_id == 2 || dist_id == 4 || dist_id == 18);
}

/**
  * Check if the Gumbel series applies for the given parameters
  * @ingroup truncated_gumbel_solutions
  *
  * The series needs mu / beta below gumbel_max_log_s0() and the tilt
  * transform of the delay at the largest tilt n_terms / beta, see
  * check_for_tilt_transform(). Whether it is accurate at a delay is decided
  * by the error estimate of primarycensored_gumbel_lcdf_from_terms().
  *
  * @param dist_id Distribution identifier for the delay distribution
  * @param params Array of delay distribution parameters
  * @param mu Location of the truncated Gumbel
  * @param beta Scale of the truncated Gumbel
  *
  * @return 1 if the series can be used, 0 otherwise
  */
int check_for_gumbel_params(int dist_id, array[] real params, real mu,
                            real beta) {
  if (!(beta > 0) || is_inf(beta) || is_nan(mu) || is_inf(mu)) return 0;
  if (mu / beta > gumbel_max_log_s0()) return 0;
  return check_for_tilt_transform(
    dist_id, gumbel_n_terms(mu / beta) / beta, params
  );
}

/**
  * Log of a sum of two terms, either of which may be `-inf`
  * @ingroup truncated_gumbel_solutions
  *
  * log_sum_exp() of two `-inf` differentiates to NaN, so a `-inf` term is
  * dropped.
  *
  * @param a Log of the first term
  * @param b Log of the second term
  *
  * @return log(exp(a) + exp(b))
  */
real gumbel_log_sum_exp(real a, real b) {
  if (a == negative_infinity()) return b;
  if (b == negative_infinity()) return a;
  return log_sum_exp(a, b);
}

/**
  * Absolute value, zero if not finite
  * @ingroup truncated_gumbel_solutions
  *
  * @param x Number
  *
  * @return abs(x), or 0 if x is infinite
  */
real gumbel_finite_abs(real x) {
  return is_inf(x) || is_nan(x) ? 0 : abs(x);
}

/**
  * Compute the truncated Gumbel terms at an endpoint
  * @ingroup truncated_gumbel_solutions
  *
  * @param t Endpoint, d or d - pwindow
  * @param dist_id Distribution identifier, see check_for_gumbel()
  * @param beta Scale of the truncated Gumbel
  * @param n_terms Number of terms, see gumbel_n_terms()
  * @param params Array of distribution parameters
  *
  * @return Vector of length 2 (n_terms + 1) with, for each tilt
  * xi = n / beta, n = 0, ..., n_terms, the pair
  * [log T_f(xi; t), log(T_f(xi; Inf) - T_f(xi; t))] of
  * log_tilt_transform_pair() at positions 2 n + 1 and 2 n + 2. The lower terms
  * are `-inf` for t <= 0 for delays on the non-negative reals. Only defined
  * where check_for_gumbel_params() is 1.
  */
vector primarycensored_gumbel_terms(real t, int dist_id, real beta,
                                    int n_terms, array[] real params) {
  vector[2 * (n_terms + 1)] terms;
  real inv_beta = inv(beta);
  for (n in 0:n_terms) {
    terms[(2 * n + 1):(2 * n + 2)] = log_tilt_transform_pair(
      t, dist_id, n * inv_beta, params
    );
  }
  return terms;
}

/**
  * Combine the truncated Gumbel terms at d and q into the log CDF
  * @ingroup truncated_gumbel_solutions
  *
  * The bracket B is the difference of the positive and negative parts of the
  * series, (1 - G(0)) (F(d) - F(q)) plus the even terms, and the odd terms.
  * Each difference T_f(xi; d) - T_f(xi; q) is taken between the lower tail
  * terms or between the upper tail terms, whichever loses less precision,
  * see primarycensored_tail_diff().
  *
  * The estimated relative error is the machine precision, times the sum of
  * the absolute terms over B, times one plus the largest magnitude of a log
  * term, plus the bound s0^(N + 1) / (N + 1)! on the truncation error times
  * the delay mass in the window over B. It is `inf` where B is not
  * positive.
  *
  * @param terms_d Terms at d from primarycensored_gumbel_terms()
  * @param terms_q Terms at q = d - pwindow
  * @param d Delay
  * @param pwindow Primary event window
  * @param mu Location of the truncated Gumbel
  * @param beta Scale of the truncated Gumbel
  * @param n_terms Number of terms, see gumbel_n_terms()
  *
  * @return Vector [log of the primary event censored CDF at d, estimated
  * relative error]
  */
vector primarycensored_gumbel_lcdf_from_terms(vector terms_d, vector terms_q,
                                              data real d,
                                              data real pwindow, real mu,
                                              real beta, int n_terms) {
  real log_s0 = mu / beta;
  real log_c = -(d - mu) / beta;
  real log_f_q = terms_q[1];
  real log_delta0 = primarycensored_tail_diff(
    terms_d[1], terms_q[1], terms_d[2], terms_q[2]
  );
  // No delay mass in the window, so every term vanishes together
  if (log_delta0 == negative_infinity()) {
    return [log_f_q, 0]';
  }
  real log_pos = log1m_exp(-exp(log_s0)) + log_delta0;
  real log_neg = negative_infinity();
  real scale = 0;
  for (n in 1:n_terms) {
    int i = 2 * n + 1;
    real log_delta = primarycensored_tail_diff(
      terms_d[i], terms_q[i], terms_d[i + 1], terms_q[i + 1]
    );
    real log_a = n * log_c + log_delta - lgamma(n + 1);
    if (n % 2 == 1) {
      log_neg = gumbel_log_sum_exp(log_neg, log_a);
    } else {
      log_pos = gumbel_log_sum_exp(log_pos, log_a);
    }
    scale = fmax(
      scale,
      n * abs(log_c) + fmax(
        fmax(gumbel_finite_abs(terms_d[i]), gumbel_finite_abs(terms_d[i + 1])),
        fmax(gumbel_finite_abs(terms_q[i]), gumbel_finite_abs(terms_q[i + 1]))
      )
    );
  }
  real log_bracket = primarycensored_log_diff_exp(log_pos, log_neg);
  // The normalisation G(w) - G(0) is exp(-s_w) (1 - exp(-delta)), with
  // delta = s0 - s_w = s_w (exp(w / beta) - 1) at most exp(mu / beta)
  real log_s_w = -(pwindow - mu) / beta;
  real log_d = -exp(log_s_w)
               + log1m_exp(-exp(log_s_w + log_diff_exp(pwindow / beta, 0)));
  real error;
  if (log_bracket == negative_infinity()) {
    error = positive_infinity();
  } else {
    int n_last = n_terms + 1;
    real log_rounding = log(machine_precision()) + log1p(scale)
                        + gumbel_log_sum_exp(log_pos, log_neg) - log_bracket;
    real log_trunc = log_delta0 + n_last * log_s0 - lgamma(n_last + 1)
                     - log1m(fmin(exp(log_s0) / (n_last + 1), 0.5))
                     - log_bracket;
    error = exp(log_rounding) + exp(log_trunc);
  }
  real log_cdf = gumbel_log_sum_exp(log_f_q, log_bracket - log_d);
  return [fmin(log_cdf, 0), error]';
}

/**
  * Largest value of u = s(z) - s(pwindow) that is integrated
  * @ingroup truncated_gumbel_solutions
  *
  * The quantile variable of the window is at most 1 - 2^-53 in double
  * precision, which is u of about 36.7, so the mass beyond 36 (2.3e-16 of
  * the window) is not resolved and is given the value at 36.
  *
  * @return 36
  */
real gumbel_numeric_u_max() {
  return 36;
}

/**
  * Limits of the numerical integral of a truncated Gumbel primary
  * @ingroup truncated_gumbel_solutions
  *
  * The window density is a narrow spike when mu / beta is large, which an
  * adaptive integrator can step over. Two variables are used, both with a
  * smooth integrand. With mu at or above the window end, s(pwindow) =
  * exp((mu - pwindow) / beta) is at least 1 and the window is integrated in
  * its upper quantile v, from the end of the window, where the density is
  * exp(-u) / (1 - exp(-Delta)) for u = s(z) - s(pwindow) and Delta =
  * s(0) - s(pwindow). The limits are the quantile of the point d for delays
  * on the non-negative reals, below which the delay CDF is 0, and 1. With
  * mu below the end, the window is integrated in z over the range that holds
  * its mass, from where s(z) - s(pwindow) is 745 or more, to 46 beta above
  * the larger of mu and 0, cut at d for delays on the non-negative reals.
  * The mass outside is below 1e-16, see gumbel_numeric_u_max(), or 1e-20.
  *
  * @param d Delay
  * @param dist_id Distribution identifier
  * @param pwindow Primary event window
  * @param mu Location of the truncated Gumbel
  * @param beta Scale of the truncated Gumbel
  *
  * @return Vector [lower, upper] of limits of the integration variable. The
  * integral is 0 if upper is not above lower.
  */
vector gumbel_numeric_limits(data real d, int dist_id, data real pwindow,
                             real mu, real beta) {
  int positive = dist_has_positive_support(dist_id);
  real log_sw = (mu - pwindow) / beta;
  if (log_sw >= 0) {
    if (!positive || d >= pwindow) {
      return [0, 1]';
    }
    real log_delta = log_sw + log_diff_exp(pwindow / beta, 0);
    real u_cap = exp(fmin(log_delta, log(gumbel_numeric_u_max())));
    real log_uk = log_sw + log_diff_exp((pwindow - d) / beta, 0);
    if (log_uk >= log(u_cap)) {
      return [1, 1]';
    }
    real c = log_delta > 4 ? 1 : -expm1(-exp(log_delta));
    return [-expm1(-exp(log_uk)) / c, 1]';
  }
  real z_low = fmax(0, mu - beta * log(745 + exp(log_sw)));
  real z_high = fmin(pwindow, fmax(mu, 0) + 46 * beta);
  if (positive) {
    z_high = fmin(z_high, d);
  }
  return [z_low, z_high]';
}

/**
  * ODE system for the truncated Gumbel primary event censored CDF
  * @ingroup truncated_gumbel_solutions
  *
  * The integration variable tau is in [0, 1] and is mapped to the limits of
  * gumbel_numeric_limits(). With mu at or above the window end it is the
  * upper quantile v, the density of v is 1, and z = pwindow - beta log(1 +
  * u / s(pwindow)). The delay is taken as d - z = (d - pwindow) + beta log(1
  * + u / s(pwindow)), which keeps its small difference from 0 where the
  * spike is within rounding of pwindow. With mu below it the integrand is
  * the delay CDF times the window density in z.
  *
  * @param tau Integration variable in [0, 1]
  * @param y State, the integral so far
  * @param d Delay
  * @param pwindow Primary event window
  * @param dist_id Distribution identifier
  * @param params Array of distribution parameters
  * @param mu Location of the truncated Gumbel
  * @param beta Scale of the truncated Gumbel
  *
  * @return The derivative of the state with respect to tau
  */
vector primarycensored_gumbel_ode(real tau, vector y, data real d,
                                  data real pwindow, int dist_id,
                                  array[] real params, real mu, real beta) {
  vector[2] limits = gumbel_numeric_limits(d, dist_id, pwindow, mu, beta);
  real width = limits[2] - limits[1];
  real x = limits[1] + width * tau;
  real log_sw = (mu - pwindow) / beta;
  if (log_sw >= 0) {
    real log_delta = log_sw + log_diff_exp(pwindow / beta, 0);
    real u_cap = exp(fmin(log_delta, log(gumbel_numeric_u_max())));
    real c = log_delta > 4 ? 1 : -expm1(-exp(log_delta));
    // x is the upper quantile v, so 1 - v c is exp(-u)
    real u = -log1m(fmin(x * c, -expm1(-u_cap)));
    real delay = (d - pwindow) + beta * log1p(u * exp(-log_sw));
    return rep_vector(width * exp(dist_lcdf(delay | params, dist_id)), 1);
  }
  return rep_vector(
    width * exp(dist_lcdf(d - x | params, dist_id)
                + tgumbel_lpdf(x | 0, pwindow, mu, beta)),
    1
  );
}

/**
  * Primary event censored CDF for a truncated Gumbel primary by numerical
  * integration
  * @ingroup truncated_gumbel_solutions
  *
  * Integrates the delay CDF against the window density with an ODE solver
  * in a variable in which the integrand is smooth for every mu and beta, see
  * gumbel_numeric_limits(). It is primarycensored_numeric_cdf() for the
  * truncated Gumbel primary, whose integration over the window would miss a
  * spike of the density. The solver tolerances are a relative 1e-12 and an
  * absolute 1e-18, and the result was within 1e-11 of a tight reference in
  * tests, relative to a CDF above 1e-6 and absolute below it. The result is
  * in [0, 1].
  *
  * @param d Delay
  * @param dist_id Distribution identifier
  * @param params Array of distribution parameters
  * @param pwindow Primary event window
  * @param mu Location of the truncated Gumbel
  * @param beta Scale of the truncated Gumbel
  *
  * @return Primary event censored CDF, not normalized for truncation
  */
real primarycensored_gumbel_numeric_cdf(data real d, int dist_id,
                                        array[] real params,
                                        data real pwindow, real mu,
                                        real beta) {
  if (dist_has_positive_support(dist_id) && d <= 0) {
    return 0;
  }
  vector[2] limits = gumbel_numeric_limits(d, dist_id, pwindow, mu, beta);
  if (!(limits[2] > limits[1])) {
    return 0;
  }
  vector[1] y0 = rep_vector(0.0, 1);
  real result = ode_rk45_tol(
    primarycensored_gumbel_ode, y0, 0.0, {1.0}, 1e-12, 1e-18, 1000000,
    d, pwindow, dist_id, params, mu, beta
  )[1, 1];
  return fmin(1, fmax(0, result));
}

/**
  * Log of the numerical primary event censored CDF for a truncated Gumbel
  * primary
  * @ingroup truncated_gumbel_solutions
  *
  * A CDF that is not above 0, which is below the solver tolerance, is `-inf`
  * and one above 1 is capped at 0, so the result is never NaN and never
  * above 0.
  *
  * @param d Delay
  * @param dist_id Distribution identifier
  * @param params Array of distribution parameters
  * @param pwindow Primary event window
  * @param mu Location of the truncated Gumbel
  * @param beta Scale of the truncated Gumbel
  *
  * @return Log of primarycensored_gumbel_numeric_cdf()
  */
real primarycensored_gumbel_numeric_lcdf(data real d, int dist_id,
                                         array[] real params,
                                         data real pwindow, real mu,
                                         real beta) {
  real cdf = primarycensored_gumbel_numeric_cdf(
    d | dist_id, params, pwindow, mu, beta
  );
  return cdf > 0 ? fmin(log(cdf), 0) : negative_infinity();
}

/**
  * Compute the primary event censored log CDF for a truncated Gumbel
  * primary
  * @ingroup truncated_gumbel_solutions
  *
  * Uses the series where its estimated relative error is at most
  * gumbel_error_tolerance(), and the ODE path, see
  * primarycensored_numeric_cdf(), for the delays where it is not. Only for
  * check_for_gumbel() is 1 and check_for_gumbel_params() is 1.
  *
  * @param d Delay
  * @param dist_id Distribution identifier: 2 (Gamma), 4 (Exponential) or 18
  *   (Normal)
  * @param params Array of distribution parameters
  * @param pwindow Primary event window
  * @param mu Location of the truncated Gumbel
  * @param beta Scale of the truncated Gumbel
  *
  * @return Log of the primary event censored CDF at d
  */
real primarycensored_gumbel_lcdf(data real d, int dist_id,
                                 array[] real params, data real pwindow,
                                 real mu, real beta) {
  if (dist_has_positive_support(dist_id) && d <= 0) {
    return negative_infinity();
  }
  int n_terms = gumbel_n_terms(mu / beta);
  vector[2] fit = primarycensored_gumbel_lcdf_from_terms(
    primarycensored_gumbel_terms(d, dist_id, beta, n_terms, params),
    primarycensored_gumbel_terms(d - pwindow, dist_id, beta, n_terms, params),
    d, pwindow, mu, beta, n_terms
  );
  if (fit[2] <= gumbel_error_tolerance()) {
    return fit[1];
  }
  return primarycensored_gumbel_numeric_lcdf(
    d | dist_id, params, pwindow, mu, beta
  );
}

/**
  * Check if the truncated Gumbel solution can be vectorised over integer
  * delays
  * @ingroup truncated_gumbel_solutions
  *
  * With an integer pwindow, q = d - pwindow is an integer delay too, so
  * primarycensored_gumbel_lcdf_vectorized() can compute the terms once per
  * delay and share them.
  *
  * @param dist_id Distribution identifier for the delay distribution
  * @param primary_id Distribution identifier for the primary distribution
  * @param pwindow Primary event window
  *
  * @return 1 if the vectorised truncated Gumbel solution applies, 0
  * otherwise
  */
int check_for_gumbel_vectorized(int dist_id, int primary_id,
                                data real pwindow) {
  return check_for_gumbel(dist_id, primary_id)
         && pwindow >= 1 && floor(pwindow) == pwindow;
}

/**
  * Compute the truncated Gumbel primary event censored log CDF at integer
  * delays
  * @ingroup truncated_gumbel_solutions
  *
  * The log CDF at d combines the terms at d and at q = d - pwindow. Both are
  * integer delays, so the terms are computed once per delay and used for
  * both, halving the transform evaluations. The values are the same as from
  * primarycensored_gumbel_lcdf() at each delay. Only for cases where
  * check_for_gumbel_vectorized() is 1 and check_for_gumbel_params() is 1.
  *
  * @param start First delay to compute
  * @param n Last delay to compute, and the length of the result
  * @param dist_id Distribution identifier: 2 (Gamma), 4 (Exponential) or 18
  *   (Normal)
  * @param params Array of distribution parameters
  * @param pwindow Primary event window, a positive integer
  * @param mu Location of the truncated Gumbel
  * @param beta Scale of the truncated Gumbel
  *
  * @return Vector whose element d is the log CDF at d, for d in start:n.
  * Elements before start are not computed.
  */
vector primarycensored_gumbel_lcdf_vectorized(data int start, data int n,
                                              data int dist_id,
                                              array[] real params,
                                              data real pwindow, real mu,
                                              real beta) {
  int pw = to_int(pwindow);
  int positive = dist_has_positive_support(dist_id);
  int n_terms = gumbel_n_terms(mu / beta);
  // Endpoints below 0 have the same terms as 0 for delays on the
  // non-negative reals, so they share the entry for 0
  int first = positive ? max(start - pw, 0) : start - pw;
  vector[n] log_cdfs;
  // terms[t - first + 1] holds the terms at endpoint t
  array[n - first + 1] vector[2 * (n_terms + 1)] terms;
  for (t in first:n) {
    terms[t - first + 1] = primarycensored_gumbel_terms(
      t, dist_id, beta, n_terms, params
    );
  }
  for (d in start:n) {
    if (positive && d <= 0) {
      log_cdfs[d] = negative_infinity();
    } else {
      int q_index = (positive ? max(d - pw, 0) : d - pw) - first + 1;
      vector[2] fit = primarycensored_gumbel_lcdf_from_terms(
        terms[d - first + 1], terms[q_index], d, pwindow, mu, beta, n_terms
      );
      if (fit[2] <= gumbel_error_tolerance()) {
        log_cdfs[d] = fit[1];
      } else {
        log_cdfs[d] = primarycensored_gumbel_numeric_lcdf(
          d | dist_id, params, pwindow, mu, beta
        );
      }
    }
  }
  return log_cdfs;
}
