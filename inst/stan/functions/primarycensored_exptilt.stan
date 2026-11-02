// Tilt transforms T_f(xi; tau) = int exp(xi u) f(u) du up to tau, from 0 for
// delays on the non-negative reals. The R equivalents are `.pcens_tilt_*()`.

/**
  * Compute the log CDF and log CCDF of a unit rate Gamma distribution from
  * the log of x
  * @ingroup tilt_transforms
  *
  * Below x = a + 1 the CCDF is `log1m_exp()` of `gamma_lcdf_logx()`. Above,
  * Q(a, x) comes from its continued fraction (modified Lentz) and the CDF
  * from `log1m_exp()` of it. The shape gradient of `gamma_lccdf()` is not
  * used as it can be wrong by 40% or `NaN` for shapes above about 100.
  *
  * @param log_x Log of the argument, log(x)
  * @param a Shape parameter of the Gamma distribution (a > 0)
  *
  * @return Vector [log P(a, exp(log_x)), log Q(a, exp(log_x))]
  */
vector gamma_lcdf_lccdf_logx(real log_x, real a) {
  if (!(a > 0) || is_inf(a)) {
    reject("gamma_lcdf_lccdf_logx: shape must be positive finite, found ", a);
  }
  real x = exp(log_x);
  if (!(x >= a + 1)) {
    real log_p = gamma_lcdf_logx(log_x, a);
    return [log_p, log1m_exp(log_p)]';
  }
  if (is_inf(x)) return [0, negative_infinity()]';
  real b = x + 1 - a;
  real c = 1e300;
  real d = 1 / b;
  real h = d;
  real delta = 1;
  // err bounds |delta - 1|. Flooring |i - a| at 1 keeps it above 0 at an
  // integer shape, where the fraction ends but its shape gradient does not.
  real err = 1;
  int i = 0;
  while (err >= 1e-15) {
    i += 1;
    real an = -i * (i - a);
    real d_prev = d;
    b += 2;
    d = an * d + b;
    if (abs(d) < 1e-300) d = 1e-300;
    c = b + an / c;
    if (abs(c) < 1e-300) c = 1e-300;
    d = 1 / d;
    err *= i * fmax(abs(i - a), 1) * abs(d * d_prev / delta);
    delta = d * c;
    h *= delta;
  }
  real log_q = a * log_x - x - lgamma(a) + log(h);
  return [log1m_exp(log_q), log_q]';
}

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
  * (rate / (rate - xi))^shape, with both tails from gamma_lcdf_lccdf_logx().
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
    real log_total = -shape * log1m(xi / rate);
    if (t <= 0) return [negative_infinity(), log_total]';
    return log_total + gamma_lcdf_lccdf_logx(log(t) + log(rate - xi), shape);
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
  * The log of G_k(t) = int_0^t (t - u)^k f(u) du for k = 1, 2, 3.
  * primarycensored_gamma_uniform_terms() gives [log(t m_0), log m_1], and
  * with the shape raised by 2 the same for m_2 and m_3 up to a constant.
  *
  * @param t Point, positive
  * @param shape Shape
  * @param rate Rate
  *
  * @return Vector [log G_1(t), log G_2(t), log G_3(t)]
  */
vector primarycensored_gamma_tilt_moments(real t, real shape, real rate) {
  real log_t = log(t);
  vector[2] m01 = primarycensored_gamma_uniform_terms(t, {shape, rate});
  vector[2] m23 = primarycensored_gamma_uniform_terms(t, {shape + 2, rate})
                  + log(shape) + log1p(shape) - 2 * log(rate);
  return primarycensored_moments_from_partial(
    log_t, [m01[1] - log_t, m01[2], m23[1] - log_t, m23[2]]'
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
 * endpoint in primarycensored_exptilt_lcdf_vectorized().
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
  return primary_id == 2 && (dist_id == 2 || dist_id == 4 || dist_id == 18);
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
  return log_sum_exp(terms_q[1], log_num - log_den);
}

/**
  * Combine the moments at d and q into the small tilt log CDF
  * @ingroup exponential_tilt_solutions
  *
  * The uniform window limit with its corrections to second order in the tilt,
  * for the small window form of primarycensored_exptilt_form(). With
  * G_k(t) = int (t - u)^k f(u) du,
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
  * Compute the exponentially tilted terms at an endpoint
  * @ingroup exponential_tilt_solutions
  *
  * @param t Endpoint, d or d - pwindow
  * @param dist_id Distribution identifier, see check_for_exptilt()
  * @param params Array of distribution parameters
  * @param rho Tilt
  *
  * @return Vector of the lower and upper terms of log_tilt_transform_pair()
  * for xi = 0 and xi = -rho
  */
vector primarycensored_exptilt_terms(real t, int dist_id, array[] real params,
                                     real rho) {
  return append_row(
    log_tilt_transform_pair(t, dist_id, 0, params),
    log_tilt_transform_pair(t, dist_id, -rho, params)
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
  if (dist_has_positive_support(dist_id) && d <= 0) {
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
  return primarycensored_exptilt_lcdf_from_terms(
    primarycensored_exptilt_terms(d, dist_id, params, rho),
    primarycensored_exptilt_terms(q, dist_id, params, rho), d, rho, pwindow
  );
}

/**
  * Compute the exponentially tilted log CDF at integer delays
  * @ingroup exponential_tilt_solutions
  *
  * The small tilt forms of primarycensored_exptilt_form() only hold below a
  * delay, so they cover a first run of delays, each computed with
  * primarycensored_exptilt_lcdf(). The direct form covers the rest, with the
  * terms at each endpoint computed once and shared by the delays t and
  * t + pwindow. The values are the same as from
  * primarycensored_exptilt_lcdf().
  *
  * @param start First delay to compute, at least 1
  * @param n Last delay to compute, and the length of the result
  * @param dist_id Distribution identifier, see check_for_exptilt()
  * @param params Array of distribution parameters
  * @param pw Primary event window, a positive integer
  * @param rho Tilt
  *
  * @return Vector whose element d is the log CDF at d, for d in start:n.
  * Elements before start are not computed.
  */
vector primarycensored_exptilt_lcdf_vectorized(data int start, data int n,
                                               data int dist_id,
                                               array[] real params,
                                               data int pw, real rho) {
  vector[n] log_cdfs;
  int d0 = start;
  while (d0 <= n && primarycensored_exptilt_form(d0, dist_id, rho, pw)) {
    log_cdfs[d0] = primarycensored_exptilt_lcdf(d0 | dist_id, params, pw, rho);
    d0 += 1;
  }
  if (d0 > n) return log_cdfs;
  // terms[t - d0 + pw + 1] holds the terms at endpoint t
  array[n - d0 + pw + 1] vector[4] terms;
  for (t in (d0 - pw):n) {
    terms[t - d0 + pw + 1] = primarycensored_exptilt_terms(
      t, dist_id, params, rho
    );
  }
  for (d in d0:n) {
    log_cdfs[d] = primarycensored_exptilt_lcdf_from_terms(
      terms[d - d0 + pw + 1], terms[d - d0 + 1], d, rho, pw
    );
  }
  return log_cdfs;
}
