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
  * Relative error of the derivative of gamma_lcdf() and gamma_lccdf() in the
  * shape
  * @ingroup truncated_gumbel_solutions
  *
  * The derivative of the regularised incomplete gamma function in the shape
  * is computed by Stan with a relative error of 1e-3 to 1e-2, for example
  * 0.8% for the upper tail at a shape of 3 and a point of 8. A sum of
  * alternating terms multiplies the error of the derivative of each term by
  * the ratio of the sum of the absolute terms to the sum, so the series is
  * not used for the gamma where that would give an error above
  * gumbel_gradient_tolerance(). The derivatives in the other parameters,
  * and of the exponential and the normal, are exact, so they are limited by
  * the rounding error in gumbel_error_tolerance().
  *
  * @return 1e-2
  */
real gumbel_gamma_shape_gradient_error() {
  return 1e-2;
}

/**
  * Largest estimated relative error of the gradient for which the Gumbel
  * series is used for the gamma
  * @ingroup truncated_gumbel_solutions
  *
  * See gumbel_gamma_shape_gradient_error(). The numerical path has the error
  * of the derivative of the delay without amplification.
  *
  * @return 2e-2
  */
real gumbel_gradient_tolerance() {
  return 2e-2;
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
  * Check if the Gumbel series is accurate enough at a delay
  * @ingroup truncated_gumbel_solutions
  *
  * The series is used if its estimated relative error is at most
  * gumbel_error_tolerance(), and for the gamma delay, whose derivative in
  * the shape is approximate, if the amplification of that error is at most
  * gumbel_gradient_tolerance() over gumbel_gamma_shape_gradient_error().
  *
  * @param dist_id Distribution identifier
  * @param fit Vector from primarycensored_gumbel_lcdf_from_terms()
  *
  * @return 1 if the series is used, 0 if the numerical path is
  */
int gumbel_series_accepted(int dist_id, vector fit) {
  if (!(fit[2] <= gumbel_error_tolerance())) return 0;
  if (dist_id == 2
      && !(fit[3] * gumbel_gamma_shape_gradient_error()
           <= gumbel_gradient_tolerance())) {
    return 0;
  }
  return 1;
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
  * The amplification is the sum of the absolute terms over B, which
  * multiplies the error of the derivative of each term in the gradient, see
  * gumbel_gamma_shape_gradient_error().
  *
  * @return Vector [log of the primary event censored CDF at d, estimated
  * relative error, amplification]
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
    return [log_f_q, 0, 1]';
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
  real amplification;
  if (log_bracket == negative_infinity()) {
    error = positive_infinity();
    amplification = positive_infinity();
  } else {
    int n_last = n_terms + 1;
    real log_rounding = log(machine_precision()) + log1p(scale)
                        + gumbel_log_sum_exp(log_pos, log_neg) - log_bracket;
    real log_trunc = log_delta0 + n_last * log_s0 - lgamma(n_last + 1)
                     - log1m(fmin(exp(log_s0) / (n_last + 1), 0.5))
                     - log_bracket;
    error = exp(log_rounding) + exp(log_trunc);
    amplification = exp(gumbel_log_sum_exp(log_pos, log_neg) - log_bracket);
  }
  real log_cdf = gumbel_log_sum_exp(log_f_q, log_bracket - log_d);
  return [fmin(log_cdf, 0), error, amplification]';
}

/**
  * Length in u of a piece of the numerical integral
  * @ingroup truncated_gumbel_solutions
  *
  * With u = s(z) - s(pwindow) and the window density exp(-u) / (1 -
  * exp(-Delta)), the quantile variable v = 1 - exp(-(u - a)) of a piece that
  * starts at u = a is at most 1 - 2^-53 in double precision, which is a
  * length of about 36.7. A piece is at most 36 long, so the mass beyond it
  * (2.3e-16 of the piece) is not resolved and is given the value at 36, and
  * the next piece takes over from there, see gumbel_numeric_spike_lcdf().
  *
  * @return 36
  */
real gumbel_numeric_u_max() {
  return 36;
}

/**
  * Limits of the numerical integral in z of a truncated Gumbel primary
  * @ingroup truncated_gumbel_solutions
  *
  * For mu below the end of the window, the window is integrated in z over the
  * range that holds its mass, from where s(z) - s(pwindow) is 745 or more, to
  * 46 beta above the larger of mu and 0, cut at d for delays on the
  * non-negative reals. The mass outside is below 1e-16 of the window.
  *
  * @param d Delay
  * @param dist_id Distribution identifier
  * @param pwindow Primary event window
  * @param mu Location of the truncated Gumbel
  * @param beta Scale of the truncated Gumbel
  *
  * @return Vector [lower, upper] of limits of z. The integral is 0 if upper
  * is not above lower.
  */
vector gumbel_numeric_z_limits(data real d, int dist_id, data real pwindow,
                               real mu, real beta) {
  real log_sw = (mu - pwindow) / beta;
  real z_low = fmax(0, mu - beta * log(745 + exp(log_sw)));
  real z_high = fmin(pwindow, fmax(mu, 0) + 46 * beta);
  if (dist_has_positive_support(dist_id)) {
    z_high = fmin(z_high, d);
  }
  return [z_low, z_high]';
}

/**
  * Log CDF of the delay for the numerical integral, with finite partials
  * @ingroup truncated_gumbel_solutions
  *
  * The same as dist_lcdf(), except that the exponential, gamma and normal
  * forms are chosen so that the partials are finite where the CDF is very
  * small, which is where the integrand starts. The solver evaluates
  * sensitivities there, and a partial that is not finite on the tape gives a
  * NaN gradient even where the result is multiplied by 0. The exponential is
  * log(1 - exp(-rate x)) as log1m_exp(), which is exact for a small rate x
  * where exponential_lcdf() is `-inf`. The gamma is `-inf` without calling
  * gamma_lcdf() where its lower tail underflows, see
  * gamma_lcdf_underflows(). The normal is primarycensored_log_std_normal_cdf()
  * of the standardised point, which has an exact derivative for any point.
  *
  * @param x Delay
  * @param params Array of distribution parameters
  * @param dist_id Distribution identifier
  *
  * @return Log of the delay CDF at x
  */
real gumbel_delay_lcdf(real x, array[] real params, int dist_id) {
  if (dist_has_positive_support(dist_id) && x <= 0) {
    return negative_infinity();
  }
  if (dist_id == 4) {
    return log1m_exp(-params[1] * x);
  } else if (dist_id == 2) {
    if (gamma_lcdf_underflows(x * params[2], params[1])) {
      return negative_infinity();
    }
    return gamma_lcdf(x | params[1], params[2]);
  } else if (dist_id == 18) {
    return primarycensored_log_std_normal_cdf((x - params[1]) / params[2]);
  }
  return dist_lcdf(x | params, dist_id);
}

/**
  * Delay at a point u of the window for mu at or above the window end
  * @ingroup truncated_gumbel_solutions
  *
  * With u = s(z) - s(pwindow) the delay is d - z = (d - pwindow) + beta
  * log(1 + u / s(pwindow)), which keeps its small difference from 0 where
  * the spike is within rounding of pwindow.
  *
  * @param u Point in u, at least 0
  * @param d Delay
  * @param pwindow Primary event window
  * @param mu Location of the truncated Gumbel, at least pwindow
  * @param beta Scale of the truncated Gumbel
  *
  * @return The delay d - z
  */
real gumbel_numeric_delay(real u, data real d, data real pwindow, real mu,
                          real beta) {
  return (d - pwindow) + beta * log1p(u * exp(-(mu - pwindow) / beta));
}

/**
  * Delay argument and log density weight of the numerical integral
  * @ingroup truncated_gumbel_solutions
  *
  * The integration variable tau is in [0, 1] and is mapped to the limits. The
  * window density is a narrow spike when mu / beta is large, which an
  * adaptive integrator can step over, so two variables are used, both with a
  * smooth integrand.
  *
  * With mu at or above the window end, s(pwindow) = exp((mu - pwindow) /
  * beta) is at least 1 and the integral is taken in the upper quantile v of
  * a piece of the window that starts at u = a, see
  * gumbel_numeric_spike_lcdf(). The limits are 0 and the largest v, the
  * density of v is 1 and the log weight is 0, and the delay is that at
  * u = a - log(1 - v), see gumbel_numeric_delay().
  *
  * With mu below it the integral is in z, with the limits of
  * gumbel_numeric_z_limits(), the delay is d - z and the log weight is the
  * log window density at z.
  *
  * @param tau Integration variable in [0, 1]
  * @param lo Lower limit of the variable
  * @param hi Upper limit of the variable
  * @param a Start of the piece in u, for mu at or above the window end
  * @param d Delay
  * @param pwindow Primary event window
  * @param mu Location of the truncated Gumbel
  * @param beta Scale of the truncated Gumbel
  *
  * @return Vector [delay, log weight]
  */
vector gumbel_numeric_point(real tau, real lo, real hi, real a,
                            data real d, data real pwindow, real mu,
                            real beta) {
  real x = lo + (hi - lo) * tau;
  if (mu >= pwindow) {
    return [gumbel_numeric_delay(a - log1m(x), d, pwindow, mu, beta), 0]';
  }
  return [d - x, tgumbel_lpdf(x | 0, pwindow, mu, beta)]';
}

/**
  * ODE system for the truncated Gumbel primary event censored CDF
  * @ingroup truncated_gumbel_solutions
  *
  * The derivative with respect to tau of the integral of the delay CDF, see
  * gumbel_numeric_point(), scaled by exp(-log_shift) so that the integrand is
  * of order 1, and so the solver tolerances are relative to a small CDF.
  *
  * @param tau Integration variable in [0, 1]
  * @param y State, the integral so far
  * @param d Delay
  * @param pwindow Primary event window
  * @param dist_id Distribution identifier
  * @param params Array of distribution parameters
  * @param mu Location of the truncated Gumbel
  * @param beta Scale of the truncated Gumbel
  * @param lo Lower limit of the integration variable, see
  *   gumbel_numeric_point()
  * @param hi Upper limit of the integration variable
  * @param a Start of the piece in u
  * @param log_shift Log of the scale of the integrand
  *
  * @return The derivative of the state with respect to tau
  */
vector primarycensored_gumbel_ode(real tau, vector y, data real d,
                                  data real pwindow, int dist_id,
                                  array[] real params, real mu, real beta,
                                  real lo, real hi, real a,
                                  real log_shift) {
  vector[2] point = gumbel_numeric_point(
    tau, lo, hi, a, d, pwindow, mu, beta
  );
  return rep_vector(
    (hi - lo)
    * exp(gumbel_delay_lcdf(point[1] | params, dist_id) + point[2] - log_shift),
    1
  );
}

/**
  * Log of the integral of one piece of the window in u
  * @ingroup truncated_gumbel_solutions
  *
  * The integral of the delay CDF against exp(-u) over u in [a, a + len], in
  * the upper quantile v = 1 - exp(-(u - a)) of the piece, so that exp(-a) is
  * carried on the log scale, and the integrand of v is the delay CDF at the
  * delay of gumbel_numeric_point(), which is increasing in u. It is scaled by
  * the delay CDF at the end of the piece, so the solver tolerances of a
  * relative 1e-12 and an absolute 1e-18 hold relative to the piece.
  *
  * The largest v is c = 1 - exp(-len). It is an argument and not computed
  * here because an argument that does not depend on the parameters is not
  * differentiated by the solver. The sensitivity of the integral to c is
  * large for c near 1, from the change of variable, so pieces of a length
  * that does not depend on the parameters are passed as data, see
  * gumbel_numeric_spike_lcdf().
  *
  * @param d Delay
  * @param dist_id Distribution identifier
  * @param params Array of distribution parameters
  * @param pwindow Primary event window
  * @param mu Location of the truncated Gumbel
  * @param beta Scale of the truncated Gumbel
  * @param a Start of the piece in u
  * @param c Largest upper quantile of the piece, 1 - exp(-len)
  * @param log_shift Log of the delay CDF at the end of the piece
  *
  * @return Log of the integral of the delay CDF times exp(-u) over the piece,
  * `-inf` if it is 0
  */
real gumbel_numeric_piece_lcdf(data real d, int dist_id, array[] real params,
                               data real pwindow, real mu, real beta, real a,
                               real c, real log_shift) {
  if (log_shift == negative_infinity()) {
    return negative_infinity();
  }
  vector[1] y0 = rep_vector(0.0, 1);
  real integral = ode_rk45_tol(
    primarycensored_gumbel_ode, y0, 0.0, {1.0}, 1e-12, 1e-18, 1000000,
    d, pwindow, dist_id, params, mu, beta, 0.0, c, a, log_shift
  )[1, 1];
  if (!(integral > 0)) {
    return negative_infinity();
  }
  return -a + log_shift + log(integral);
}

/**
  * Log of the numerical CDF for a truncated Gumbel primary at or above the
  * window end
  * @ingroup truncated_gumbel_solutions
  *
  * With mu at or above the window end the window is integrated in u = s(z) -
  * s(pwindow), where the density is exp(-u) / (1 - exp(-Delta)) on [0,
  * Delta], with Delta = s(0) - s(pwindow), in pieces of at most
  * gumbel_numeric_u_max(), see gumbel_numeric_piece_lcdf(). The first piece
  * starts at u = 0, or at the point where the delay CDF leaves 0 for a delay
  * on the non-negative reals with d below pwindow, and the weight of a piece
  * is on the log scale so a CDF of 1e-300 or less is not lost.
  *
  * The integrand of a piece in its quantile variable is the delay CDF, which
  * is smooth if the CDF changes by a factor of at most exp(3) over the second
  * half of the piece. A piece is halved until it does, so where the delay CDF
  * rises faster than exp(-u) falls, as for a normal in its lower tail and a
  * mass that is far from the start of the window, the pieces are short where
  * it rises and the mass is integrated where it is. The log delay CDF is
  * concave in u, so after a piece whose log CDF rises by less than its
  * length, the rest is at most the end value over one minus that slope. The
  * pieces stop where this is below 1e-15 of the total, and at the end of the
  * window. At most 300 pieces are used, and more is an error.
  *
  * The pieces are 36 long while more than 48 of the window is left, and
  * otherwise 12 long until at most 12 is left, where the last piece ends at
  * the end of the window. Its length depends on the parameters, and the
  * others do not, see gumbel_numeric_piece_lcdf().
  *
  * @param d Delay
  * @param dist_id Distribution identifier
  * @param params Array of distribution parameters
  * @param pwindow Primary event window
  * @param mu Location of the truncated Gumbel, at least pwindow
  * @param beta Scale of the truncated Gumbel
  *
  * @return Log of the CDF, `-inf` if it is 0
  */
real gumbel_numeric_spike_lcdf(data real d, int dist_id, array[] real params,
                               data real pwindow, real mu, real beta) {
  real log_delta = (mu - pwindow) / beta + log_diff_exp(pwindow / beta, 0);
  real a = 0;
  if (dist_has_positive_support(dist_id) && d < pwindow) {
    // The delay CDF is 0 for u below the kink, where d - z is 0
    real log_uk = (mu - pwindow) / beta
                  + log_diff_exp((pwindow - d) / beta, 0);
    // Beyond u of 1e13 the log CDF is below -1e13, which is 0 in double
    // precision, and u + 36 is not distinguishable from u from 1e16
    if (log_uk >= log_delta || log_uk > 30) {
      return negative_infinity();
    }
    a = exp(log_uk);
  }
  real top = exp(log_delta);
  real log_total = negative_infinity();
  int finished = 0;
  for (k in 1:300) {
    // The length of a piece that is not the last is 36 or 12 over 2^j
    int big = top - a > 48;
    int j = 0;
    int last = top - a <= 12;
    real len = last ? top - a : (big ? 36.0 : 12.0);
    real log_f_end = gumbel_delay_lcdf(
      gumbel_numeric_delay(a + len, d, pwindow, mu, beta) | params, dist_id
    );
    real log_f_mid = gumbel_delay_lcdf(
      gumbel_numeric_delay(a + len / 2, d, pwindow, mu, beta)
      | params, dist_id
    );
    for (i in 1:30) {
      if (!(log_f_end - log_f_mid > 3)) {
        break;
      }
      if (last) {
        // Continue with pieces that are at most half of what is left
        last = 0;
        while (12.0 / pow(2, j) > len / 2) {
          j += 1;
        }
      } else {
        j += 1;
      }
      len = (big ? 36.0 : 12.0) / pow(2, j);
      log_f_end = gumbel_delay_lcdf(
        gumbel_numeric_delay(a + len, d, pwindow, mu, beta) | params, dist_id
      );
      log_f_mid = gumbel_delay_lcdf(
        gumbel_numeric_delay(a + len / 2, d, pwindow, mu, beta)
        | params, dist_id
      );
    }
    if (last) {
      len = top - a;
      log_total = gumbel_log_sum_exp(
        log_total,
        gumbel_numeric_piece_lcdf(
          d | dist_id, params, pwindow, mu, beta, a, -expm1(-(top - a)),
          log_f_end
        )
      );
    } else {
      log_total = gumbel_log_sum_exp(
        log_total,
        gumbel_numeric_piece_lcdf(
          d | dist_id, params, pwindow, mu, beta, a,
          -expm1(-((big ? 36.0 : 12.0) / pow(2, j))), log_f_end
        )
      );
    }
    a += len;
    real slope = 2 * (log_f_end - log_f_mid) / len;
    if (last || a >= top
        || (slope < 0.99
            && -a + log_f_end - log1m(slope) < log_total + log(1e-15))) {
      finished = 1;
      break;
    }
  }
  if (!finished) {
    reject(
      "The truncated Gumbel numerical integral did not finish for d ", d,
      ", pwindow ", pwindow, ", mu ", mu, " and beta ", beta
    );
  }
  return log_total - tgumbel_log1m_exp_neg_exp(log_delta);
}

/**
  * Log of the integrand of the numerical CDF in z
  * @ingroup truncated_gumbel_solutions
  *
  * The log of the delay CDF at d - z plus the log window density at z.
  *
  * @param z Point in the window
  * @param d Delay
  * @param dist_id Distribution identifier
  * @param params Array of distribution parameters
  * @param pwindow Primary event window
  * @param mu Location of the truncated Gumbel
  * @param beta Scale of the truncated Gumbel
  *
  * @return Log of the integrand at z
  */
real gumbel_z_log_integrand(real z, data real d, int dist_id,
                            array[] real params, data real pwindow, real mu,
                            real beta) {
  return gumbel_delay_lcdf(d - z | params, dist_id)
         + tgumbel_lpdf(z | 0, pwindow, mu, beta);
}

/**
  * Integral of the numerical CDF in z from the peak to one end of the range
  * @ingroup truncated_gumbel_solutions
  *
  * The integrand scaled by exp(-log_shift), from z_from to z_to, with the
  * sign of z_to - z_from.
  *
  * @param d Delay
  * @param dist_id Distribution identifier
  * @param params Array of distribution parameters
  * @param pwindow Primary event window
  * @param mu Location of the truncated Gumbel
  * @param beta Scale of the truncated Gumbel
  * @param z_from Start of the integral
  * @param z_to End of the integral
  * @param log_shift Log of the scale of the integrand
  *
  * @return The scaled integral
  */
real gumbel_numeric_z_half(data real d, int dist_id, array[] real params,
                           data real pwindow, real mu, real beta,
                           real z_from, real z_to, real log_shift) {
  vector[1] y0 = rep_vector(0.0, 1);
  return ode_rk45_tol(
    primarycensored_gumbel_ode, y0, 0.0, {1.0}, 1e-12, 1e-18, 1000000,
    d, pwindow, dist_id, params, mu, beta, z_from, z_to, 0.0, log_shift
  )[1, 1];
}

/**
  * Log of the numerical CDF for a truncated Gumbel primary below the window
  * end
  * @ingroup truncated_gumbel_solutions
  *
  * With mu below the window end the window is integrated in z over the limits
  * of gumbel_numeric_z_limits(). The log integrand, see
  * gumbel_z_log_integrand(), is concave in z for the delays with a log
  * concave CDF, so it has one peak, and it can be far narrower than the range
  * and far below the log of 1, for a delay CDF and a window density that are
  * large in different places. The peak is found on 17 equally spaced points
  * and then by golden section search within the two intervals around the
  * best point. The integral is taken from the peak to each end of the range
  * as two integrals, each starting from the peak, where the integrand is
  * exp(log_shift) and decays, scaled by it so that the solver tolerances of a
  * relative 1e-12 and an absolute 1e-18 hold relative to the peak.
  *
  * @param d Delay
  * @param dist_id Distribution identifier
  * @param params Array of distribution parameters
  * @param pwindow Primary event window
  * @param mu Location of the truncated Gumbel, below pwindow
  * @param beta Scale of the truncated Gumbel
  *
  * @return Log of the CDF, `-inf` if it is 0
  */
real gumbel_numeric_z_lcdf(data real d, int dist_id, array[] real params,
                           data real pwindow, real mu, real beta) {
  vector[2] z_limits = gumbel_numeric_z_limits(d, dist_id, pwindow, mu, beta);
  if (!(z_limits[2] > z_limits[1])) {
    return negative_infinity();
  }
  real z_low = z_limits[1];
  real z_high = z_limits[2];
  real step = (z_high - z_low) / 16;
  real z_peak = z_low;
  real log_shift = gumbel_z_log_integrand(
    z_low, d, dist_id, params, pwindow, mu, beta
  );
  for (i in 1:16) {
    real z = z_low + step * i;
    real g = gumbel_z_log_integrand(z, d, dist_id, params, pwindow, mu, beta);
    if (g > log_shift) {
      log_shift = g;
      z_peak = z;
    }
  }
  if (log_shift == negative_infinity()) {
    return negative_infinity();
  }
  // Golden section search in the two intervals around the best point
  real inv_phi = 0.6180339887498949;
  real lo = fmax(z_low, z_peak - step);
  real hi = fmin(z_high, z_peak + step);
  real z1 = hi - inv_phi * (hi - lo);
  real z2 = lo + inv_phi * (hi - lo);
  real g1 = gumbel_z_log_integrand(z1, d, dist_id, params, pwindow, mu, beta);
  real g2 = gumbel_z_log_integrand(z2, d, dist_id, params, pwindow, mu, beta);
  for (i in 1:40) {
    if (g1 > g2) {
      hi = z2;
      z2 = z1;
      g2 = g1;
      z1 = hi - inv_phi * (hi - lo);
      g1 = gumbel_z_log_integrand(z1, d, dist_id, params, pwindow, mu, beta);
    } else {
      lo = z1;
      z1 = z2;
      g1 = g2;
      z2 = lo + inv_phi * (hi - lo);
      g2 = gumbel_z_log_integrand(z2, d, dist_id, params, pwindow, mu, beta);
    }
  }
  if (g1 > log_shift && g1 >= g2) {
    log_shift = g1;
    z_peak = z1;
  } else if (g2 > log_shift) {
    log_shift = g2;
    z_peak = z2;
  }
  real integral = 0;
  if (z_peak > z_low) {
    integral += -gumbel_numeric_z_half(
      d, dist_id, params, pwindow, mu, beta, z_peak, z_low, log_shift
    );
  }
  if (z_high > z_peak) {
    integral += gumbel_numeric_z_half(
      d, dist_id, params, pwindow, mu, beta, z_peak, z_high, log_shift
    );
  }
  if (!(integral > 0)) {
    return negative_infinity();
  }
  return log_shift + log(integral);
}

/**
  * Log of the numerical primary event censored CDF for a truncated Gumbel
  * primary, not capped at 0
  * @ingroup truncated_gumbel_solutions
  *
  * The integral of the delay CDF against the window density with an ODE
  * solver, in a variable in which the integrand is smooth for every mu and
  * beta, see gumbel_numeric_spike_lcdf() and gumbel_numeric_z_lcdf(). It is
  * evaluated on the log scale, so a CDF much below 1e-300 is not lost.
  *
  * @param d Delay
  * @param dist_id Distribution identifier
  * @param params Array of distribution parameters
  * @param pwindow Primary event window
  * @param mu Location of the truncated Gumbel
  * @param beta Scale of the truncated Gumbel
  *
  * @return Log of the CDF, which is above 0 by at most the error of the
  * solver, and `-inf` where the CDF is 0
  */
real gumbel_numeric_log_cdf(data real d, int dist_id, array[] real params,
                            data real pwindow, real mu, real beta) {
  if (dist_has_positive_support(dist_id) && d <= 0) {
    return negative_infinity();
  }
  if (mu >= pwindow) {
    return gumbel_numeric_spike_lcdf(d | dist_id, params, pwindow, mu, beta);
  }
  return gumbel_numeric_z_lcdf(d | dist_id, params, pwindow, mu, beta);
}

/**
  * Primary event censored CDF for a truncated Gumbel primary by numerical
  * integration
  * @ingroup truncated_gumbel_solutions
  *
  * Integrates the delay CDF against the window density with an ODE solver
  * in a variable in which the integrand is smooth for every mu and beta, see
  * gumbel_numeric_log_cdf(). It is primarycensored_numeric_cdf() for the
  * truncated Gumbel primary, whose integration over the window would miss a
  * spike of the density. The result was within 1e-11 of a tight reference in
  * tests, relative to the CDF. It is in [0, 1].
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
  real log_cdf = gumbel_numeric_log_cdf(d | dist_id, params, pwindow, mu, beta);
  if (is_nan(log_cdf)) {
    return 0;
  }
  return exp(fmin(log_cdf, 0));
}

/**
  * Log of the numerical primary event censored CDF for a truncated Gumbel
  * primary
  * @ingroup truncated_gumbel_solutions
  *
  * Evaluated on the log scale, see gumbel_numeric_log_cdf(). An integral that
  * is not above 0, which is below the solver tolerance, gives `-inf` and one
  * above what makes the CDF 1 is capped at 0, so the result is never NaN and
  * never above 0.
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
  real log_cdf = gumbel_numeric_log_cdf(d | dist_id, params, pwindow, mu, beta);
  if (is_nan(log_cdf)) {
    return negative_infinity();
  }
  return fmin(log_cdf, 0);
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
  vector[3] fit = primarycensored_gumbel_lcdf_from_terms(
    primarycensored_gumbel_terms(d, dist_id, beta, n_terms, params),
    primarycensored_gumbel_terms(d - pwindow, dist_id, beta, n_terms, params),
    d, pwindow, mu, beta, n_terms
  );
  if (gumbel_series_accepted(dist_id, fit)) {
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
      vector[3] fit = primarycensored_gumbel_lcdf_from_terms(
        terms[d - first + 1], terms[q_index], d, pwindow, mu, beta, n_terms
      );
      if (gumbel_series_accepted(dist_id, fit)) {
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
