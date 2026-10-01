/*
 * Truncated Gumbel primary event window
 *
 * With the density of tgumbel.stan on [0, w], G(z) = exp(-exp(-(z - mu) /
 * beta)), D = G(w) - G(0), q = d - w, F the delay CDF and T_f(xi; t) the tilt
 * transform of tilt_transform.stan,
 *   F_G(d) = F(q) + B(d) / D,
 *   B(d) = (1 - G(0)) (F(d) - F(q))
 *          + sum_{n >= 1} (-c)^n / n! (T_f(n / beta; d) - T_f(n / beta; q)),
 * with c = exp(-(d - mu) / beta).
 * Every term depends on d or q alone, so the terms are computed once per
 * endpoint in primarycensored_analytical_lcdf_vectorized().
 * The series alternates, so the terms are summed on the log scale with the
 * odd and even terms apart, and it is used for normal delays where mu / beta
 * is below log(15) and its estimated relative error is below 1e-8.
 * Otherwise the numerical path of primarycensored_gumbel_numeric_cdf() is
 * used, which integrates in a variable where a narrow window density is
 * smooth.
 */

/**
  * Number of terms of the Gumbel series
  * @ingroup truncated_gumbel_solutions
  *
  * The smallest N, at least 2, with s0^n / n! below 1e-18 max(s0, 1) after
  * N, where s0 = exp(mu / beta).
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
  * The truncated Gumbel primary is primary_id 4 with primary_params =
  * [mu, beta]. The series is for the normal delay. Whether it applies for
  * given parameters is check_for_analytical_params().
  *
  * @param dist_id Distribution identifier for the delay distribution
  * @param primary_id Distribution identifier for the primary distribution
  *
  * @return 1 if the delay has a truncated Gumbel solution and the primary is
  *   the truncated Gumbel, 0 otherwise
  */
int check_for_gumbel(int dist_id, int primary_id) {
  return primary_id == 4 && dist_id == 18;
}

/**
  * Check if the Gumbel series applies for the given parameters
  * @ingroup truncated_gumbel_solutions
  *
  * Needs a normal delay, mu / beta below log(15) and a rounding estimate of
  * at most 1e-7, from the size of the log transform at the largest tilt.
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
  if (!check_for_gumbel(dist_id, 4)) return 0;
  if (!(beta > 0) || is_inf(beta) || is_nan(mu) || is_inf(mu)) return 0;
  real log_s0 = mu / beta;
  if (log_s0 > log(15)) return 0;
  real xi = gumbel_n_terms(log_s0) / beta;
  real size = xi * abs(params[1]) + 0.5 * square(xi * params[2]);
  return log(machine_precision()) + log1p(size) + exp(log_s0)
         <= log(1e-7);
}

/**
  * Check if the Gumbel series is accurate enough at a delay
  * @ingroup truncated_gumbel_solutions
  *
  * Needs an estimated relative error of at most 1e-8.
  *
  * @param fit Vector from primarycensored_gumbel_lcdf_from_terms()
  *
  * @return 1 if the series is used, 0 if the numerical path is
  */
int gumbel_series_accepted(vector fit) {
  return fit[2] <= 1e-8;
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
  * @return abs(x), or 0 if x is not finite
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
  * @return Vector of length 2 (n_terms + 1) holding, for each tilt
  *   xi = n / beta, n = 0, ..., n_terms, the pair [log T_f(xi; t),
  *   log(T_f(xi; Inf) - T_f(xi; t))] of log_tilt_transform_pair(). Only
  *   defined where check_for_gumbel_params() is 1.
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
  * Each difference of transforms is taken between the lower or the upper
  * tail terms, see primarycensored_tail_diff(). The estimated relative error
  * is the machine precision times the sum of the absolute terms over the
  * bracket and the size of the log terms, plus the truncation bound.
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
  *   relative error]
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
  // log D = -s_w + log(1 - exp(-delta)), delta = s_w (exp(w / beta) - 1)
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
  * Limits of the numerical integral in z for mu below the window end
  * @ingroup truncated_gumbel_solutions
  *
  * The range that holds the mass of the window, cut at d for delays on the
  * non-negative reals.
  *
  * @param d Delay
  * @param dist_id Distribution identifier
  * @param pwindow Primary event window
  * @param mu Location of the truncated Gumbel
  * @param beta Scale of the truncated Gumbel
  *
  * @return Vector [lower, upper] of limits of z
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
  * Upper limit of the numerical integral in z that is data
  * @ingroup truncated_gumbel_solutions
  *
  * The window end, cut at d for delays on the non-negative reals. It is
  * passed to the solver in place of the limit to keep it out of the
  * sensitivities, which are singular at d for a delay CDF with an infinite
  * slope at 0.
  *
  * @param d Delay
  * @param dist_id Distribution identifier
  * @param pwindow Primary event window
  *
  * @return The upper limit of z that is data
  */
real gumbel_numeric_z_data_limit(data real d, int dist_id,
                                 data real pwindow) {
  if (dist_has_positive_support(dist_id)) {
    return fmin(pwindow, d);
  }
  return pwindow;
}

/**
  * Log CDF of the delay for the numerical integral, with finite partials
  * @ingroup truncated_gumbel_solutions
  *
  * The same as dist_lcdf(), with forms of the lognormal, exponential, gamma
  * and normal whose partials are finite where the CDF is very small.
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
  if (dist_id == 1) {
    if (lognormal_lcdf_underflows(x, params[1], params[2])) {
      return negative_infinity();
    }
    return primarycensored_log_std_normal_cdf(
      (log(x) - params[1]) / params[2]
    );
  } else if (dist_id == 4) {
    return log1m_exp(-params[1] * x);
  } else if (dist_id == 2) {
    if (gamma_lcdf_underflows(x * params[2], params[1])) {
      return negative_infinity();
    }
    return gamma_lcdf(x | params[1], params[2]);
  } else if (dist_id == 18) {
    return primarycensored_log_std_normal_cdf((x - params[1]) / params[2]);
  } else if (dist_id == 3) {
    if (params[1] * (log(x) - log(params[2])) < -700) {
      return negative_infinity();
    }
  } else if (dist_id == 13) {
    if (gamma_lcdf_underflows(x / 2, params[1] / 2)) {
      return negative_infinity();
    }
    return gamma_lcdf(x | params[1] / 2, 0.5);
  } else if (dist_id == 16 || dist_id == 19 || dist_id == 22) {
    // Inverse gamma, with the shape and scale of each family
    real shape = dist_id == 16 ? params[1] : params[1] / 2;
    real scale = dist_id == 16 ? params[2]
                 : (dist_id == 19 ? 0.5 : params[1] * square(params[2]) / 2);
    if (gamma_lccdf_underflows(scale / x, shape)) {
      return negative_infinity();
    }
    return inv_gamma_lcdf(x | shape, scale);
  } else if (dist_id == 9) {
    if (x >= 1) {
      return 0;
    }
    if (x < params[1] / (params[1] + params[2])
        && params[1] * log(x) - log(params[1]) - lbeta(params[1], params[2])
           < -700) {
      return negative_infinity();
    }
  }
  return dist_lcdf(x | params, dist_id);
}

/**
  * Delay at a point u of the window for mu at or above the window end
  * @ingroup truncated_gumbel_solutions
  *
  * The delay d - z = (d - pwindow) + beta log(1 + u / s(pwindow)) keeps its
  * small difference from 0 where the spike is within rounding of pwindow.
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
  * The integration variable tau in [0, 1] is mapped to [lo, hi]. For mu at
  * or above the window end the integral is in the upper quantile of a piece
  * that starts at u = a, see gumbel_numeric_spike_lcdf(), with log weight 0.
  * For mu below it the integral is in z, weighted by the log window density.
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
  * The derivative in tau of the integral of the delay CDF, see
  * gumbel_numeric_point(), scaled by exp(-log_shift) so that the integrand
  * is of order 1 and the solver tolerances are relative.
  *
  * @param tau Integration variable in [0, 1]
  * @param y State, the integral so far
  * @param d Delay
  * @param pwindow Primary event window
  * @param dist_id Distribution identifier
  * @param params Array of distribution parameters
  * @param mu Location of the truncated Gumbel
  * @param beta Scale of the truncated Gumbel
  * @param lo Lower limit of the integration variable
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
  * the upper quantile v = 1 - exp(-(u - a)), scaled by the delay CDF at the
  * end of the piece. The largest v, c = 1 - exp(-len), is an argument so
  * that a length that does not depend on the parameters is not
  * differentiated.
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
  * @return Log of the integral of the delay CDF times exp(-u) over the
  *   piece, `-inf` if it is 0
  */
real gumbel_numeric_piece_lcdf(data real d, int dist_id, array[] real params,
                               data real pwindow, real mu, real beta, real a,
                               real c, real log_shift) {
  if (log_shift == negative_infinity()) {
    return negative_infinity();
  }
  vector[1] y0 = rep_vector(0.0, 1);
  real integral = ode_rk45_tol(
    primarycensored_gumbel_ode, y0, 0.0, {1.0}, 1e-10, 1e-14, 200000,
    d, pwindow, dist_id, params, mu, beta, 0.0, c, a, log_shift
  )[1, 1];
  if (!(integral > 0)) {
    return negative_infinity();
  }
  return -a + log_shift + log(integral);
}

/**
  * Log of the numerical CDF for mu at or above the window end
  * @ingroup truncated_gumbel_solutions
  *
  * The window is integrated in u = s(z) - s(pwindow), with density
  * exp(-u) / (1 - exp(-Delta)) on [0, Delta], in pieces of at most
  * 36, from the first u with a delay CDF above 0.
  * A piece is halved until the delay CDF changes by a factor of at most
  * exp(3) over its second half. The log delay CDF is concave in u, so the
  * pieces stop where the rest is below 1e-15 of the total. At most 300 are
  * used. The result is `-inf` where the delay CDF is 0 for u below 1000, as
  * the log CDF is then below -1000.
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
    // The CDF is below exp(-1000) beyond u of 1000
    if (log_uk >= log_delta || log_uk > log(1000)) {
      return negative_infinity();
    }
    a = exp(log_uk);
  }
  real top = exp(log_delta);
  // A piece of at most 36 ends within 2^-53 of 1 in the upper quantile.
  // The lengths are literals so that they are not differentiated.
  real log_total = negative_infinity();
  int finished = 0;
  for (k in 1:300) {
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
    // The CDF is below exp(-1000) beyond u of 1000
    if (log_total == negative_infinity() && a > 1000) {
      return negative_infinity();
    }
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
      ", pwindow ", pwindow, ", mu ", mu, " and beta ", beta,
      ". The log CDF is extremely small, below about -1000, so the delay ",
      "distribution and the window are far apart."
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
    primarycensored_gumbel_ode, y0, 0.0, {1.0}, 1e-10, 1e-14, 200000,
    d, pwindow, dist_id, params, mu, beta, z_from, z_to, 0.0, log_shift
  )[1, 1];
}

/**
  * Log of the numerical CDF for mu below the window end
  * @ingroup truncated_gumbel_solutions
  *
  * The window is integrated in z over the limits of
  * gumbel_numeric_z_limits(). The log integrand is concave in z for delays
  * with a log concave CDF, and can be far narrower than the range. The peak
  * is found on a grid and by golden section search, and the integral is
  * taken from it to each end of the range, scaled by the integrand at the
  * peak.
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
  // Golden section search around the best grid point
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
    // An upper limit at the window end or d is data
    if (z_high >= gumbel_numeric_z_data_limit(d, dist_id, pwindow)) {
      integral += gumbel_numeric_z_half(
        d, dist_id, params, pwindow, mu, beta, z_peak,
        gumbel_numeric_z_data_limit(d, dist_id, pwindow), log_shift
      );
    } else {
      integral += gumbel_numeric_z_half(
        d, dist_id, params, pwindow, mu, beta, z_peak, z_high, log_shift
      );
    }
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
  * Integrates the delay CDF against the window density with an ODE solver in
  * a variable in which the integrand is smooth, on the log scale so a CDF
  * below 1e-300 is not lost. The solver tolerances are 1e-10 relative, as
  * tighter ones do not converge with the sensitivities of a delay CDF.
  *
  * @param d Delay
  * @param dist_id Distribution identifier
  * @param params Array of distribution parameters
  * @param pwindow Primary event window
  * @param mu Location of the truncated Gumbel
  * @param beta Scale of the truncated Gumbel
  *
  * @return Log of the CDF, above 0 by at most the solver error, `-inf` if it
  *   is 0
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
  * Used by primarycensored_cdf() for the truncated Gumbel, whose density can
  * be a spike that integration over the window would miss.
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
  * Never NaN and never above 0.
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
  * Compute the primary event censored log CDF for a truncated Gumbel primary
  * @ingroup truncated_gumbel_solutions
  *
  * Uses the series where gumbel_series_accepted() is 1 and the numerical
  * path otherwise. Only for check_for_gumbel() and
  * check_for_gumbel_params() of 1.
  *
  * @param d Delay
  * @param dist_id Distribution identifier, 18 (Normal)
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
  int n_terms = gumbel_n_terms(mu / beta);
  vector[2] fit = primarycensored_gumbel_lcdf_from_terms(
    primarycensored_gumbel_terms(d, dist_id, beta, n_terms, params),
    primarycensored_gumbel_terms(d - pwindow, dist_id, beta, n_terms, params),
    d, pwindow, mu, beta, n_terms
  );
  if (gumbel_series_accepted(fit)) {
    return fit[1];
  }
  return primarycensored_gumbel_numeric_lcdf(
    d | dist_id, params, pwindow, mu, beta
  );
}
