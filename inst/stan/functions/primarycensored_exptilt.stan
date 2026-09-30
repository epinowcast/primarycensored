/*
 * Exponentially tilted primary event window
 *
 * For a primary window of width w with density rho exp(rho z) / (exp(rho w) -
 * 1) on [0, w], the primary event censored CDF at delay d is
 *   F_rho(d) = F(q) + (exp(rho d) (J(d) - J(q)) - (F(d) - F(q))) /
 *              (exp(rho w) - 1),
 * with q = d - w, F the delay CDF and J(x) = T_f(-rho; x) the tilt
 * transform of the delay, see tilt_transform.stan. F(x) = J(x) = 0 for
 * x <= 0 for delays on the non-negative reals. Every term depends on d or q
 * alone, so the terms of primarycensored_exptilt_terms() are computed once
 * per endpoint and shared between neighbouring delays in
 * primarycensored_exptilt_lcdf_vectorized().
 *
 * A window of another shape plugs in the same way. It needs the transform at
 * its own tilts and its own combination of the terms. A delay distribution
 * plugs in through tilt_transform.stan and is then supported by every
 * window.
 */

/**
  * Check if the analytical solution is the exponentially tilted solution
  * @ingroup exponential_tilt_solutions
  *
  * The delay needs a tilt transform, see check_for_tilt_transform(), and the
  * tilted primary is primary_id 2 with primary_params = [rho].
  *
  * @param dist_id Distribution identifier for the delay distribution
  * @param primary_id Distribution identifier for the primary distribution
  *
  * @return 1 if the delay has an exponentially tilted solution and the
  * primary is exponentially tilted, 0 otherwise. Whether it applies for
  * given parameters is check_for_analytical_params().
  */
int check_for_exptilt(int dist_id, int primary_id) {
  return primary_id == 2 && (dist_id == 2 || dist_id == 4 || dist_id == 18);
}

/**
  * Check if the small tilt forms replace the direct form
  * @ingroup exponential_tilt_solutions
  *
  * The direct form cancels as rho goes to zero and loses about
  * 1e-14 / (|rho| w) relative precision. The form for small |rho| w
  * (primarycensored_exptilt_small_window_lcdf_from_terms()) has a truncation
  * error of about (|rho| w)^2 / 12. They are both below 1e-9 at the
  * threshold 1e-4. The derivative of the small form in rho is less accurate,
  * with a relative error of about |rho| w / 6, up to 2e-5 at the threshold,
  * as the second order term is not included. Its derivatives in the delay
  * parameters have the error of its value.
  *
  * @param rho Tilt
  * @param pwindow Primary event window
  *
  * @return 1 if |rho| * pwindow is below 1e-4, 0 otherwise
  */
int exptilt_is_small_window(real rho, data real pwindow) {
  return abs(rho) * pwindow < 1e-4;
}

/**
  * Check if the small delay form replaces the direct form
  * @ingroup exponential_tilt_solutions
  *
  * For delays on the non-negative reals with d < pwindow the terms at
  * d - pwindow are zero. The direct form then cancels as |rho| d goes to
  * zero, and primarycensored_exptilt_small_delay_lcdf_from_terms() keeps
  * precision for d close to zero.
  *
  * @param dist_id Distribution identifier
  * @param rho Tilt
  * @param d Delay
  * @param pwindow Primary event window
  *
  * @return 1 if the delay has non-negative support, d < pwindow and
  * |rho| * d is below 1e-4, 0 otherwise
  */
int exptilt_is_small_delay(int dist_id, real rho, data real d,
                           data real pwindow) {
  return dist_has_positive_support(dist_id) && d < pwindow
         && abs(rho) * d < 1e-4;
}

/**
  * Compute the exponentially tilted terms at an endpoint
  * @ingroup exponential_tilt_solutions
  *
  * @param t Endpoint, d or d - pwindow
  * @param dist_id Distribution identifier, see check_for_exptilt()
  * @param rho Tilt
  * @param params Array of distribution parameters
  *
  * @return Vector [log F(t), log(1 - F(t)), log J(t), log(J(Inf) - J(t))].
  * The lower tail terms are `-inf` for t <= 0 for delays on the non-negative
  * reals. Only defined where check_for_tilt_transform() is 1 for -rho.
  */
vector primarycensored_exptilt_terms(real t, int dist_id, real rho,
                                     array[] real params) {
  return append_row(
    log_tilt_transform_pair(t, dist_id, 0, params),
    log_tilt_transform_pair(t, dist_id, -rho, params)
  );
}

/**
  * Log of a difference between two points from either tail
  * @ingroup exponential_tilt_solutions
  *
  * For a lower tail quantity L and an upper tail quantity U with L + U
  * constant, L(d) - L(q) equals U(q) - U(d). The relative precision of a
  * difference is best when the subtracted term is small, so this takes
  * the representation with the smaller ratio of the terms. This avoids
  * cancellation in the upper tail, where for positive tilts
  * exp(rho d) J(d) is much larger than the result.
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
  // Lower tail terms that both underflow give NaN, taken as a zero
  // difference
  if (is_nan(lower_q - lower_d) || lower_q - lower_d <= upper_d - upper_q) {
    return primarycensored_log_diff_exp(lower_d, lower_q);
  }
  return primarycensored_log_diff_exp(upper_q, upper_d);
}

/**
  * Combine the exponentially tilted terms at d and q into the log CDF
  * @ingroup exponential_tilt_solutions
  *
  * The direct form. Each of F(d) - F(q) and J(d) - J(q) is taken between the
  * lower tail terms or between the upper tail terms, whichever loses less
  * precision, see primarycensored_tail_diff(). For rho > 0 the numerator
  * exp(rho d) (J(d) - J(q)) - (F(d) - F(q)) is positive and for rho < 0 it
  * is negative, as is exp(rho w) - 1. Use the small tilt forms where
  * exptilt_is_small_window() or exptilt_is_small_delay() is 1.
  *
  * @param terms_d Terms at d from primarycensored_exptilt_terms()
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
  // Deep enough into the lower tail every term underflows together. Both
  // are then `-inf` and log_sum_exp would differentiate to NaN.
  if (terms_q[1] == negative_infinity() && log_num == negative_infinity()) {
    return negative_infinity();
  }
  // The CDF is at most 1. Rounding can put the log CDF a little above 0 in
  // the upper tail, where log_diff_exp() of two such values is NaN.
  return fmin(log_sum_exp(terms_q[1], log_num - log_den), 0);
}

/**
  * Combine the moments at d and q into the small tilt log CDF
  * @ingroup exponential_tilt_solutions
  *
  * The uniform window limit with its first order correction in the tilt,
  * for exptilt_is_small_window() is 1. With G_k(t) = int (t - u)^k f(u) du,
  *   F_rho(d) = (G_1(d) - G_1(q)) / w
  *     + rho (G_2(d) - w G_1(d) - G_2(q) - w G_1(q)) / (2 w) + O((rho w)^2)
  * which is the uniform window solution at rho = 0. Terms are scaled by
  * G_1(d) so nothing underflows when the CDF is small. The derivative in rho
  * has a relative error of about |rho| w / 6 from the O((rho w)^2) term, see
  * exptilt_is_small_window().
  *
  * @param moments_d Moments [log G_1, log G_2] at d from
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
  real relative = (1 - g1_q)
                  + 0.5 * rho * (g2_d - pwindow - g2_q - pwindow * g1_q);
  if (relative <= 0) return negative_infinity();
  return fmin(scale + log(relative) - log(pwindow), 0);
}

/**
  * Compute the small delay log CDF from the moments at d
  * @ingroup exponential_tilt_solutions
  *
  * For exptilt_is_small_delay() is 1 the terms at d - pwindow are zero and
  * F_rho(d) = rho (G_1(d) + rho G_2(d) / 2) / (exp(rho w) - 1) + O((rho d)^2).
  * The derivative in rho has a relative error of at most about |rho| d / 6.
  *
  * @param moments_d Moments [log G_1, log G_2] at d from
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
    + log1p(0.5 * rho * exp(moments_d[2] - moments_d[1])),
    0
  );
}

/**
  * Compute the primary event censored log CDF for an exponentially tilted
  * primary
  * @ingroup exponential_tilt_solutions
  *
  * Chooses the direct form or the small tilt forms, see
  * exptilt_is_small_window() and exptilt_is_small_delay(). A zero width
  * window has no primary event uncertainty, and the log CDF is the delay
  * log CDF, which all forms divide by the width to reach. Only for
  * check_for_exptilt() is 1 and check_for_tilt_transform() is 1 for -rho.
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
  if (exptilt_is_small_window(rho, pwindow)) {
    return primarycensored_exptilt_small_window_lcdf_from_terms(
      primarycensored_tilt_moments(d, dist_id, params),
      primarycensored_tilt_moments(q, dist_id, params), rho, pwindow
    );
  }
  if (exptilt_is_small_delay(dist_id, rho, d, pwindow)) {
    return primarycensored_exptilt_small_delay_lcdf_from_terms(
      primarycensored_tilt_moments(d, dist_id, params), rho, pwindow
    );
  }
  return primarycensored_exptilt_lcdf_from_terms(
    primarycensored_exptilt_terms(d, dist_id, rho, params),
    primarycensored_exptilt_terms(q, dist_id, rho, params), d, rho, pwindow
  );
}

/**
  * Check if the exponentially tilted solution can be vectorised over integer
  * delays
  * @ingroup exponential_tilt_solutions
  *
  * With an integer pwindow, q = d - pwindow is an integer delay too, so
  * primarycensored_exptilt_lcdf_vectorized() can compute the terms once per
  * delay and share them.
  *
  * @param dist_id Distribution identifier for the delay distribution
  * @param primary_id Distribution identifier for the primary distribution
  * @param pwindow Primary event window
  *
  * @return 1 if the vectorised exponentially tilted solution applies, 0
  * otherwise
  */
int check_for_exptilt_vectorized(int dist_id, int primary_id,
                                 data real pwindow) {
  return check_for_exptilt(dist_id, primary_id)
         && pwindow >= 1 && floor(pwindow) == pwindow;
}

/**
  * Compute the exponentially tilted primary event censored log CDF at
  * integer delays
  * @ingroup exponential_tilt_solutions
  *
  * The log CDF at d combines the terms at d and at q = d - pwindow. Both are
  * integer delays, so the terms are computed once per delay and used for
  * both, halving the transform evaluations. The values are the same as from
  * primarycensored_exptilt_lcdf() at each delay. Only for cases where
  * check_for_exptilt_vectorized() is 1 and check_for_tilt_transform() is 1
  * for -rho.
  *
  * @param start First delay to compute
  * @param n Last delay to compute, and the length of the result
  * @param dist_id Distribution identifier: 2 (Gamma), 4 (Exponential) or 18
  *   (Normal)
  * @param params Array of distribution parameters
  * @param pwindow Primary event window, a positive integer
  * @param rho Tilt, the exponential growth rate of the primary
  *
  * @return Vector whose element d is the log CDF at d, for d in start:n.
  * Elements before start are not computed.
  */
vector primarycensored_exptilt_lcdf_vectorized(data int start, data int n,
                                               data int dist_id,
                                               array[] real params,
                                               data real pwindow,
                                               real rho) {
  int pw = to_int(pwindow);
  int positive = dist_has_positive_support(dist_id);
  // Endpoints below 0 have the same terms as 0 for delays on the
  // non-negative reals, so they share the entry for 0
  int first = positive ? max(start - pw, 0) : start - pw;
  vector[n] log_cdfs;
  if (exptilt_is_small_window(rho, pwindow)) {
    // moments[t - first + 1] holds the moments at endpoint t
    array[n - first + 1] vector[2] moments;
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
      terms[t - first + 1] = primarycensored_exptilt_terms(
        t, dist_id, rho, params
      );
    }
    for (d in start:n) {
      if (exptilt_is_small_delay(dist_id, rho, d, pwindow)) {
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
