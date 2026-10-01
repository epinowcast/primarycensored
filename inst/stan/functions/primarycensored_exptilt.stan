/*
 * Exponentially tilted primary event window
 *
 * For a window of width w with density rho exp(rho z) / (exp(rho w) - 1) on
 * [0, w], with q = d - w, F the delay CDF and J(x) = T_f(-rho; x) the tilt
 * transform of tilt_transform.stan,
 *   F_rho(d) = F(q) + (exp(rho d) (J(d) - J(q)) - (F(d) - F(q))) /
 *              (exp(rho w) - 1).
 * Every term depends on d or q alone, so the terms are computed once per
 * endpoint in primarycensored_analytical_lcdf_vectorized().
 *
 * The direct form cancels as rho goes to zero, see ?pcens_cdf_exptilt. For
 * |rho| * pwindow below 1e-2 (1e-5 for delays on the reals) the small window
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
  if (abs(rho) * pwindow < (positive ? 1e-2 : 1e-5)) {
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
