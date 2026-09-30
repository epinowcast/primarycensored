/*
 * Truncated logistic primary event window
 *
 * The primary event censored CDF is F(a) + Phi / D_L with a = d - w and
 * D_L the mass of the window, where Phi is a sum of series of the transforms
 * of tilt_transform.stan at the tilts +/- n / s. The series are summed with
 * binomial taper weights, with the number of terms chosen by
 * tlogis_series_terms() to bound the truncation error. Each transform
 * depends on one endpoint, a, u or d, so the terms of
 * primarycensored_tlogis_terms() are computed once per endpoint and shared
 * in primarycensored_tlogis_lcdf_vectorized(). The series are given in the
 * R documentation of pcens_cdf_tlogis.
 */

/**
  * Check if the analytical solution is the truncated logistic solution
  * @ingroup truncated_logistic_solutions
  *
  * The delay needs a tilt transform, see check_for_tilt_transform(), and the
  * truncated logistic primary is primary_id 3 with
  * primary_params = [location, scale].
  *
  * @param dist_id Distribution identifier for the delay distribution
  * @param primary_id Distribution identifier for the primary distribution
  *
  * @return 1 if the delay has a truncated logistic solution and the primary
  * is truncated logistic, 0 otherwise. Whether it applies for given
  * parameters is check_for_tlogis_params().
  */
int check_for_tlogis(int dist_id, int primary_id) {
  return primary_id == 3 && (dist_id == 2 || dist_id == 4 || dist_id == 18);
}

/**
  * Log of the largest truncation error of a weighted geometric series
  * @ingroup truncated_logistic_solutions
  *
  * The series sum (-x)^n is summed with weights that are 1 for n < n0 and
  * taper to 0 at n0 + M, see tlogis_weights(). For one x the error is
  * x^n0 ((1 - x) / 2)^M / (1 + x). Its largest value for x in
  * [exp(log_lo), exp(log_hi)], without the last factor, is at
  * n0 / (n0 + M) clipped to the range, and at the lower end for n0 = 0.
  *
  * @param n0 Number of terms with weight 1
  * @param M Number of terms that taper to 0
  * @param log_lo Log of the smallest x
  * @param log_hi Log of the largest x, at most 0
  *
  * @return Log of the bound of the error
  */
real tlogis_log_bound(int n0, int M, real log_lo, real log_hi) {
  if (M == 0) return n0 * log_hi;
  if (n0 == 0) return M * (log1m_exp(log_lo) - log2());
  real log_x = fmin(fmax(log(1.0 * n0 / (n0 + M)), log_lo), log_hi);
  return M * (log1m_exp(log_x) - log2()) + n0 * log_x;
}

/**
  * Number of terms and taper of a series of the truncated logistic window
  * @ingroup truncated_logistic_solutions
  *
  * Searches the plain partial sum, a taper after 35% of the terms and a taper
  * from the first term for the fewest terms, at most 64, with the bound of
  * tlogis_log_bound() below the tolerance. The same rule as R.
  *
  * @param log_lo Log of the smallest x of the series
  * @param log_hi Log of the largest x of the series
  * @param log_tol Log of the tolerance on the bound of the error
  *
  * @return Array {n0, M}, or {-1, -1} if no rule meets the tolerance
  */
array[] int tlogis_series_terms(real log_lo, real log_hi, real log_tol) {
  if (is_inf(log_tol) || is_nan(log_tol)) return {-1, -1};
  for (k in 1:64) {
    int taper = to_int(floor(0.35 * k));
    if (tlogis_log_bound(k, 0, log_lo, log_hi) <= log_tol) return {k, 0};
    if (tlogis_log_bound(taper, k - taper, log_lo, log_hi) <= log_tol) {
      return {taper, k - taper};
    }
    if (tlogis_log_bound(0, k, log_lo, log_hi) <= log_tol) return {0, k};
  }
  return {-1, -1};
}

/**
  * Weights of the terms of a series of the truncated logistic window
  * @ingroup truncated_logistic_solutions
  *
  * The weights of the average of the partial sums S_n0, ..., S_{n0 + M} with
  * Binomial(M, 1/2) weights. They are 1 for n < n0, and the probability that
  * the Binomial is more than n - n0 for the rest.
  *
  * @param n0 Number of terms with weight 1
  * @param M Number of terms that taper to 0
  *
  * @return Vector of the n0 + M weights
  */
vector tlogis_weights(int n0, int M) {
  vector[n0 + M] weights = rep_vector(1, n0 + M);
  real tail = 0;
  // The upper tail from the top, so that small weights keep their precision
  for (i in 1:M) {
    int top = M - i + 1;
    tail += exp(lchoose(M, top) - M * log2());
    weights[n0 + top] = tail;
  }
  return weights;
}

/**
  * Plan the series of the truncated logistic window
  * @ingroup truncated_logistic_solutions
  *
  * Finds the series that the CDF needs for the location and the scale, with
  * n0 and M of each, see tlogis_series_terms(). A location before the window
  * needs one series of positive tilts, n >= 1, over the whole window. A
  * location inside the window needs a series of positive tilts n >= 0 over
  * the part above it, and of negative tilts n >= 1 over the part below. A
  * location after the window needs only negative tilts. The tolerance is
  * relative to the mass of the window.
  *
  * @param location Location of the logistic distribution
  * @param scale Scale of the logistic distribution
  * @param pwindow Primary event window
  * @param log_mass Log of the mass of the window, see tlogis_log_mass()
  *
  * @return Array {ok, pos_form, pos_n0, pos_M, neg_n0, neg_M} where ok is 0
  * if a series needed has no rule and then the rest are 0, pos_form is 0 for
  * no series of positive tilts, 1 for n >= 1 over the whole window and 2 for
  * n >= 0 over the part above the location, and neg_n0 + neg_M is 0 for no
  * series of negative tilts
  */
array[] int tlogis_plan(real location, real scale, data real pwindow,
                        real log_mass) {
  real log_tol = log(1e-10) + log_mass;
  array[2] int pos = {0, 0};
  array[2] int neg = {0, 0};
  int pos_form = 0;
  if (location < 0) {
    pos_form = 1;
    pos = tlogis_series_terms(
      (location - pwindow) / scale, location / scale, log_tol
    );
    if (pos[1] < 0) return {0, 0, 0, 0, 0, 0};
  } else {
    if (location < pwindow) {
      pos_form = 2;
      pos = tlogis_series_terms(-(pwindow - location) / scale, 0, log_tol);
      if (pos[1] < 0) return {0, 0, 0, 0, 0, 0};
    }
    if (location > 0) {
      neg = tlogis_series_terms(
        -location / scale, fmin(0, (pwindow - location) / scale), log_tol
      );
      if (neg[1] < 0) return {0, 0, 0, 0, 0, 0};
    }
  }
  return {1, pos_form, pos[1], pos[2], neg[1], neg[2]};
}

/**
  * Number of positive tilts of a plan beyond the tilt 0
  * @ingroup truncated_logistic_solutions
  *
  * @param plan Plan from tlogis_plan()
  *
  * @return The positive tilts n / s of the series with n >= 1, which is all
  * of the terms of the series over the whole window and all but the first
  * of the series over the part above the location
  */
int tlogis_n_pos(array[] int plan) {
  if (plan[2] == 0) return 0;
  return plan[3] + plan[4] - (plan[2] == 2 ? 1 : 0);
}

/**
  * Number of negative tilts of a plan
  * @ingroup truncated_logistic_solutions
  *
  * @param plan Plan from tlogis_plan()
  *
  * @return The negative tilts -n / s for n >= 1 of the plan
  */
int tlogis_n_neg(array[] int plan) {
  return plan[5] + plan[6];
}

/**
  * Check if the truncated logistic solution applies for the parameters
  * @ingroup truncated_logistic_solutions
  *
  * The window needs to be a positive finite number, the scale positive and
  * the location finite. A truncation rule needs to exist for each series,
  * see tlogis_plan(), and the delay needs the largest positive and the
  * largest negative tilt of the series, see check_for_tilt_transform(). The
  * exponential and gamma delays need the largest positive tilt to be below
  * their rate, which fails unless the location is after the window or the
  * rate is large. Where this is 0 the numerical path is used.
  *
  * @param dist_id Distribution identifier for the delay distribution
  * @param params Array of delay distribution parameters
  * @param primary_params Array [location, scale] of the primary
  * @param pwindow Primary event window
  *
  * @return 1 if the truncated logistic solution applies, 0 otherwise
  */
int check_for_tlogis_params(int dist_id, array[] real params,
                            array[] real primary_params,
                            data real pwindow) {
  real location = primary_params[1];
  real scale = primary_params[2];
  if (!(pwindow > 0) || is_inf(pwindow)) return 0;
  if (!(scale > 0) || is_inf(scale) || is_nan(location) || is_inf(location)) {
    return 0;
  }
  real log_mass = tlogis_log_mass(0, pwindow, location, scale);
  array[6] int plan = tlogis_plan(location, scale, pwindow, log_mass);
  if (plan[1] == 0) return 0;
  // The series run to the tilts n / s for n up to the number of positive
  // tilts, and -n / s for n up to the number of negative tilts
  if (tlogis_n_pos(plan) > 0) {
    if (!check_for_tilt_transform(
      dist_id, tlogis_n_pos(plan) / scale, params
    )) {
      return 0;
    }
  }
  if (tlogis_n_neg(plan) > 0) {
    if (!check_for_tilt_transform(
      dist_id, -tlogis_n_neg(plan) / scale, params
    )) {
      return 0;
    }
  }
  return 1;
}

/**
  * Compute the truncated logistic terms at an endpoint
  * @ingroup truncated_logistic_solutions
  *
  * @param t Endpoint, d, d - pwindow or the split point
  * @param dist_id Distribution identifier, see check_for_tlogis()
  * @param params Array of delay distribution parameters
  * @param scale Scale of the logistic distribution
  * @param n_pos Number of positive tilts n / s, n = 1, ..., n_pos
  * @param n_neg Number of negative tilts -n / s, n = 1, ..., n_neg
  * @param use_pos 1 to compute the positive tilts, 0 to leave them 0
  * @param use_neg 1 to compute the negative tilts, 0 to leave them 0
  *
  * @return A matrix with 2 rows, the log of the transform over the lower and
  * the upper part of the support, and 1 + n_pos + n_neg columns, for the tilt
  * 0 (the delay CDF), then the positive tilts and then the negative tilts.
  * The columns that are not computed are 0. Only defined where
  * check_for_tlogis_params() is 1.
  */
matrix primarycensored_tlogis_terms(real t, int dist_id, array[] real params,
                                    real scale, int n_pos, int n_neg,
                                    int use_pos, int use_neg) {
  matrix[2, 1 + n_pos + n_neg] terms = rep_matrix(0, 2, 1 + n_pos + n_neg);
  terms[, 1] = log_tilt_transform_pair(t, dist_id, 0, params);
  if (use_pos) {
    for (n in 1:n_pos) {
      terms[, 1 + n] = log_tilt_transform_pair(t, dist_id, n / scale, params);
    }
  }
  if (use_neg) {
    for (n in 1:n_neg) {
      terms[, 1 + n_pos + n] = log_tilt_transform_pair(
        t, dist_id, -n / scale, params
      );
    }
  }
  return terms;
}

/**
  * Log of the difference of a tilt transform between two endpoints
  * @ingroup truncated_logistic_solutions
  *
  * @param terms_hi Terms at the upper endpoint
  * @param terms_lo Terms at the lower endpoint
  * @param column Column of the tilt, 1 for the delay CDF
  *
  * @return Log of T_f(xi; hi) - T_f(xi; lo), from the lower or the upper
  * transforms, see primarycensored_tail_diff()
  */
real tlogis_terms_diff(matrix terms_hi, matrix terms_lo, int column) {
  return primarycensored_tail_diff(
    terms_hi[1, column], terms_lo[1, column],
    terms_hi[2, column], terms_lo[2, column]
  );
}

/**
  * Weighted alternating sum of a series with its scale
  * @ingroup truncated_logistic_solutions
  *
  * @param log_c Log of the terms c_n
  * @param weights Weights of the terms, see tlogis_weights()
  * @param log_scale Log scale to divide the terms by
  *
  * @return The sum of (-1)^n W_n c_n divided by exp(log_scale)
  */
real tlogis_weighted_sum(vector log_c, vector weights, real log_scale) {
  int n = rows(weights);
  real total = 0;
  for (i in 1:n) {
    total += (i % 2 == 1 ? 1 : -1) * weights[i] * exp(log_c[i] - log_scale);
  }
  return total;
}

/**
  * Part of the truncated logistic integral for a location before the window
  * @ingroup truncated_logistic_solutions
  *
  * For m < 0 the whole window is above the location and
  * Phi = sum_{n >= 1} (-1)^{n - 1} exp(n m / s)
  *   {dF(a, d) - exp(-n d / s) dT(n / s; a, d)}, with the terms weighted.
  *
  * @param terms_a Terms at a = d - pwindow
  * @param terms_d Terms at d
  * @param d Delay
  * @param location Location of the logistic distribution
  * @param scale Scale of the logistic distribution
  * @param weights Weights of the terms, see tlogis_weights()
  *
  * @return Vector [log scale, sum], the part is exp(log scale) times the sum,
  * and both are 0 where the terms are all 0
  */
vector tlogis_part_before_window(matrix terms_a, matrix terms_d,
                                 data real d, real location, real scale,
                                 vector weights) {
  int n = rows(weights);
  real log_df = tlogis_terms_diff(terms_d, terms_a, 1);
  vector[n] log_c;
  for (k in 1:n) {
    log_c[k] = k * location / scale + primarycensored_log_diff_exp(
      log_df, -k * d / scale + tlogis_terms_diff(terms_d, terms_a, 1 + k)
    );
  }
  real log_scale = max(log_c);
  if (is_inf(log_scale)) return [0, 0]';
  return [log_scale, tlogis_weighted_sum(log_c, weights, log_scale)]';
}

/**
  * Part of the truncated logistic integral above the location
  * @ingroup truncated_logistic_solutions
  *
  * For 0 <= m < w the part of the window with the primary event time above
  * the location, from a to u = d - m, is
  * Phi_P = sum_{n >= 0} (-1)^n exp(-n (d - m) / s) dT(n / s; a, u)
  *   - L(0) dF(a, u), with the terms weighted.
  *
  * @param terms_a Terms at a = d - pwindow
  * @param terms_u Terms at u = d - location
  * @param d Delay
  * @param location Location of the logistic distribution
  * @param scale Scale of the logistic distribution
  * @param weights Weights of the terms, see tlogis_weights()
  *
  * @return Vector [log scale, sum], see tlogis_part_before_window()
  */
vector tlogis_part_above_location(matrix terms_a, matrix terms_u,
                                  data real d, real location, real scale,
                                  vector weights) {
  int n = rows(weights);
  real log_df = tlogis_terms_diff(terms_u, terms_a, 1);
  vector[n] log_c;
  for (k in 0:(n - 1)) {
    log_c[k + 1] = (k == 0 ? log_df
                           : tlogis_terms_diff(terms_u, terms_a, 1 + k))
                   - k * (d - location) / scale;
  }
  // The constant term of the series, L(0) times the mass of the delay
  real log_constant = log_inv_logit(-location / scale) + log_df;
  real log_scale = fmax(max(log_c), log_constant);
  if (is_inf(log_scale)) return [0, 0]';
  return [
    log_scale,
    tlogis_weighted_sum(log_c, weights, log_scale)
    - exp(log_constant - log_scale)
  ]';
}

/**
  * Part of the truncated logistic integral below the location
  * @ingroup truncated_logistic_solutions
  *
  * For m > 0 the part of the window with the primary event time below the
  * location, from u = d - min(m, w) to d, is
  * Phi_N = sum_{n >= 1} (-1)^{n - 1} exp(-n m / s)
  *   {exp(n d / s) dT(-n / s; u, d) - dF(u, d)}, with the terms weighted.
  *
  * @param terms_u Terms at u = d - min(location, pwindow)
  * @param terms_d Terms at d
  * @param n_pos Number of positive tilt columns of the terms
  * @param d Delay
  * @param location Location of the logistic distribution
  * @param scale Scale of the logistic distribution
  * @param weights Weights of the terms, see tlogis_weights()
  *
  * @return Vector [log scale, sum], see tlogis_part_before_window()
  */
vector tlogis_part_below_location(matrix terms_u, matrix terms_d, int n_pos,
                                  data real d, real location, real scale,
                                  vector weights) {
  int n = rows(weights);
  real log_df = tlogis_terms_diff(terms_d, terms_u, 1);
  vector[n] log_c;
  for (k in 1:n) {
    log_c[k] = -k * location / scale + primarycensored_log_diff_exp(
      k * d / scale + tlogis_terms_diff(terms_d, terms_u, 1 + n_pos + k),
      log_df
    );
  }
  real log_scale = max(log_c);
  if (is_inf(log_scale)) return [0, 0]';
  return [log_scale, tlogis_weighted_sum(log_c, weights, log_scale)]';
}

/**
  * Combine the terms at a, u and d into the truncated logistic log CDF
  * @ingroup truncated_logistic_solutions
  *
  * The direct form. Each series gives a part of Phi on its own scale, and
  * F_L(d) = F(a) + Phi / D_L. Use the small delay form where
  * tlogis_is_small_delay() is 1.
  *
  * @param terms_a Terms at a = d - pwindow
  * @param terms_u Terms at u = d - min(max(location, 0), pwindow), which are
  *   those at d for a location at most 0 and at a for a location at least
  *   pwindow
  * @param terms_d Terms at d
  * @param d Delay
  * @param location Location of the logistic distribution
  * @param scale Scale of the logistic distribution
  * @param log_mass Log of the mass of the window
  * @param plan Plan from tlogis_plan()
  * @param w_pos Weights of the series of positive tilts
  * @param w_neg Weights of the series of negative tilts
  *
  * @return Log of the primary event censored CDF at d
  */
real primarycensored_tlogis_lcdf_from_terms(
  matrix terms_a, matrix terms_u, matrix terms_d, data real d, real location,
  real scale, real log_mass, array[] int plan, vector w_pos, vector w_neg
) {
  real log_phi = negative_infinity();
  if (plan[2] == 1) {
    vector[2] part = tlogis_part_before_window(
      terms_a, terms_d, d, location, scale, w_pos
    );
    if (part[2] > 0) log_phi = part[1] + log(part[2]);
  } else if (plan[2] == 2) {
    vector[2] part = tlogis_part_above_location(
      terms_a, terms_u, d, location, scale, w_pos
    );
    if (part[2] > 0) log_phi = part[1] + log(part[2]);
  }
  if (tlogis_n_neg(plan) > 0) {
    vector[2] part = tlogis_part_below_location(
      terms_u, terms_d, tlogis_n_pos(plan), d, location, scale, w_neg
    );
    if (part[2] > 0) log_phi = log_sum_exp(log_phi, part[1] + log(part[2]));
  }
  real log_f_a = terms_a[1, 1];
  // Deep enough into the lower tail every term underflows together. Both
  // are then `-inf` and log_sum_exp would differentiate to NaN.
  if (log_f_a == negative_infinity() && log_phi == negative_infinity()) {
    return negative_infinity();
  }
  // The CDF is at most 1. Rounding can put the log CDF a little above 0.
  return fmin(log_sum_exp(log_f_a, log_phi - log_mass), 0);
}

/**
  * Check if the small delay form replaces the direct form
  * @ingroup truncated_logistic_solutions
  *
  * For delays on the non-negative reals with d below the window and d / scale
  * below 1e-4 only the primary event times up to d contribute and the direct
  * form cancels. Where this is 1 use
  * primarycensored_tlogis_small_delay_lcdf().
  *
  * @param dist_id Distribution identifier
  * @param d Delay
  * @param pwindow Primary event window
  * @param scale Scale of the logistic distribution
  *
  * @return 1 if the delay has non-negative support, 0 < d < pwindow and
  * d / scale is below 1e-4, 0 otherwise
  */
int tlogis_is_small_delay(int dist_id, data real d, data real pwindow,
                          real scale) {
  return dist_has_positive_support(dist_id) && d > 0 && d < pwindow
         && d < 1e-4 * scale;
}

/**
  * Small delay form of the truncated logistic log CDF
  * @ingroup truncated_logistic_solutions
  *
  * With G_k(d) = int_0^d (d - u)^k f(u) du,
  * F_L(d) = {L'(0) G_1(d) + L''(0) G_2(d) / 2} / D_L + O((d / s)^2), with
  * L'(0) = L(0) (1 - L(0)) / s and L''(0) = L'(0) (1 - 2 L(0)) / s.
  *
  * @param d Delay, positive
  * @param dist_id Distribution identifier, see check_for_tlogis()
  * @param params Array of delay distribution parameters
  * @param location Location of the logistic distribution
  * @param scale Scale of the logistic distribution
  * @param log_mass Log of the mass of the window
  *
  * @return Log of the primary event censored CDF at d
  */
real primarycensored_tlogis_small_delay_lcdf(data real d, int dist_id,
                                             array[] real params,
                                             real location, real scale,
                                             real log_mass) {
  vector[2] moments = primarycensored_tilt_moments(d, dist_id, params);
  if (moments[1] == negative_infinity()) return negative_infinity();
  real log_l1 = log_inv_logit(-location / scale)
                + log_inv_logit(location / scale) - log(scale);
  real ratio = (1 - 2 * inv_logit(-location / scale)) / scale;
  return fmin(
    moments[1] + log_l1 - log_mass
    + log1p(0.5 * ratio * exp(moments[2] - moments[1])),
    0
  );
}

/**
  * Compute the primary event censored log CDF for a truncated logistic
  * primary
  * @ingroup truncated_logistic_solutions
  *
  * Chooses the direct form or the small delay form, see
  * tlogis_is_small_delay(). Only for check_for_tlogis() is 1 and
  * check_for_tlogis_params() is 1.
  *
  * The gamma shape gradient comes from the tilt transforms, see
  * log_tilt_transform_pair(), not from Stan's `gamma_lccdf()` and
  * `gamma_lcdf()`, whose shape gradients are too inaccurate for the series.
  *
  * @param d Delay
  * @param dist_id Distribution identifier: 2 (Gamma), 4 (Exponential) or 18
  *   (Normal)
  * @param params Array of delay distribution parameters
  * @param pwindow Primary event window
  * @param location Location of the logistic primary
  * @param scale Scale of the logistic primary
  *
  * @return Log of the primary event censored CDF at d
  */
real primarycensored_tlogis_lcdf(data real d, int dist_id,
                                 array[] real params, data real pwindow,
                                 real location, real scale) {
  if (dist_has_positive_support(dist_id) && d <= 0) {
    return negative_infinity();
  }
  real log_mass = tlogis_log_mass(0, pwindow, location, scale);
  if (tlogis_is_small_delay(dist_id, d, pwindow, scale)) {
    return primarycensored_tlogis_small_delay_lcdf(
      d | dist_id, params, location, scale, log_mass
    );
  }
  array[6] int plan = tlogis_plan(location, scale, pwindow, log_mass);
  int n_pos = tlogis_n_pos(plan);
  int n_neg = tlogis_n_neg(plan);
  int has_pos = plan[2] > 0;
  int has_neg = n_neg > 0;
  vector[plan[3] + plan[4]] w_pos = tlogis_weights(plan[3], plan[4]);
  vector[plan[5] + plan[6]] w_neg = tlogis_weights(plan[5], plan[6]);
  real a = d - pwindow;
  // The split point is d for a location at most 0 and a at least pwindow
  int u_is_d = location <= 0;
  int u_is_a = location >= pwindow;
  matrix[2, 1 + n_pos + n_neg] terms_a = primarycensored_tlogis_terms(
    a, dist_id, params, scale, n_pos, n_neg, has_pos, has_neg && u_is_a
  );
  matrix[2, 1 + n_pos + n_neg] terms_d = primarycensored_tlogis_terms(
    d, dist_id, params, scale, n_pos, n_neg, has_pos && u_is_d, has_neg
  );
  if (u_is_d) {
    return primarycensored_tlogis_lcdf_from_terms(
      terms_a, terms_d, terms_d, d, location, scale, log_mass, plan, w_pos,
      w_neg
    );
  } else if (u_is_a) {
    return primarycensored_tlogis_lcdf_from_terms(
      terms_a, terms_a, terms_d, d, location, scale, log_mass, plan, w_pos,
      w_neg
    );
  }
  return primarycensored_tlogis_lcdf_from_terms(
    terms_a,
    primarycensored_tlogis_terms(
      d - location, dist_id, params, scale, n_pos, n_neg, has_pos, has_neg
    ),
    terms_d, d, location, scale, log_mass, plan, w_pos, w_neg
  );
}

/**
  * Check if the truncated logistic solution can be vectorised over integer
  * delays
  * @ingroup truncated_logistic_solutions
  *
  * With an integer pwindow, a = d - pwindow is an integer delay too, so
  * primarycensored_tlogis_lcdf_vectorized() can compute the terms once per
  * delay and share them.
  *
  * @param dist_id Distribution identifier for the delay distribution
  * @param primary_id Distribution identifier for the primary distribution
  * @param pwindow Primary event window
  *
  * @return 1 if the vectorised truncated logistic solution applies, 0
  * otherwise
  */
int check_for_tlogis_vectorized(int dist_id, int primary_id,
                                data real pwindow) {
  return check_for_tlogis(dist_id, primary_id)
         && pwindow >= 1 && floor(pwindow) == pwindow;
}

/**
  * Compute the truncated logistic primary event censored log CDF at integer
  * delays
  * @ingroup truncated_logistic_solutions
  *
  * The log CDF at d combines the terms at d, at a = d - pwindow and, for a
  * location inside the window, at the split point u = d - location. The
  * first two are integer delays, so the terms are computed once per delay and
  * shared, halving the transform evaluations. The split points are integer
  * delays only for an integer location, where they are shared too. For a
  * location inside the window that is not an integer each split point is
  * different, so the delays use primarycensored_tlogis_lcdf() one at a time.
  * The values are the same as from primarycensored_tlogis_lcdf() at each
  * delay. Only for cases where check_for_tlogis_vectorized() is 1 and
  * check_for_tlogis_params() is 1.
  *
  * @param start First delay to compute
  * @param n Last delay to compute, and the length of the result
  * @param dist_id Distribution identifier: 2 (Gamma), 4 (Exponential) or 18
  *   (Normal)
  * @param params Array of delay distribution parameters
  * @param pwindow Primary event window, a positive integer
  * @param location Location of the logistic primary
  * @param scale Scale of the logistic primary
  *
  * @return Vector whose element d is the log CDF at d, for d in start:n.
  * Elements before start are not computed.
  */
vector primarycensored_tlogis_lcdf_vectorized(data int start, data int n,
                                              data int dist_id,
                                              array[] real params,
                                              data real pwindow,
                                              real location, real scale) {
  int pw = to_int(pwindow);
  int positive = dist_has_positive_support(dist_id);
  vector[n] log_cdfs;
  real log_mass = tlogis_log_mass(0, pwindow, location, scale);
  array[6] int plan = tlogis_plan(location, scale, pwindow, log_mass);
  int n_pos = tlogis_n_pos(plan);
  int n_neg = tlogis_n_neg(plan);
  int has_pos = plan[2] > 0;
  int has_neg = n_neg > 0;
  // The integer split point, for a location that is an integer or outside
  // the window
  int sp = -1;
  if (location <= 0) {
    sp = 0;
  } else if (location >= pwindow) {
    sp = pw;
  } else {
    for (i in 1:(pw - 1)) {
      if (location == i) sp = i;
    }
  }
  if (sp < 0) {
    for (d in start:n) {
      log_cdfs[d] = primarycensored_tlogis_lcdf(
        d | dist_id, params, pwindow, location, scale
      );
    }
    return log_cdfs;
  }
  vector[plan[3] + plan[4]] w_pos = tlogis_weights(plan[3], plan[4]);
  vector[plan[5] + plan[6]] w_neg = tlogis_weights(plan[5], plan[6]);
  // Endpoints below 0 have the same terms as 0 for delays on the
  // non-negative reals, so they share the entry for 0
  int first = positive ? max(start - pw, 0) : start - pw;
  // terms[t - first + 1] holds the terms at endpoint t
  array[n - first + 1] matrix[2, 1 + n_pos + n_neg] terms;
  for (t in first:n) {
    terms[t - first + 1] = primarycensored_tlogis_terms(
      t, dist_id, params, scale, n_pos, n_neg, has_pos, has_neg
    );
  }
  for (d in start:n) {
    if (tlogis_is_small_delay(dist_id, d, pwindow, scale)) {
      log_cdfs[d] = primarycensored_tlogis_small_delay_lcdf(
        d | dist_id, params, location, scale, log_mass
      );
    } else {
      int a_index = (positive ? max(d - pw, 0) : d - pw) - first + 1;
      int u_index = (positive ? max(d - sp, 0) : d - sp) - first + 1;
      log_cdfs[d] = primarycensored_tlogis_lcdf_from_terms(
        terms[a_index], terms[u_index], terms[d - first + 1], d, location,
        scale, log_mass, plan, w_pos, w_neg
      );
    }
  }
  return log_cdfs;
}
