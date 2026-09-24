// Test-only reference: primarycensored_sone_lpmf_vectorized() and the
// analytical uniform primary log CDFs as they were before the vectorised
// function shared terms between nodes (primarycensored 1.5.2). Used by
// test-stan-sone-vectorized-equivalence.R to check the current functions
// return the same values and gradients.

real ref_primarycensored_gamma_uniform_lcdf(data real d, real q, array[] real params, data real pwindow) {
  real shape = params[1];
  real rate = params[2];
  real log_window = log(pwindow);
  // log E where E = k * theta = shape / rate is the mean of the delay
  real log_E = log(shape) - log(rate);

  // F_T(d; k) and the recursion to F_T(d; k+1):
  // P(k+1, y) = P(k, y) - y^k e^{-y} / Gamma(k+1), with y = rate * d
  real log_F_T_d_k = gamma_lcdf(d | shape, rate);
  real gamma_kp1_pdf_log_d
    = shape * log(rate * d) - rate * d - lgamma(shape + 1);
  real log_F_T_d_kp1 = log_diff_exp(log_F_T_d_k, gamma_kp1_pdf_log_d);

  // q-dependent terms. Final algebra is unified; only a guard to avoid
  // log_diff_exp(-inf, -inf) and log(0) when q == 0 (q is data, so autodiff
  // is unaffected by this branch).
  real log_q_F_T_q;    // log(q * F_T(q; k))
  real log_E_tF_T_q;   // log(E * F_T(q; k+1))
  if (q > 0) {
    real log_F_T_q_k = gamma_lcdf(q | shape, rate);
    real gamma_kp1_pdf_log_q
      = shape * log(rate * q) - rate * q - lgamma(shape + 1);
    real log_F_T_q_kp1 = log_diff_exp(log_F_T_q_k, gamma_kp1_pdf_log_q);
    log_q_F_T_q = log(q) + log_F_T_q_k;
    log_E_tF_T_q = log_E + log_F_T_q_kp1;
  } else {
    log_q_F_T_q = negative_infinity();
    log_E_tF_T_q = negative_infinity();
  }

  // Unified form: F_{S+}(d) = (A - B) / w_P with A, B sums of positives:
  //   A = d * F_T(d; k)   + E * F_T(q; k+1)
  //   B = q * F_T(q; k)   + E * F_T(d; k+1)
  // Ordering A >= B is guaranteed by F_{S+}(d) >= 0.
  real log_A = log_sum_exp(log(d) + log_F_T_d_k, log_E_tF_T_q);
  real log_B = log_sum_exp(log_q_F_T_q, log_E + log_F_T_d_kp1);

  return log_diff_exp(log_A, log_B) - log_window;
}
real ref_primarycensored_lognormal_uniform_lcdf(data real d, real q, array[] real params, data real pwindow) {
  real mu = params[1];
  real sigma = params[2];
  real mu_sigma2 = mu + square(sigma);
  real log_window = log(pwindow);
  // log E where E = exp(mu + sigma^2/2) is the mean of the delay
  real log_E = mu + 0.5 * square(sigma);

  // Each term is formed whole and dropped whole. Adding a `-inf` log CDF to
  // the parameter-dependent `log(d)` or `log_E` first would leave an edge
  // back to the parameters that `log_sum_exp` differentiates to
  // `exp(-inf - -inf)`. `q <= 0` underflows on the same test, so it needs no
  // separate branch.
  real log_d_F_T_d = lognormal_lcdf_underflows(d, mu, sigma)
                     ? negative_infinity()
                     : log(d) + lognormal_lcdf(d | mu, sigma);
  real log_E_tF_T_d = lognormal_lcdf_underflows(d, mu_sigma2, sigma)
                      ? negative_infinity()
                      : log_E + lognormal_lcdf(d | mu_sigma2, sigma);
  real log_q_F_T_q = lognormal_lcdf_underflows(q, mu, sigma)
                     ? negative_infinity()
                     : log(q) + lognormal_lcdf(q | mu, sigma);
  real log_E_tF_T_q = lognormal_lcdf_underflows(q, mu_sigma2, sigma)
                      ? negative_infinity()
                      : log_E + lognormal_lcdf(q | mu_sigma2, sigma);

  // Unified form: F_{S+}(d) = (A - B) / w_P with
  //   A = d * F_T(d) + E * tilde F_T(q)
  //   B = q * F_T(q) + E * tilde F_T(d)
  // Ordering A >= B is guaranteed by F_{S+}(d) >= 0.
  real log_A = log_sum_exp(log_d_F_T_d, log_E_tF_T_q);
  real log_B = log_sum_exp(log_q_F_T_q, log_E_tF_T_d);

  // Deep enough into the lower tail every term underflows together. Both
  // are then constant `-inf` and `log_diff_exp` would give NaN, so return
  // the limit directly.
  if (is_inf(log_A)) {
    return negative_infinity();
  }

  return log_diff_exp(log_A, log_B) - log_window;
}
real ref_primarycensored_weibull_uniform_lcdf(data real d, real q, array[] real params, data real pwindow) {
  real shape = params[1];
  real scale = params[2];
  real log_window = log(pwindow);
  real log_scale = log(scale);

  // For Weibull: E = scale (lambda) and tilde F_T(t) = g(t; lambda, k), so
  // log(E * tilde F_T(t)) = log(scale) + log_weibull_g(t, shape, scale).
  real log_F_T_d = weibull_lcdf(d | shape, scale);
  real log_E_tF_T_d = log_scale + log_weibull_g(d, shape, scale);

  // q-dependent terms (guard only to avoid log(0); final algebra is unified).
  real log_q_F_T_q;    // log(q * F_T(q))
  real log_E_tF_T_q;   // log(E * tilde F_T(q)) = log(scale * g(q; lambda, k))
  if (q > 0) {
    log_q_F_T_q = log(q) + weibull_lcdf(q | shape, scale);
    log_E_tF_T_q = log_scale + log_weibull_g(q, shape, scale);
  } else {
    log_q_F_T_q = negative_infinity();
    log_E_tF_T_q = negative_infinity();
  }

  // Unified form: F_{S+}(d) = (A - B) / w_P with
  //   A = d * F_T(d)    + scale * g(q; lambda, k)
  //   B = q * F_T(q)    + scale * g(d; lambda, k)
  // Ordering A >= B is guaranteed by F_{S+}(d) >= 0.
  real log_A = log_sum_exp(log(d) + log_F_T_d, log_E_tF_T_q);
  real log_B = log_sum_exp(log_q_F_T_q, log_E_tF_T_d);

  return log_diff_exp(log_A, log_B) - log_window;
}
real ref_primarycensored_gengamma_uniform_lcdf(data real d, real q, array[] real params, data real pwindow) {
  real shape = params[1];
  real scale = params[2];
  real k = params[3];
  real k_shift = k + inv(shape);
  real log_window = log(pwindow);
  // log E where E = scale * Gamma(k + 1/shape) / Gamma(k) is the mean of the
  // delay
  real log_E = log(scale) + lgamma(k_shift) - lgamma(k);

  real log_F_T_d = gengamma_lcdf(d | shape, scale, k);
  real log_tF_T_d = gengamma_lcdf(d | shape, scale, k_shift);

  // q-dependent terms (guard only to avoid log(0); final algebra is unified).
  real log_q_F_T_q;    // log(q * F_T(q))
  real log_E_tF_T_q;   // log(E * tilde F_T(q))
  if (q > 0) {
    log_q_F_T_q = log(q) + gengamma_lcdf(q | shape, scale, k);
    log_E_tF_T_q = log_E + gengamma_lcdf(q | shape, scale, k_shift);
  } else {
    log_q_F_T_q = negative_infinity();
    log_E_tF_T_q = negative_infinity();
  }

  // Unified form: F_{S+}(d) = (A - B) / w_P with
  //   A = d * F_T(d) + E * tilde F_T(q)
  //   B = q * F_T(q) + E * tilde F_T(d)
  // Ordering A >= B is guaranteed by F_{S+}(d) >= 0.
  real log_A = log_sum_exp(log(d) + log_F_T_d, log_E_tF_T_q);
  real log_B = log_sum_exp(log_q_F_T_q, log_E + log_tF_T_d);

  return log_diff_exp(log_A, log_B) - log_window;
}

// primarycensored_lcdf() with internal bounds [0 or -inf, inf] as the old
// vectorised function called it, routing the analytical uniform primary
// cases through the reference copies above.
real ref_primarycensored_lcdf(data real d, data int dist_id,
                              array[] real params, data real pwindow,
                              data int primary_id,
                              array[] real primary_params) {
  if (primary_id == 1 && d > 0) {
    real q = max({d - pwindow, 0});
    if (dist_id == 1) {
      return ref_primarycensored_lognormal_uniform_lcdf(d | q, params, pwindow);
    } else if (dist_id == 2) {
      return ref_primarycensored_gamma_uniform_lcdf(d | q, params, pwindow);
    } else if (dist_id == 3) {
      return ref_primarycensored_weibull_uniform_lcdf(d | q, params, pwindow);
    } else if (dist_id == 5) {
      return ref_primarycensored_gengamma_uniform_lcdf(d | q, params, pwindow);
    }
  }
  return primarycensored_lcdf(
    d | dist_id, params, pwindow,
    dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),
    positive_infinity(), primary_id, primary_params
  );
}

vector ref_primarycensored_sone_lpmf_vectorized(
  data int max_delay, data real L, data real D, data int dist_id,
  array[] real params, data real pwindow,
  data int primary_id, array[] real primary_params
) {

  int upper_interval = max_delay + 1;
  vector[upper_interval] log_pmfs;
  vector[upper_interval] log_cdfs;
  real log_normalizer;

  // Check if D is at least max_delay + 1
  if (D < upper_interval) {
    reject("D must be at least max_delay + 1");
  }

  // Compute log CDFs (without truncation normalization). The internal lower
  // bound below is 0 for positive-support delays and -inf otherwise; it is
  // inlined rather than bound to a local so Stan's type checker treats it as
  // data-only.
  // Start from max(1, floor(L)) to avoid computing unused CDFs when L > 0;
  // for L <= 0 (including -inf) start at 1 since F(d) = 0 for d <= 0.
  int start_idx = (!is_inf(L) && L > 0) ? max(1, to_int(floor(L))) : 1;
  for (d in start_idx:upper_interval) {
    log_cdfs[d] = ref_primarycensored_lcdf(
      d | dist_id, params, pwindow, primary_id, primary_params
    );
  }

  // Get CDF at lower truncation point L
  real log_cdf_L;
  if (is_inf(L)) {
    // No left truncation (L = -inf sentinel)
    log_cdf_L = negative_infinity();
  } else if (L >= 1 && L <= upper_interval && floor(L) == L) {
    // L is a positive integer within computed range, reuse cached value
    log_cdf_L = log_cdfs[to_int(L)];
  } else {
    // L is outside computed range or non-integer, compute directly
    log_cdf_L = ref_primarycensored_lcdf(
      L | dist_id, params, pwindow, primary_id, primary_params
    );
  }

  // Compute log normalizer: log(F(D) - F(L))
  real log_cdf_D;
  if (D > upper_interval) {
    if (is_inf(D)) {
      log_cdf_D = 0; // log(1) = 0 for infinite D
    } else {
      log_cdf_D = ref_primarycensored_lcdf(
        D | dist_id, params, pwindow, primary_id, primary_params
      );
    }
  } else {
    log_cdf_D = log_cdfs[upper_interval];
  }

  log_normalizer = primarycensored_log_normalizer(log_cdf_D, log_cdf_L, L);

  // Compute log PMFs: log((F(d) - F(d-1)) / (F(D) - F(L)))
  for (d in 1:upper_interval) {
    if (d <= L) {
      // Delay interval [d-1, d) is entirely at or below L
      log_pmfs[d] = negative_infinity();
    } else if (d - 1 < L) {
      // L falls within interval [d-1, d), so compute mass in [L, d)
      log_pmfs[d] = log_diff_exp(log_cdfs[d], log_cdf_L) - log_normalizer;
    } else if (d == 1 && dist_has_positive_support(dist_id)) {
      // First interval [0, 1) with L <= 0 and positive-support delay:
      // F(0) = 0, so PMF = F(1) / normalizer
      log_pmfs[d] = log_cdfs[d] - log_normalizer;
    } else if (d == 1) {
      // First interval [0, 1) with L <= 0 and real-support delay: F(0) is
      // non-zero in general, so compute it explicitly.
      real log_cdf_0 = ref_primarycensored_lcdf(
        0.0 | dist_id, params, pwindow, primary_id, primary_params
      );
      log_pmfs[d] = log_diff_exp(log_cdfs[d], log_cdf_0) - log_normalizer;
    } else {
      // Standard case: PMF = (F(d) - F(d-1)) / normalizer
      log_pmfs[d] = log_diff_exp(log_cdfs[d], log_cdfs[d-1]) - log_normalizer;
    }
  }

  return log_pmfs;
}
