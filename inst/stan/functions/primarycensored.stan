/**
  * Primary event censored distribution functions
  */

/**
  * Compute the log normalizer for truncation: log(F(D) - F(L))
  * @ingroup truncation_helpers
  *
  * @param log_cdf_D Log CDF at upper truncation point D
  * @param log_cdf_L Log CDF at lower truncation point L (negative_infinity if
  *   L = -inf, i.e. no lower truncation)
  * @param L Lower truncation point (-inf indicates no lower truncation)
  *
  * @return Log normalizer for truncation
  */
real primarycensored_log_normalizer(real log_cdf_D, real log_cdf_L, real L) {
  if (!is_inf(L)) {
    return log_diff_exp(log_cdf_D, log_cdf_L);
  } else {
    return log_cdf_D;
  }
}

/**
  * Apply truncation normalization to a log CDF value
  * @ingroup truncation_helpers
  *
  * Computes log((F(x) - F(L)) / (F(D) - F(L)))
  *
  * @param log_cdf Log CDF value to normalize
  * @param log_cdf_L Log CDF at lower truncation point L (negative_infinity if
  *   L = -inf, i.e. no lower truncation)
  * @param log_normalizer Log normalizer from primarycensored_log_normalizer
  * @param L Lower truncation point (-inf indicates no lower truncation)
  *
  * @return Normalized log CDF value
  */
real primarycensored_apply_truncation(real log_cdf, real log_cdf_L,
                                      real log_normalizer, real L) {
  if (!is_inf(L)) {
    return log_diff_exp(log_cdf, log_cdf_L) - log_normalizer;
  } else {
    return log_cdf - log_normalizer;
  }
}

/**
  * Compute log CDFs at both truncation bounds L and D
  * @ingroup truncation_helpers
  *
  * @param L Lower truncation point (-inf indicates no lower truncation)
  * @param D Upper truncation point (+inf indicates no upper truncation)
  * @param dist_id Distribution identifier
  * @param params Array of distribution parameters
  * @param pwindow Primary event window
  * @param primary_id Primary distribution identifier
  * @param primary_params Primary distribution parameters
  *
  * @return 2-element vector: [log_cdf_L, log_cdf_D]
  *
  * @note F(L) for finite L is computed via primarycensored_lcdf with internal
  *   bounds [0, +inf]. This is correct for the analytical solutions currently
  *   supported, all of which assume the delay distribution has non-negative
  *   support, so F(L) = 0 for any L <= 0. If future analytical solutions add
  *   distributions with negative support, this internal lower bound will need
  *   to be relaxed to -inf so the underlying CDF is evaluated directly.
  */
vector primarycensored_truncation_bounds(
  data real L, data real D,
  data int dist_id, array[] real params, data real pwindow,
  data int primary_id, array[] real primary_params
) {
  vector[2] result;
  // Internal lower bound for the un-truncated distribution: 0 lets the
  // `d <= L` early-exit in primarycensored_lcdf return -inf for delays below
  // the natural support of positive-support distributions; -inf disables that
  // short-circuit so distributions with support on the reals are integrated.
  // Expression is inlined (rather than bound to a local) so Stan's data-flow
  // checker recognises it as data-only.

  // Get CDF at lower truncation point L
  if (is_inf(L)) {
    result[1] = negative_infinity();
  } else {
    result[1] = primarycensored_lcdf(
      L | dist_id, params, pwindow,
      dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),
      positive_infinity(), primary_id, primary_params
    );
  }

  // Get CDF at upper truncation point D
  if (is_inf(D)) {
    result[2] = 0;
  } else {
    result[2] = primarycensored_lcdf(
      D | dist_id, params, pwindow,
      dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),
      positive_infinity(), primary_id, primary_params
    );
  }

  return result;
}

/**
  * Compute the primary event censored CDF for a single delay
  * @ingroup primary_censored_single
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
real primarycensored_cdf(data real d, data int dist_id, array[] real params,
                               data real pwindow, data real L, data real D,
                               data int primary_id,
                               array[] real primary_params) {
  real result;
  if (d <= L) {
    return 0;
  }

  if (d >= D) {
    return 1;
  }

  // Check if an analytical solution exists
  if (check_for_analytical(dist_id, primary_id)) {
    // Use analytical solution
    result = primarycensored_analytical_cdf(
      d | dist_id, params, pwindow, L, D, primary_id, primary_params
    );
  } else {
    // Use numerical integration for other cases. The integration variable
    // ranges over the primary-event time, so the natural lower bound is
    // d - pwindow. For positive-support delays the integrand `F_delay(t)` is
    // 0 for t <= 0, so an unclipped lower bound just adds a flat zero region
    // for negative t. Distributions with support on the reals also accept the
    // unclipped lower bound directly.
    real lower_bound = d - pwindow;
    int n_params = num_elements(params);
    int n_primary_params = num_elements(primary_params);
    array[n_params + n_primary_params] real theta = append_array(params, primary_params);
    array[4] int ids = {dist_id, primary_id, n_params, n_primary_params};

    vector[1] y0 = rep_vector(0.0, 1);
    result = ode_rk45(primarycensored_ode, y0, lower_bound, {d}, theta, {d, pwindow}, ids)[1, 1];

    // Apply truncation normalization on log scale for numerical stability
    if (!is_inf(D) || !is_inf(L)) {
      real log_result = log(result);
      vector[2] bounds = primarycensored_truncation_bounds(
        L, D, dist_id, params, pwindow, primary_id, primary_params
      );
      real log_cdf_L = bounds[1];
      real log_cdf_D = bounds[2];

      real log_normalizer = primarycensored_log_normalizer(log_cdf_D, log_cdf_L, L);
      log_result = primarycensored_apply_truncation(
        log_result, log_cdf_L, log_normalizer, L
      );
      result = exp(log_result);
    }
  }

  return result;
}

/**
  * Compute the primary event censored log CDF for a single delay
  * @ingroup primary_censored_single
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
  *
  * @code
  * // Example: Weibull delay distribution with uniform primary distribution
  * real d = 3.0;
  * int dist_id = 3; // Weibull
  * array[2] real params = {2.0, 1.5}; // shape and scale
  * real pwindow = 1.0;
  * real L = 0.0;
  * real D = positive_infinity();
  * int primary_id = 1; // Uniform
  * array[0] real primary_params = rep_array(0.0, 0);
  * real log_cdf = primarycensored_lcdf(
  *   d, dist_id, params, pwindow, L, D, primary_id, primary_params
  * );
  * @endcode
  */
real primarycensored_lcdf(data real d, data int dist_id, array[] real params,
                                data real pwindow, data real L, data real D,
                                data int primary_id,
                                array[] real primary_params) {
  real result;

  if (d <= L) {
    return negative_infinity();
  }

  if (d >= D) {
    return 0;
  }

  // Check if an analytical solution exists. The internal lower bound is 0 for
  // positive-support delays (lets the d <= L early-exit return -inf for d <= 0)
  // and -inf for distributions with support on the reals.
  if (check_for_analytical(dist_id, primary_id)) {
    result = primarycensored_analytical_lcdf(
      d | dist_id, params, pwindow,
      dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),
      positive_infinity(), primary_id, primary_params
    );
  } else {
    // Use numerical integration
    result = log(primarycensored_cdf(
      d | dist_id, params, pwindow,
      dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),
      positive_infinity(), primary_id, primary_params
    ));
  }

  // Handle truncation normalization
  if (!is_inf(D) || !is_inf(L)) {
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
  * Compute the primary event censored log PMF for a single delay
  * @ingroup primary_censored_single
  *
  * @param d Delay (integer)
  * @param dist_id Distribution identifier
  * @param params Array of distribution parameters
  * @param pwindow Primary event window
  * @param d_upper Upper bound for the delay interval
  * @param L Minimum delay (lower truncation point)
  * @param D Maximum delay (upper truncation point)
  * @param primary_id Primary distribution identifier
  * @param primary_params Primary distribution parameters
  *
  * @return Primary event censored log PMF, normalized over [L, D] if truncation
  * is applied
  *
  * @code
  * // Example: Weibull delay distribution with uniform primary distribution
  * int d = 3;
  * int dist_id = 3; // Weibull
  * array[2] real params = {2.0, 1.5}; // shape and scale
  * real pwindow = 1.0;
  * real d_upper = 4.0;
  * real L = 0.0;
  * real D = positive_infinity();
  * int primary_id = 1; // Uniform
  * array[0] real primary_params = {};
  * real log_pmf = primarycensored_lpmf(
  *   d, dist_id, params, pwindow, d_upper, L, D, primary_id, primary_params
  * );
  * @endcode
  */
real primarycensored_lpmf(data int d, data int dist_id, array[] real params,
                                data real pwindow, data real d_upper,
                                data real L, data real D, data int primary_id,
                                array[] real primary_params) {
  if (d_upper > D) {
    reject("Upper truncation point is greater than D. It is ", d_upper,
           " and D is ", D, ". Resolve this by increasing D to be greater or equal to d + swindow or decreasing swindow.");
  }
  if (d_upper <= d) {
    reject("Upper truncation point is less than or equal to d. It is ", d_upper,
           " and d is ", d, ". Resolve this by increasing d to be less than d_upper.");
  }
  if (d < L) {
    return negative_infinity();
  }
  real log_cdf_upper = primarycensored_lcdf(
    d_upper | dist_id, params, pwindow,
    dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),
    positive_infinity(), primary_id, primary_params
  );
  real log_cdf_lower = primarycensored_lcdf(
    d | dist_id, params, pwindow,
    dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),
    positive_infinity(), primary_id, primary_params
  );

  // Apply truncation normalization: log((F(d_upper) - F(d)) / (F(D) - F(L)))
  if (!is_inf(D) || !is_inf(L)) {
    real log_cdf_D;
    real log_cdf_L;

    // Get CDF at lower truncation point L
    if (is_inf(L)) {
      // No left truncation (L = -inf sentinel)
      log_cdf_L = negative_infinity();
    } else if (d == L) {
      // Reuse already computed CDF at d
      log_cdf_L = log_cdf_lower;
    } else {
      // Compute CDF at L directly
      log_cdf_L = primarycensored_lcdf(
        L | dist_id, params, pwindow,
        dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),
        positive_infinity(), primary_id, primary_params
      );
    }

    // Get CDF at upper truncation point D
    if (d_upper == D) {
      log_cdf_D = log_cdf_upper;
    } else if (is_inf(D)) {
      log_cdf_D = 0;
    } else {
      log_cdf_D = primarycensored_lcdf(
        D | dist_id, params, pwindow,
        dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),
        positive_infinity(), primary_id, primary_params
      );
    }

    real log_normalizer = primarycensored_log_normalizer(log_cdf_D, log_cdf_L, L);
    return log_diff_exp(log_cdf_upper, log_cdf_lower) - log_normalizer;
  } else {
    return log_diff_exp(log_cdf_upper, log_cdf_lower);
  }
}

/**
  * Compute the primary event censored PMF for a single delay
  * @ingroup primary_censored_single
  *
  * @param d Delay (integer)
  * @param dist_id Distribution identifier
  * @param params Array of distribution parameters
  * @param pwindow Primary event window
  * @param d_upper Upper bound for the delay interval
  * @param L Minimum delay (lower truncation point)
  * @param D Maximum delay (upper truncation point)
  * @param primary_id Primary distribution identifier
  * @param primary_params Primary distribution parameters
  *
  * @return Primary event censored PMF, normalized over [L, D] if truncation
  * is applied
  *
  * @code
  * // Example: Weibull delay distribution with uniform primary distribution
  * int d = 3;
  * int dist_id = 3; // Weibull
  * array[2] real params = {2.0, 1.5}; // shape and scale
  * real pwindow = 1.0;
  * real swindow = 1.0;
  * real d_upper = d + swindow; // = 4.0
  * real L = 0.0;
  * real D = positive_infinity();
  * int primary_id = 1; // Uniform
  * array[0] real primary_params = {};
  * real pmf = primarycensored_pmf(d, dist_id, params, pwindow, d_upper, L, D, primary_id, primary_params);
  * @endcode
  */
real primarycensored_pmf(data int d, data int dist_id, array[] real params,
                               data real pwindow, data real d_upper,
                               data real L, data real D, data int primary_id,
                               array[] real primary_params) {
  return exp(
    primarycensored_lpmf(
      d | dist_id, params, pwindow, d_upper, L, D, primary_id, primary_params
    )
  );
}

/**
  * Compute the primary event censored log PMF for integer delays up to max_delay
  * @ingroup primary_censored_vectorized
  *
  * @param max_delay Maximum delay to compute PMF for
  * @param L Minimum delay (lower truncation point)
  * @param D Maximum delay (upper truncation point), must be at least max_delay + 1
  * @param dist_id Distribution identifier
  * @param params Array of distribution parameters
  * @param pwindow Primary event window
  * @param primary_id Primary distribution identifier
  * @param primary_params Primary distribution parameters
  *
  * @return Vector of primary event censored log PMFs for delays \[0, 1\] to
  * \[max_delay, max_delay + 1\].
  *
  * This function differs from primarycensored_lpmf in that it:
  * 1. Computes PMFs for all integer delays from \[0, 1\] to \[max_delay,
  *    max_delay + 1\] in one call.
  * 2. Assumes integer delays (swindow = 1)
  * 3. Is more computationally efficient for multiple delay calculation as it
  *    reduces the number of integration calls.
  *
  * @code
  * // Example: Weibull delay distribution with uniform primary distribution
  * int max_delay = 10;
  * real L = 0.0;
  * real D = 15.0;
  * int dist_id = 3; // Weibull
  * array[2] real params = {2.0, 1.5}; // shape and scale
  * real pwindow = 7.0;
  * int primary_id = 1; // Uniform
  * array[0] real primary_params = {};

  * vector[max_delay] log_pmf =
  *   primarycensored_sone_lpmf_vectorized(
  *      max_delay, L, D, dist_id, params, pwindow, primary_id,
  *      primary_params
  *   );
  * @endcode
  */
vector primarycensored_sone_lpmf_vectorized(
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
    log_cdfs[d] = primarycensored_lcdf(
      d | dist_id, params, pwindow,
      dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),
      positive_infinity(), primary_id, primary_params
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
    log_cdf_L = primarycensored_lcdf(
      L | dist_id, params, pwindow,
      dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),
      positive_infinity(), primary_id, primary_params
    );
  }

  // Compute log normalizer: log(F(D) - F(L))
  real log_cdf_D;
  if (D > upper_interval) {
    if (is_inf(D)) {
      log_cdf_D = 0; // log(1) = 0 for infinite D
    } else {
      log_cdf_D = primarycensored_lcdf(
        D | dist_id, params, pwindow,
        dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),
        positive_infinity(), primary_id, primary_params
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
    } else if (d == 1) {
      // First interval [0, 1) with L <= 0: F(0) = 0, so PMF = F(1) / normalizer
      log_pmfs[d] = log_cdfs[d] - log_normalizer;
    } else {
      // Standard case: PMF = (F(d) - F(d-1)) / normalizer
      log_pmfs[d] = log_diff_exp(log_cdfs[d], log_cdfs[d-1]) - log_normalizer;
    }
  }

  return log_pmfs;
}

/**
  * Compute the primary event censored PMF for integer delays up to max_delay
  * @ingroup primary_censored_vectorized
  *
  * @param max_delay Maximum delay to compute PMF for
  * @param L Minimum delay (lower truncation point)
  * @param D Maximum delay (upper truncation point), must be at least max_delay + 1
  * @param dist_id Distribution identifier
  * @param params Array of distribution parameters
  * @param pwindow Primary event window
  * @param primary_id Primary distribution identifier
  * @param primary_params Primary distribution parameters
  *
  * @return Vector of primary event censored PMFs for integer delays 1 to
  * max_delay
  *
  * This function differs from primarycensored_pmf in that it:
  * 1. Computes PMFs for all integer delays from \[0, 1\] to \[max_delay,
  *    max_delay + 1\] in one call.
  * 2. Assumes integer delays (swindow = 1)
  * 3. Is more computationally efficient for multiple delay calculations
  *
  * @code
  * // Example: Weibull delay distribution with uniform primary distribution
  * int max_delay = 10;
  * real L = 0.0;
  * real D = 15.0;
  * int dist_id = 3; // Weibull
  * array[2] real params = {2.0, 1.5}; // shape and scale
  * real pwindow = 7.0;
  * int primary_id = 1; // Uniform
  * array[0] real primary_params = {};
  * vector[max_delay] pmf =
  *   primarycensored_sone_lpmf_vectorized(
  *      max_delay, L, D, dist_id, params, pwindow, primary_id, primary_params
  *   );
  * @endcode
  */
vector primarycensored_sone_pmf_vectorized(
  data int max_delay, data real L, data real D, data int dist_id,
  array[] real params, data real pwindow,
  data int primary_id,
  array[] real primary_params
) {
  return exp(
    primarycensored_sone_lpmf_vectorized(
      max_delay, L, D, dist_id, params, pwindow, primary_id, primary_params
    )
  );
}
