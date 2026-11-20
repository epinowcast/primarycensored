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
  * Log of the incomplete gamma series for the tilt transform of a Stacy
  * family delay
  * @ingroup tilt_transforms
  *
  * The Weibull (k = 1) and generalised gamma have CDF P(k, (t / scale)^shape)
  * and transform
  *   T_f(xi; t) = sum_n (xi scale)^n / n! Gamma(k + n / shape) / Gamma(k)
  *                P(k + n / shape, (t / scale)^shape).
  * The terms for negative xi alternate in sign, and the sum of their absolute
  * values is the transform at -xi. Positive and negative terms are summed
  * apart until the remainder, bounded by a geometric series once
  * n + 1 > |xi| t, is below 1e-17 of that sum.
  *
  * @param t Point, positive
  * @param xi Tilt, not zero
  * @param shape Power of the Stacy family
  * @param scale Scale
  * @param k Shape of the underlying gamma
  *
  * @return Vector [log T_f(xi; t), log of the sum of the absolute terms over
  *   the result], or `nan` for both where the sum of the absolute terms is
  *   more than 1e5 times the result, |xi| t is 400 or more, or 400 terms do
  *   not converge
  */
vector primarycensored_stacy_series(real t, real xi, real shape, real scale,
                                    real k) {
  real log_x = shape * (log(t) - log(scale));
  int alternating = xi < 0;
  real log_z = log(abs(xi)) + log(scale);
  real log_w = log(abs(xi)) + log(t);
  real w = exp(log_w);
  vector[2] unreliable = rep_vector(not_a_number(), 2);
  if (w >= 400) return unreliable;
  real log_k = lgamma(k);
  real log_pos = 0;
  real log_neg = negative_infinity();
  real log_first = 0;
  for (n in 0:400) {
    real s = k + n / shape;
    real log_term = n * log_z - lgamma(n + 1) + lgamma(s) - log_k
                    + gamma_lcdf_logx(log_x, s);
    if (n == 0) {
      log_first = log_term;
      log_pos = log_term;
    } else if (alternating && n % 2 == 1) {
      log_neg = log_sum_exp(log_neg, log_term);
    } else {
      log_pos = log_sum_exp(log_pos, log_term);
    }
    real log_abs = log_sum_exp(log_pos, log_neg);
    // The transform is at most the CDF for negative xi
    if (alternating && log_abs - log_first > log(1e5)) return unreliable;
    if (n + 1 > w
        && log_term + log_w - log(n + 1 - w) < log_abs + log(1e-17)) {
      if (log_pos <= log_neg) return unreliable;
      real log_sum = primarycensored_log_diff_exp(log_pos, log_neg);
      if (log_abs - log_sum > log(1e5)) return unreliable;
      return [log_sum, log_abs - log_sum]';
    }
  }
  return unreliable;
}

/**
  * Check if the tilt transform is closed form for a delay and tilt
  * @ingroup tilt_transforms
  *
  * The exponential (4) and gamma (2) forms need the delay with the rate
  * lowered by xi, which exists if rate - xi > 0. The normal (18) form and the
  * Weibull (3) and generalised gamma (5) series have no restriction. The
  * lognormal (1) form needs xi sigma^2 exp(mu) below exp(690) in magnitude.
  * Callers use the numerical path when this is 0.
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
  if (dist_id == 18 || dist_id == 3 || dist_id == 5) return 1;
  if (dist_id == 1) {
    if (params[2] <= 0) return 0;
    if (xi >= 0) return 1;
    return log(-xi) + 2 * log(params[2]) + params[1] < 690;
  }
  return 0;
}

/**
  * Check if the tilt transform is the path to use at a point
  * @ingroup tilt_transforms
  *
  * This is check_for_tilt_transform() and, for the lognormal (1) with
  * xi > 0, that the series is the path to use at t. The series costs about
  * xi t + 9 sqrt(xi t) + 30 terms and is limited to 20000 terms. It is used
  * up to xi t of 60, and beyond that only where xi w is above 2, as the ODE
  * is inaccurate in the lower tail there. The cut-off of 60 is set for the
  * ODE and differs from the one in R, which is set for the R numerical path.
  *
  * @param dist_id Distribution identifier for the delay distribution
  * @param xi Tilt. The exponentially tilted window with tilt rho needs
  *   xi = -rho
  * @param params Array of distribution parameters, as for dist_lcdf()
  * @param t Point, the largest of the points if there are several
  * @param pwindow Primary event window
  *
  * @return 1 if the transform is closed form and is the path to use at t,
  * 0 otherwise
  */
int check_for_tilt_transform_at(int dist_id, real xi, array[] real params,
                                data real t, data real pwindow) {
  if (!check_for_tilt_transform(dist_id, xi, params)) return 0;
  if (dist_id == 1 && xi > 0 && t > 0) {
    if (xi * t + 9 * sqrt(xi * t) + 30 > 20000) return 0;
    return xi * t <= 60 || xi * pwindow > 2;
  }
  return 1;
}

/**
  * Nodes of the 32 point Gauss-Legendre rule on [-1, 1]
  * @ingroup tilt_transforms
  *
  * @return The 16 positive nodes, the others are their negatives
  */
vector primarycensored_gauss_legendre_nodes() {
  return [
    0.048307665687738324, 0.14447196158279649, 0.23928736225213706,
    0.33186860228212767, 0.42135127613063533, 0.50689990893222936,
    0.5877157572407623, 0.66304426693021523, 0.73218211874028971,
    0.79448379596794239, 0.84936761373256997, 0.89632115576605209,
    0.93490607593773967, 0.96476225558750639, 0.98561151154526838,
    0.99726386184948157
  ]';
}

/**
  * Log weights of the 32 point Gauss-Legendre rule on [-1, 1]
  * @ingroup tilt_transforms
  *
  * @return The log weights of the 16 positive nodes, shared by the negatives
  */
vector primarycensored_gauss_legendre_log_weights() {
  return [
    -2.3377969318791947, -2.3471775191742132, -2.3661171972104831,
    -2.3949868407368724, -2.4343797887852916, -2.4851635920714812,
    -2.5485636934665021, -2.626297960160425, -2.7207977655059796,
    -2.8355865687925124, -2.9759677006607701, -3.1503787891021764,
    -3.3733722294882642, -3.6733185431605069, -4.1181622816923937,
    -4.9591900848988502
  ]';
}

/**
  * Log of a Gaussian weighted integral of the tilt over one panel
  * @ingroup tilt_transforms
  *
  * The log of the integral of exp(-rho exp(mu + sigma z)) phi(z) over [a, b]
  * by the 32 point Gauss-Legendre rule, where phi is the standard normal
  * density.
  *
  * @param a Lower limit
  * @param b Upper limit
  * @param mu Mean of the log of the delay
  * @param sigma Standard deviation of the log of the delay
  * @param rho Tilt, at least 0
  *
  * @return Log of the integral, `-inf` if b <= a
  */
real primarycensored_lognormal_tilt_one_panel(real a, real b, real mu,
                                              real sigma, real rho) {
  if (b <= a) return negative_infinity();
  real half = 0.5 * (b - a);
  real mid = 0.5 * (a + b);
  vector[16] nodes = primarycensored_gauss_legendre_nodes();
  vector[16] log_weights = primarycensored_gauss_legendre_log_weights();
  vector[32] terms;
  for (i in 1:16) {
    real z_up = mid + half * nodes[i];
    real z_down = mid - half * nodes[i];
    terms[i] = log_weights[i] - rho * exp(mu + sigma * z_up)
               - 0.5 * square(z_up);
    terms[16 + i] = log_weights[i] - rho * exp(mu + sigma * z_down)
                    - 0.5 * square(z_down);
  }
  return log_sum_exp(terms) + log(half) - 0.5 * log(2 * pi());
}

/**
  * Log of a Gaussian weighted integral of the tilt over a range
  * @ingroup tilt_transforms
  *
  * As primarycensored_lognormal_tilt_one_panel(), with the range split into
  * ceil(sigma / 1.8) equal panels. One panel loses accuracy for sigma above
  * 1.8.
  *
  * @param a Lower limit
  * @param b Upper limit
  * @param mu Mean of the log of the delay
  * @param sigma Standard deviation of the log of the delay
  * @param rho Tilt, at least 0
  *
  * @return Log of the integral, `-inf` if b <= a
  */
real primarycensored_lognormal_tilt_panel(real a, real b, real mu,
                                          real sigma, real rho) {
  if (b <= a) return negative_infinity();
  int n_panels = 1;
  while (n_panels * 1.8 < sigma) {
    n_panels += 1;
  }
  if (n_panels == 1) {
    return primarycensored_lognormal_tilt_one_panel(a, b, mu, sigma, rho);
  }
  real width = (b - a) / n_panels;
  vector[n_panels] panels;
  for (j in 1:n_panels) {
    panels[j] = primarycensored_lognormal_tilt_one_panel(
      a + (j - 1) * width, a + j * width, mu, sigma, rho
    );
  }
  return log_sum_exp(panels);
}

/**
  * Mode, limits and total of the tilted lognormal integrand
  * @ingroup tilt_transforms
  *
  * The log integrand in z = (log u - mu) / sigma is
  * -rho exp(mu + sigma z) - z^2 / 2. It is concave with a single mode
  * z0 = -W0(rho sigma^2 exp(mu)) / sigma. Beyond the limits lo and hi it is
  * below exp(-40) of its peak. None of these depends on the point, so they
  * are computed once per call, see log_tilt_transform_context().
  *
  * @param mu Mean of the log of the delay
  * @param sigma Standard deviation of the log of the delay
  * @param rho Tilt, above 0
  *
  * @return Vector [z0, lo, hi, log total]
  */
vector primarycensored_lognormal_tilt_bump(real mu, real sigma, real rho) {
  real w0 = lambert_w0(exp(log(rho) + 2 * log(sigma) + mu));
  real z0 = -w0 / sigma;
  real u0 = exp(mu - w0);
  real hi = fmin(
    z0 + sqrt(80 / (1 + w0)), fmax((log(u0 + 40 / rho) - mu) / sigma, -z0)
  );
  real lo = fmax(
    z0 - sqrt(80), -sqrt(square(z0) + 2 * w0 / square(sigma) + 80)
  );
  real log_total = log_sum_exp(
    primarycensored_lognormal_tilt_panel(lo, z0, mu, sigma, rho),
    primarycensored_lognormal_tilt_panel(z0, hi, mu, sigma, rho)
  );
  return [z0, lo, hi, log_total]';
}

/**
  * Log of the lognormal tilt transform by quadrature for a negative tilt
  * @ingroup tilt_transforms
  *
  * The lower and upper transform for xi = -rho < 0. The transform on the
  * side of the mode that t is on is one panel that starts or ends at t. The
  * other is the difference from the total, which does not cancel as it is at
  * least the mass on the far side of the mode.
  *
  * @param t Point, positive
  * @param mu Mean of the log of the delay
  * @param sigma Standard deviation of the log of the delay
  * @param rho Tilt, above 0
  * @param bump Output of primarycensored_lognormal_tilt_bump()
  *
  * @return Vector [log T_f(-rho; t), log(T_f(-rho; Inf) - T_f(-rho; t))]
  */
vector primarycensored_lognormal_tilt_quadrature(real t, real mu,
                                                 real sigma, real rho,
                                                 vector bump) {
  real z = (log(t) - mu) / sigma;
  // The slope of the log integrand at t, positive left of the mode
  real slope = -rho * sigma * t - z;
  if (z < bump[1]) {
    real width = 80 / (slope + sqrt(square(slope) + 80));
    real log_lower = primarycensored_lognormal_tilt_panel(
      z - width, z, mu, sigma, rho
    );
    return [
      log_lower, primarycensored_log_diff_exp(bump[4], log_lower)
    ]';
  }
  real curvature = 1 + rho * square(sigma) * t;
  real width = 80 / (abs(slope) + sqrt(square(slope) + 80 * curvature));
  // The integrand is below exp(-40) of its value at t beyond the point where
  // the tilt alone has decayed by 40
  width = fmin(
    width, fmax(abs(z), (log(t + 40 / rho) - mu) / sigma) - z
  );
  real log_upper = primarycensored_lognormal_tilt_panel(
    z, z + width, mu, sigma, rho
  );
  return [primarycensored_log_diff_exp(bump[4], log_upper), log_upper]';
}

/**
  * Log of the lognormal tilt transform by series for a positive tilt
  * @ingroup tilt_transforms
  *
  * The lower transform for xi > 0 as the sum of the positive terms
  * xi^k m_k(t) / k! over the partial moments
  * m_k(t) = exp(k mu + k^2 sigma^2 / 2) Phi((log t - mu) / sigma - k sigma).
  * The sum stops once past k = xi t a term is below exp(-41) of the sum. The
  * total diverges, so there is no upper transform. It rejects past 20000
  * terms, see check_for_tilt_transform_at().
  *
  * @param t Point, positive
  * @param mu Mean of the log of the delay
  * @param sigma Standard deviation of the log of the delay
  * @param xi Tilt, above 0
  *
  * @return log T_f(xi; t)
  */
real primarycensored_lognormal_tilt_series(real t, real mu, real sigma,
                                           real xi) {
  real z = (log(t) - mu) / sigma;
  real log_xi_mu = log(xi) + mu;
  real log_factorial = 0;
  real total = primarycensored_log_std_normal_cdf(z);
  for (k in 1:20000) {
    log_factorial += log(k);
    real log_term = k * log_xi_mu + 0.5 * square(k * sigma) - log_factorial
                    + primarycensored_log_std_normal_cdf(z - k * sigma);
    total = log_sum_exp(total, log_term);
    if (k > xi * t && log_term < total - 41) return total;
  }
  reject(
    "The lognormal tilt transform needs more than 20000 terms for t = ", t,
    " and xi = ", xi, ". Use the numerical path."
  );
}

/**
  * Log of the lognormal tilt transform over the lower and upper support
  * @ingroup tilt_transforms
  *
  * The upper transform is `inf` for xi > 0, as the total diverges.
  *
  * @param t Point
  * @param mu Mean of the log of the delay
  * @param sigma Standard deviation of the log of the delay
  * @param xi Tilt
  * @param bump Output of primarycensored_lognormal_tilt_bump() for -xi,
  *   used for xi < 0
  *
  * @return Vector [log T_f(xi; t), log(T_f(xi; Inf) - T_f(xi; t))]
  */
vector primarycensored_lognormal_tilt_pair(real t, real mu, real sigma,
                                           real xi, vector bump) {
  if (xi == 0) {
    if (t <= 0) return [negative_infinity(), 0]';
    real z = (log(t) - mu) / sigma;
    return [
      primarycensored_log_std_normal_cdf(z),
      primarycensored_log_std_normal_cdf(-z)
    ]';
  }
  if (xi > 0) {
    if (t <= 0) return [negative_infinity(), positive_infinity()]';
    return [
      primarycensored_lognormal_tilt_series(t, mu, sigma, xi),
      positive_infinity()
    ]';
  }
  if (t <= 0) return [negative_infinity(), bump[4]]';
  return primarycensored_lognormal_tilt_quadrature(t, mu, sigma, -xi, bump);
}

/**
  * Quantities shared by the tilt transforms of a delay at every point
  * @ingroup tilt_transforms
  *
  * Computed once per call and passed to log_tilt_transform_pair_shared(). It
  * has length 4 and is 0 for delays whose transform has no shared setup.
  *
  * @param dist_id Distribution identifier, see check_for_tilt_transform()
  * @param xi Tilt
  * @param params Array of distribution parameters, as for dist_lcdf()
  *
  * @return Vector of 4 shared quantities
  */
vector log_tilt_transform_context(int dist_id, real xi, array[] real params) {
  if (dist_id == 1 && xi < 0) {
    return primarycensored_lognormal_tilt_bump(params[1], params[2], -xi);
  }
  return rep_vector(0, 4);
}

/**
  * Log of the tilt transform over the lower and the upper part of the support
  * @ingroup tilt_transforms
  *
  * The lower transform is T_f(xi; t), the delay CDF for xi = 0, and `-inf`
  * for t <= 0 for delays on the non-negative reals. The upper transform is
  * T_f(xi; Inf) - T_f(xi; t). The third element is the log of the factor by
  * which the terms of the lower transform exceed it, 0 unless the transform is
  * a series. Only defined where check_for_tilt_transform() is 1.
  *
  * The gamma is the gamma CDF with the rate lowered by xi times the total
  * (rate / (rate - xi))^shape, with both tails from gamma_lcdf_lccdf_logx().
  *
  * The lognormal (1) has no closed form, see
  * primarycensored_lognormal_tilt_pair(). Its upper transform is `inf` for
  * xi > 0, where the total diverges.
  *
  * @param t Point
  * @param dist_id Distribution identifier: 1 (Lognormal), 2 (Gamma),
  *   3 (Weibull), 4 (Exponential), 5 (Generalised gamma) or 18 (Normal), see
  *   check_for_tilt_transform()
  * @param xi Tilt
  * @param params Array of distribution parameters, as for dist_lcdf()
  * @param context Output of log_tilt_transform_context() for dist_id, xi
  *   and params
  *
  * @return Vector [log T_f(xi; t), log(T_f(xi; Inf) - T_f(xi; t)), log loss].
  *   The Weibull and generalised gamma have `nan` for the lower transform
  *   where the series is not reliable, and for the upper transform if xi is
  *   not 0.
  */
vector log_tilt_transform_pair_shared(real t, int dist_id, real xi,
                                      array[] real params, vector context) {
  if (dist_id == 2) {
    real shape = params[1];
    real rate = params[2];
    real log_total = -shape * log1m(xi / rate);
    if (t <= 0) return [negative_infinity(), log_total, 0]';
    return append_row(
      log_total + gamma_lcdf_lccdf_logx(log(t) + log(rate - xi), shape), 0
    );
  } else if (dist_id == 4) {
    real rate = params[1];
    real tilted_rate = rate - xi;
    real log_total = -log1m(xi / rate);
    if (t <= 0) return [negative_infinity(), log_total, 0]';
    return [
      log_total + log1m_exp(-tilted_rate * t), log_total - tilted_rate * t, 0
    ]';
  } else if (dist_id == 18) {
    // The upper tail is the lower tail of the reflected normal
    real mu = params[1];
    real sigma = params[2];
    real z = (t - mu - xi * square(sigma)) / sigma;
    real log_total = xi * mu + 0.5 * square(xi * sigma);
    return [
      log_total + primarycensored_log_std_normal_cdf(z),
      log_total + primarycensored_log_std_normal_cdf(-z), 0
    ]';
  } else if (dist_id == 3 || dist_id == 5) {
    real shape = params[1];
    real scale = params[2];
    real k = dist_id == 3 ? 1.0 : params[3];
    if (xi == 0) {
      if (t <= 0) return [negative_infinity(), 0, 0]';
      return append_row(
        gamma_lcdf_lccdf_logx(shape * (log(t) - log(scale)), k), 0
      );
    }
    if (t <= 0) return [negative_infinity(), not_a_number(), 0]';
    vector[2] series = primarycensored_stacy_series(t, xi, shape, scale, k);
    return [series[1], not_a_number(), series[2]]';
  } else if (dist_id == 1) {
    vector[2] pair = primarycensored_lognormal_tilt_pair(
      t, params[1], params[2], xi, context
    );
    return append_row(pair, 0);
  }
  reject("Invalid distribution identifier: ", dist_id);
}

/**
  * Log of the tilt transform over the lower and the upper part of the support
  * @ingroup tilt_transforms
  *
  * As log_tilt_transform_pair_shared() with the shared quantities computed
  * for this point alone.
  *
  * @param t Point
  * @param dist_id Distribution identifier, see log_tilt_transform_pair_shared()
  * @param xi Tilt
  * @param params Array of distribution parameters, as for dist_lcdf()
  *
  * @return Vector [log T_f(xi; t), log(T_f(xi; Inf) - T_f(xi; t)), log loss]
  */
vector log_tilt_transform_pair(real t, int dist_id, real xi,
                               array[] real params) {
  return log_tilt_transform_pair_shared(
    t, dist_id, xi, params, log_tilt_transform_context(dist_id, xi, params)
  );
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
  * Log moments of a Stacy family delay about a point
  * @ingroup tilt_transforms
  *
  * As primarycensored_gamma_tilt_moments(), from the partial moments
  * int_0^t u^j f(u) du = scale^j Gamma(k + j / shape) / Gamma(k)
  * P(k + j / shape, (t / scale)^shape).
  *
  * @param t Point, positive
  * @param shape Power of the Stacy family
  * @param scale Scale
  * @param k Shape of the underlying gamma
  *
  * @return Vector [log G_1(t), log G_2(t), log G_3(t)]
  */
vector primarycensored_stacy_tilt_moments(real t, real shape, real scale,
                                          real k) {
  real log_x = shape * (log(t) - log(scale));
  real k1 = k + 1 / shape;
  real k2 = k + 2 / shape;
  real k3 = k + 3 / shape;
  real log_k = lgamma(k);
  real log_t = log(t);
  real log_m0 = gamma_lcdf_logx(log_x, k);
  real log_m1 = log(scale) + lgamma(k1) - log_k
                + gamma_lcdf_logx(log_x, k1);
  real log_m2 = 2 * log(scale) + lgamma(k2) - log_k
                + gamma_lcdf_logx(log_x, k2);
  real log_m3 = 3 * log(scale) + lgamma(k3) - log_k
                + gamma_lcdf_logx(log_x, k3);
  return primarycensored_moments_from_partial(
    log_t, [log_m0, log_m1, log_m2, log_m3]'
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
  * G_3 = sigma^3 (z (z^2 + 3) Phi(z) + (z^2 + 2) phi(z)). The lognormal
  * uses its partial moments
  * m_k(t) = exp(k mu + k^2 sigma^2 / 2) Phi(z - k sigma).
  *
  * @param t Point
  * @param dist_id Distribution identifier: 1 (Lognormal), 2 (Gamma),
  *   3 (Weibull), 4 (Exponential), 5 (Generalised gamma) or 18 (Normal), see
  *   check_for_tilt_transform()
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
  } else if (dist_id == 3 || dist_id == 5) {
    if (t <= 0) return rep_vector(negative_infinity(), 3);
    return primarycensored_stacy_tilt_moments(
      t, params[1], params[2], dist_id == 3 ? 1.0 : params[3]
    );
  } else if (dist_id == 1) {
    if (t <= 0) return rep_vector(negative_infinity(), 3);
    real mu = params[1];
    real sigma = params[2];
    real log_t = log(t);
    real z = (log_t - mu) / sigma;
    vector[4] log_m;
    for (k in 0:3) {
      log_m[k + 1] = k * mu + 0.5 * square(k * sigma)
                     + primarycensored_log_std_normal_cdf(z - k * sigma);
    }
    return primarycensored_moments_from_partial(log_t, log_m);
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
  return primary_id == 2
         && (dist_id == 1 || dist_id == 2 || dist_id == 3 || dist_id == 4
             || dist_id == 5 || dist_id == 18);
}

/**
  * Compute the exponentially tilted terms at an endpoint
  * @ingroup exponential_tilt_solutions
  *
  * @param t Endpoint, d or d - pwindow
  * @param dist_id Distribution identifier, see check_for_exptilt()
  * @param rho Tilt
  * @param params Array of distribution parameters
  * @param context Output of log_tilt_transform_context() for -rho
  *
  * @return Vector [log F(t), log(1 - F(t)), log J(t), log(J(Inf) - J(t)),
  * log loss of J], see log_tilt_transform_pair_shared(). The upper tail term
  * of J is `inf` where the total diverges, as for a lognormal delay with
  * rho < 0.
  */
vector primarycensored_exptilt_terms(real t, int dist_id, real rho,
                                     array[] real params, vector context) {
  return append_row(
    log_tilt_transform_pair(t, dist_id, 0, params)[1:2],
    log_tilt_transform_pair_shared(t, dist_id, -rho, params, context)
  );
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
  * @param upper_d Log upper tail quantity at d, `nan` if there is none
  * @param upper_q Log upper tail quantity at q, `nan` if there is none
  *
  * @return Log of the difference, `-inf` if it is zero to rounding. The
  * lower tail terms are used where the upper tail terms are `inf`.
  */
real primarycensored_tail_diff(real lower_d, real lower_q, real upper_d,
                               real upper_q) {
  if (is_nan(upper_d) || is_nan(upper_q)) {
    return primarycensored_log_diff_exp(lower_d, lower_q);
  }
  // NaN from terms that both underflow is a zero difference
  if (is_nan(lower_q - lower_d) || is_inf(upper_q)
      || lower_q - lower_d <= upper_d - upper_q) {
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
  * The result is `nan` if J at d or q is. Where J has no upper tail, it is
  * also `nan` if the terms of the numerator over the denominator exceed 1e5.
  *
  * @param terms_d Terms at d, from primarycensored_exptilt_terms()
  * @param terms_q Terms at q = d - pwindow
  * @param d Delay
  * @param rho Tilt, not zero
  * @param pwindow Primary event window
  *
  * @return Log of the primary event censored CDF at d, or `nan`
  */
real primarycensored_exptilt_lcdf_from_terms(vector terms_d, vector terms_q,
                                             data real d, real rho,
                                             data real pwindow) {
  if (is_nan(terms_d[3]) || is_nan(terms_q[3])) return not_a_number();
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
  real result = fmin(log_sum_exp(terms_q[1], log_num - log_den), 0);
  int no_tail = is_nan(terms_d[4]) || is_nan(terms_q[4]);
  real log_size = rho * d + terms_d[3] + fmax(terms_d[5], terms_q[5]);
  if (no_tail && result > negative_infinity()
      && log_size - log_den > log(1e5)) {
    return not_a_number();
  }
  return result;
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
  * Compute the exponentially tilted log CDF with shared tilt quantities
  * @ingroup exponential_tilt_solutions
  *
  * As primarycensored_exptilt_lcdf(), with the output of
  * log_tilt_transform_context() passed in so that calls at several delays
  * can share it. The context is only used by the direct form.
  *
  * @param d Delay
  * @param dist_id Distribution identifier, see primarycensored_exptilt_lcdf()
  * @param params Array of distribution parameters
  * @param pwindow Primary event window
  * @param rho Tilt, the exponential growth rate of the primary
  * @param context Output of log_tilt_transform_context() for -rho
  *
  * @return Log of the primary event censored CDF at d
  */
real primarycensored_exptilt_lcdf_shared(data real d, int dist_id,
                                         array[] real params,
                                         data real pwindow, real rho,
                                         vector context) {
  int positive = dist_has_positive_support(dist_id);
  if (positive && d <= 0) {
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
  real log_cdf = primarycensored_exptilt_lcdf_from_terms(
    primarycensored_exptilt_terms(d, dist_id, rho, params, context),
    primarycensored_exptilt_terms(q, dist_id, rho, params, context),
    d, rho, pwindow
  );
  if (is_nan(log_cdf)) {
    return log(
      primarycensored_numeric_cdf(d | dist_id, params, pwindow, 2, {rho})
    );
  }
  return log_cdf;
}

/**
  * Compute the primary event censored log CDF for an exponentially tilted
  * primary
  * @ingroup exponential_tilt_solutions
  *
  * Chooses the direct form or the small tilt forms. A zero width window
  * gives the delay log CDF. Where the direct form is `nan` it uses the ODE
  * path of primarycensored_numeric_cdf(). Only for check_for_exptilt() is 1
  * and check_for_tilt_transform() is 1 for -rho.
  *
  * @param d Delay
  * @param dist_id Distribution identifier: 1 (Lognormal), 2 (Gamma),
  *   3 (Weibull), 4 (Exponential), 5 (Generalised gamma) or 18 (Normal)
  * @param params Array of distribution parameters
  * @param pwindow Primary event window
  * @param rho Tilt, the exponential growth rate of the primary
  *
  * @return Log of the primary event censored CDF at d
  */
real primarycensored_exptilt_lcdf(data real d, int dist_id,
                                  array[] real params, data real pwindow,
                                  real rho) {
  vector[4] context = rep_vector(0, 4);
  if (primarycensored_exptilt_form(d, dist_id, rho, pwindow) == 0) {
    context = log_tilt_transform_context(dist_id, -rho, params);
  }
  return primarycensored_exptilt_lcdf_shared(
    d, dist_id, params, pwindow, rho, context
  );
}

/**
  * Compute the exponentially tilted log CDF at the two ends of an interval
  * @ingroup exponential_tilt_solutions
  *
  * As primarycensored_exptilt_lcdf() at d and d_upper, with
  * log_tilt_transform_context() computed once for both. Used by
  * primarycensored_lpmf().
  *
  * @param d Lower end of the interval
  * @param d_upper Upper end of the interval
  * @param dist_id Distribution identifier, see primarycensored_exptilt_lcdf()
  * @param params Array of distribution parameters
  * @param pwindow Primary event window
  * @param rho Tilt, the exponential growth rate of the primary
  *
  * @return Vector of the log CDF at d and at d_upper
  */
vector primarycensored_exptilt_lcdf_pair(data real d, data real d_upper,
                                         int dist_id, array[] real params,
                                         data real pwindow, real rho) {
  vector[4] context = rep_vector(0, 4);
  if (primarycensored_exptilt_form(d, dist_id, rho, pwindow) == 0
      || primarycensored_exptilt_form(d_upper, dist_id, rho, pwindow) == 0) {
    context = log_tilt_transform_context(dist_id, -rho, params);
  }
  return [
    primarycensored_exptilt_lcdf_shared(
      d, dist_id, params, pwindow, rho, context
    ),
    primarycensored_exptilt_lcdf_shared(
      d_upper, dist_id, params, pwindow, rho, context
    )
  ]';
}

/**
  * Compute the exponentially tilted log CDF at integer delays
  * @ingroup exponential_tilt_solutions
  *
  * The small tilt forms of primarycensored_exptilt_form() only hold below a
  * delay, so they cover a first run of delays, each computed with
  * primarycensored_exptilt_lcdf(). The direct form covers the rest, with the
  * terms at each endpoint computed once and shared by the delays t and
  * t + pwindow. Where the direct form uses the ODE, its values are shifted by
  * their difference from the analytical value at the delay before, so that
  * the PMF across the switch is a difference of ODE values.
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
  vector[4] context = log_tilt_transform_context(dist_id, -rho, params);
  // terms[t - d0 + pw + 1] holds the terms at endpoint t
  array[n - d0 + pw + 1] vector[5] terms;
  for (t in (d0 - pw):n) {
    terms[t - d0 + pw + 1] = primarycensored_exptilt_terms(
      t, dist_id, rho, params, context
    );
  }
  // Offset of the ODE CDF from the analytical CDF at the delay before
  real shift = 0;
  int is_ode = d0 == start;
  for (d in d0:n) {
    log_cdfs[d] = primarycensored_exptilt_lcdf_from_terms(
      terms[d - d0 + pw + 1], terms[d - d0 + 1], d, rho, pw
    );
    if (is_nan(log_cdfs[d])) {
      if (!is_ode) {
        shift = exp(log_cdfs[d - 1]) - primarycensored_numeric_cdf(
          d - 1 | dist_id, params, pw, 2, {rho}
        );
      }
      real ode_cdf = primarycensored_numeric_cdf(
        d | dist_id, params, pw, 2, {rho}
      );
      // The ODE CDF increases, so only solver error makes this negative
      log_cdfs[d] = log(shift + ode_cdf > 0 ? shift + ode_cdf : ode_cdf);
      is_ode = 1;
    } else {
      is_ode = 0;
    }
  }
  return log_cdfs;
}
