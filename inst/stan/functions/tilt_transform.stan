/*
 * Truncated exponential-moment transforms of delay distributions
 *
 * The primary event censored CDF for several non-uniform primary event
 * windows is a sum of terms built from
 *   T_f(xi; tau) = int_{-inf}^{tau} exp(xi u) f(u) du,
 * with the lower limit 0 for delays on the non-negative reals. The primary
 * event window fixes the tilts xi and the coefficients. The delay
 * distribution fixes whether T_f is closed form, and these functions are
 * the extension points for new delay distributions. A delay is added by a
 * branch in check_for_tilt_transform(), log_tilt_transform_pair_shared() and,
 * for the small tilt forms of the exponentially tilted window,
 * primarycensored_tilt_moments(). A delay whose transform has a setup that is
 * the same at every point also gives it in log_tilt_transform_context().
 * The R equivalents are the `.pcens_tilt_*()` generics.
 */

/**
  * Log of the difference of two exponentials, zero when it would be negative
  * @ingroup tilt_transforms
  *
  * log_diff_exp() is NaN if the first argument is smaller than the second
  * and its derivative is not finite when they are equal. Here the difference
  * of positive integrals that is zero to rounding is zero, and a zero
  * subtrahend returns the first argument.
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
  * The derivative of std_normal_lcdf() is an approximation with a relative
  * error of about 1e-5. Differences of terms for tilts close to zero
  * amplify that by 1e3 or more, so the normal tilt terms use log of the
  * CDF, whose derivative is the exact density over the CDF. Phi() is
  * 0.5 * (1 + erf(z / sqrt(2))) for z between -5 and 8.25, which loses
  * relative precision for negative z (1.3e-6 at z = -4.9) and the direct
  * form for a small tilt amplifies that by 1 / (|rho| w). For negative z
  * this uses 0.5 * erfc(-z / sqrt(2)), which is accurate in the tail and
  * whose derivative is exact. It underflows below -37.5. Below -37 this uses
  * the asymptotic series Phi(z) = phi(z) / (-z) *
  * (1 - 1 / z^2 + 3 / z^4 - 15 / z^6 + ... ), whose terms are
  * (-1)^k (2k - 1)!! / z^(2k). Seven terms past the first have a truncation
  * error below 1e-16 for z <= -37. Both the value and the derivative are
  * from elementary operations, so the derivative is exact, where the
  * derivative of std_normal_lcdf() has a relative error of about 1e-5. For z
  * at least 0 log(Phi(z)) is accurate and saturates at 0 above 8.25.
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
  * Log of the regularised lower incomplete gamma function from its series
  * @ingroup tilt_transforms
  *
  * The shape derivative of gamma_lcdf() is inaccurate well below the
  * shape, with a relative error of 1.7e-2 for shape 20 at 2 and of 0.5 for
  * shape 100 at 30, where the value is correct. It is also NaN for shapes
  * of about 200 or more in the tail and of about 600 or more over the bulk.
  * This uses the lower series
  * P(shape, x) = x^shape exp(-x) / Gamma(shape + 1) *
  *   sum_k x^k / ((shape + 1) ... (shape + k)),
  * built from elementary operations so that autodiff is exact. It is for
  * x < shape + 1, where each term is less than the previous one. Near
  * x = shape it needs about sqrt(80 shape) terms to converge to double
  * precision, and far fewer for x well below the shape. It rejects if it
  * has not converged after 10 sqrt(shape) + 150 terms, which does not happen
  * for x < shape + 1.
  *
  * @param x Point, positive, below shape + 1
  * @param shape Shape, positive
  *
  * @return log P(shape, x)
  */
real primarycensored_log_gamma_p_series(real x, real shape) {
  real term = 1;
  real total = 1;
  real max_terms = 10 * sqrt(shape) + 150;
  int k = 0;
  while (term >= 1e-17 * total) {
    k += 1;
    if (k > max_terms) {
      reject("The gamma lower series did not converge for x = ", x,
             " and shape = ", shape);
    }
    term *= x / (shape + k);
    total += term;
  }
  return shape * log(x) - x - lgamma(shape + 1) + log(total);
}

/**
  * Log of the regularised upper incomplete gamma function from its
  * continued fraction
  * @ingroup tilt_transforms
  *
  * The shape derivative of gamma_lccdf() is inaccurate in the bulk, with
  * relative errors of 1e-3 to 1e-2 (for example shape 2.5 at 7 or shape 20 at
  * 24), and is NaN in the tail for shapes of about 200 or more. This uses
  * the Legendre continued fraction
  * Q(shape, x) = x^shape exp(-x) / Gamma(shape) *
  *   1 / (x + 1 - shape - 1 (1 - shape) / (x + 3 - shape -
  *   2 (2 - shape) / (x + 5 - shape - ... ))),
  * evaluated by the modified Lentz method from elementary operations, so
  * that autodiff is exact and the tail is held on the log scale to any
  * depth. It is for x >= shape + 1, where it converges in at most about
  * 90 terms for shapes up to 1e3 and about sqrt(shape) terms beyond. It
  * rejects if it has not converged after 10 sqrt(shape) + 150 terms, which
  * does not happen for x >= shape + 1.
  *
  * @param x Point, positive, at least shape + 1
  * @param shape Shape, positive
  *
  * @return log Q(shape, x)
  */
real primarycensored_log_gamma_q_fraction(real x, real shape) {
  real b = x + 1 - shape;
  real d = 1 / b;
  real c = 1e300;
  real h = d;
  real max_terms = 10 * sqrt(shape) + 150;
  real del = 0;
  int i = 0;
  while (abs(del - 1) > 1e-15) {
    i += 1;
    if (i > max_terms) {
      reject("The gamma upper fraction did not converge for x = ", x,
             " and shape = ", shape);
    }
    real an = -i * (i - shape);
    b += 2;
    d = an * d + b;
    if (abs(d) < 1e-300) d = 1e-300;
    c = b + an / c;
    if (abs(c) < 1e-300) c = 1e-300;
    d = 1 / d;
    del = d * c;
    h *= del;
  }
  return shape * log(x) - x - lgamma(shape) + log(h);
}

/**
  * Log of the regularised lower and upper incomplete gamma functions
  * @ingroup tilt_transforms
  *
  * The tail that is not close to 1 is evaluated directly, from the series
  * for x < shape + 1 and from the continued fraction beyond, see
  * primarycensored_log_gamma_p_series() and
  * primarycensored_log_gamma_q_fraction(). The other tail is log1m_exp() of
  * it, which is accurate in value and derivative, including for an upper
  * tail far below the smallest double, where it is 0 with a zero
  * derivative. The derivatives in the shape and in x are exact in both
  * tails for any shape.
  *
  * @param x Point, positive
  * @param shape Shape, positive
  *
  * @return Vector [log P(shape, x), log Q(shape, x)]
  */
vector primarycensored_log_gamma_pq(real x, real shape) {
  if (x < shape + 1) {
    real log_lower = primarycensored_log_gamma_p_series(x, shape);
    return [log_lower, log1m_exp(log_lower)]';
  }
  real log_upper = primarycensored_log_gamma_q_fraction(x, shape);
  return [log1m_exp(log_upper), log_upper]';
}

/**
  * Log of the regularised lower incomplete gamma function
  * @ingroup tilt_transforms
  *
  * See primarycensored_log_gamma_pq(), which gives both tails.
  *
  * @param x Point, positive
  * @param shape Shape, positive
  *
  * @return log P(shape, x)
  */
real primarycensored_log_gamma_p(real x, real shape) {
  if (x < shape + 1) return primarycensored_log_gamma_p_series(x, shape);
  return log1m_exp(primarycensored_log_gamma_q_fraction(x, shape));
}

/**
  * Test whether the gamma lower tail underflows at these arguments
  * @ingroup tilt_transforms
  *
  * The tails are evaluated on the log scale, see
  * primarycensored_log_gamma_pq(), so they are represented below the smallest
  * double. A term below 1e-300 of a probability is dropped as `-inf`
  * instead, which saves its evaluation. For x < shape + 1 this uses the bound
  * P(shape, x) <= x^shape exp(-x) / (Gamma(shape + 1) (1 - x / (shape + 1)))
  * with a margin below the smallest double. The transforms of a gamma with a
  * large total can be representable below that, and are then dropped with it.
  *
  * @param x Point, positive
  * @param shape Shape
  *
  * @return 1 if the regularised lower incomplete gamma function may underflow,
  *   0 otherwise
  */
int gamma_lcdf_underflows(real x, real shape) {
  if (x >= shape + 1) return 0;
  return shape * log(x) - lgamma(shape + 1) - x - log1m(x / (shape + 1))
         < -700;
}

/**
  * Test whether the gamma upper tail underflows at these arguments
  * @ingroup tilt_transforms
  *
  * As for gamma_lcdf_underflows(), for the upper tail. It uses the asymptotic
  * form Q(shape, x) = x^(shape - 1) exp(-x) / Gamma(shape) x / (x - shape + 1)
  * for x beyond shape + 1, where the tail is not close to 1.
  *
  * @param x Point, positive
  * @param shape Shape
  *
  * @return 1 if the regularised upper incomplete gamma function may underflow,
  *   0 otherwise
  */
int gamma_lccdf_underflows(real x, real shape) {
  if (x <= shape + 1) return 0;
  return (shape - 1) * log(x) - x - lgamma(shape) + log(x / (x - shape + 1))
         < -700;
}

/**
  * Check if the tilt transform is closed form for a delay and tilt
  * @ingroup tilt_transforms
  *
  * The exponential (4) and gamma (2) forms are the total times the CDF of the
  * tilted delay, a delay with the rate lowered by xi. That delay exists only
  * if rate - xi > 0. The normal (18) form has no restriction. The lognormal
  * (1) form needs xi sigma^2 exp(mu) to stay below exp(690) in magnitude.
  * Callers use the numerical path when this is 0. The ODE is less accurate
  * there for the lower tail of a gamma with a shape below 1 (a relative
  * error of about 4e-2 for shape 0.3 at 1e-3 and tilt -1 with rate 1).
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
  * xi > 0, that the series is the faster accurate path at t. The series
  * costs about xi t + 9 sqrt(xi t) + 30 terms. It is used up to xi t of 60,
  * where it matches an ODE at tolerance 1e-10 in speed, and beyond that only
  * where xi w is above 2, as the ODE is inaccurate in the lower tail there.
  * The series is limited to 20000 terms.
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
  * The lower transform is T_f(xi; t). For xi = 0 it is the log CDF of the
  * delay. For delays on the non-negative reals it is `-inf` for t <= 0. The
  * upper transform is T_f(xi; Inf) - T_f(xi; t), the log of the total for
  * t <= 0 for those delays. Evaluating the tail directly keeps precision
  * where the lower transform is close to its total, see
  * primarycensored_tail_diff(). Only defined where check_for_tilt_transform()
  * is 1.
  *
  * For the gamma the tilted density is a gamma density with the rate lowered
  * by xi, times the total (rate / (rate - xi))^shape. Both tails come from
  * primarycensored_log_gamma_pq(), which evaluates one from a series or a
  * continued fraction and the other as log1m_exp() of it. This halves the
  * cost of the incomplete gamma function and of its derivative in the shape,
  * and the derivatives are exact for any shape. Stan's gamma_lcdf() and
  * gamma_lccdf() are not used, as their shape derivatives are inaccurate in
  * parts of the bulk and NaN for shapes of about 200 or more. For the
  * exponential T_f = rate / (rate - xi) (1 - exp(-(rate - xi) t)). For the
  * normal, completing the square gives a normal density with mean
  * mu + xi sigma^2.
  *
  * The log of the total is -shape log1m(xi / rate) (for the exponential,
  * the shape is 1), rather than shape (log(rate) - log(rate - xi)). The
  * derivative in the rate is then one term, which is exactly 0 for xi = 0.
  * Two terms of size shape / rate that cancel would swamp the derivative of
  * a tail that is far below 1, and the derivative in the rate of a log CDF
  * close to 0 would be lost, for example 0 instead of 2e-18 for shape 200,
  * rate 20 at 17.
  *
  * @param t Point
  * @param dist_id Distribution identifier: 1 (Lognormal), 2 (Gamma),
  *   4 (Exponential) or 18 (Normal), see check_for_tilt_transform()
  * @param xi Tilt
  * @param params Array of distribution parameters, as for dist_lcdf()
  * @param context Output of log_tilt_transform_context() for dist_id, xi
  *   and params
  *
  * @return Vector [log T_f(xi; t), log(T_f(xi; Inf) - T_f(xi; t))]
  */
vector log_tilt_transform_pair_shared(real t, int dist_id, real xi,
                                      array[] real params, vector context) {
  if (dist_id == 2) {
    real shape = params[1];
    real rate = params[2];
    real tilted_rate = rate - xi;
    real log_total = -shape * log1m(xi / rate);
    if (t <= 0) return [negative_infinity(), log_total]';
    real x = t * tilted_rate;
    if (gamma_lcdf_underflows(x, shape)) {
      return [negative_infinity(), log_total]';
    }
    // Beyond the point where the upper tail underflows the CDF is 1
    if (gamma_lccdf_underflows(x, shape)) {
      return [log_total, negative_infinity()]';
    }
    vector[2] log_tails = primarycensored_log_gamma_pq(x, shape);
    return [log_total + log_tails[1], log_total + log_tails[2]]';
  } else if (dist_id == 4) {
    real rate = params[1];
    real tilted_rate = rate - xi;
    real log_total = -log1m(xi / rate);
    if (t <= 0) return [negative_infinity(), log_total]';
    return [
      log_total + log1m_exp(-tilted_rate * t), log_total - tilted_rate * t
    ]';
  } else if (dist_id == 18) {
    // The upper tail is the lower tail of the reflected normal, as
    // normal_lccdf() is -inf beyond 8.25 standard deviations, where the
    // standard normal upper tail is still representable
    real mu = params[1];
    real sigma = params[2];
    real z = (t - mu - xi * square(sigma)) / sigma;
    real log_total = xi * mu + 0.5 * square(xi * sigma);
    return [
      log_total + primarycensored_log_std_normal_cdf(z),
      log_total + primarycensored_log_std_normal_cdf(-z)
    ]';
  } else if (dist_id == 1) {
    return primarycensored_lognormal_tilt_pair(
      t, params[1], params[2], xi, context
    );
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
  * @param dist_id Distribution identifier: 1 (Lognormal), 2 (Gamma),
  *   4 (Exponential) or 18 (Normal), see check_for_tilt_transform()
  * @param xi Tilt
  * @param params Array of distribution parameters, as for dist_lcdf()
  *
  * @return Vector [log T_f(xi; t), log(T_f(xi; Inf) - T_f(xi; t))]
  */
vector log_tilt_transform_pair(real t, int dist_id, real xi,
                               array[] real params) {
  return log_tilt_transform_pair_shared(
    t, dist_id, xi, params, log_tilt_transform_context(dist_id, xi, params)
  );
}

/**
  * Log of the tilt transform over the lower part of the support
  * @ingroup tilt_transforms
  *
  * For xi = 0 this is the log CDF of the delay. For delays on the
  * non-negative reals it is `-inf` for t <= 0. See
  * log_tilt_transform_pair(), which also gives the upper part.
  *
  * @param t Upper limit of the transform
  * @param dist_id Distribution identifier: 1 (Lognormal), 2 (Gamma),
  *   4 (Exponential) or 18 (Normal), see check_for_tilt_transform()
  * @param xi Tilt
  * @param params Array of distribution parameters, as for dist_lcdf()
  *
  * @return log T_f(xi; t)
  */
real log_tilt_transform(real t, int dist_id, real xi, array[] real params) {
  return log_tilt_transform_pair(t, dist_id, xi, params)[1];
}

/**
  * Log of the tilt transform over the upper part of the support
  * @ingroup tilt_transforms
  *
  * The log of T_f(xi; Inf) - T_f(xi; t). Evaluating the tail directly keeps
  * precision where log_tilt_transform() is close to its total. For delays on
  * the non-negative reals it is the log of the total for t <= 0. Only
  * defined where check_for_tilt_transform() is 1. See
  * log_tilt_transform_pair().
  *
  * @param t Lower limit of the transform
  * @param dist_id Distribution identifier: 1 (Lognormal), 2 (Gamma),
  *   4 (Exponential) or 18 (Normal), see check_for_tilt_transform()
  * @param xi Tilt
  * @param params Array of distribution parameters, as for dist_lcdf()
  *
  * @return log(T_f(xi; Inf) - T_f(xi; t))
  */
real log_tilt_transform_upper(real t, int dist_id, real xi,
                              array[] real params) {
  return log_tilt_transform_pair(t, dist_id, xi, params)[2];
}

/**
  * Log moments of a gamma delay about a point
  * @ingroup tilt_transforms
  *
  * For t > 0 the log of G_1(t) = int_0^t (t - u) f(u) du and
  * G_2(t) = int_0^t (t - u)^2 f(u) du. They come from the CDFs of gamma
  * distributions with the shape raised by one and two, which give the partial
  * moments of the delay, and every difference is of positive integrals.
  *
  * @param t Point, positive
  * @param shape Shape
  * @param rate Rate
  *
  * @return Vector [log G_1(t), log G_2(t)]
  */
vector primarycensored_gamma_tilt_moments(real t, real shape, real rate) {
  // The moments underflow with the CDFs of the raised shapes
  if (gamma_lcdf_underflows(t * rate, shape)
      || gamma_lcdf_underflows(t * rate, shape + 1)
      || gamma_lcdf_underflows(t * rate, shape + 2)) {
    return rep_vector(negative_infinity(), 2);
  }
  // Beyond the point where the upper tails underflow the CDFs are 1 and the
  // moments are those of the whole distribution
  if (gamma_lccdf_underflows(t * rate, shape + 2)) {
    real mean_delay = shape / rate;
    real second_moment = shape * (shape + 1) / square(rate);
    return [
      log(t - mean_delay),
      log(square(t) - 2 * t * mean_delay + second_moment)
    ]';
  }
  real log_t = log(t);
  real x = t * rate;
  real log_m0 = primarycensored_log_gamma_p(x, shape);
  // Partial first and second moments of the delay
  real log_m1 = log(shape) - log(rate)
                + primarycensored_log_gamma_p(x, shape + 1);
  real log_m2 = log(shape) + log(shape + 1) - 2 * log(rate)
                + primarycensored_log_gamma_p(x, shape + 2);
  real log_g1 = primarycensored_log_diff_exp(log_t + log_m0, log_m1);
  real log_h = primarycensored_log_diff_exp(log_t + log_m1, log_m2);
  real log_g2 = primarycensored_log_diff_exp(log_t + log_g1, log_h);
  return [log_g1, log_g2]';
}

/**
  * Log moments of a delay about a point
  * @ingroup tilt_transforms
  *
  * The log of G_1(t) = int (t - u) f(u) du and
  * G_2(t) = int (t - u)^2 f(u) du over the support up to t. They are used
  * by the small tilt forms of the exponentially tilted primary, which would
  * otherwise cancel as the tilt goes to zero. `-inf` for both for t <= 0
  * for delays on the non-negative reals.
  *
  * The exponential is the gamma with shape 1. Its own closed forms cancel
  * when the rate times t is small. The normal with z = (t - mu) / sigma has
  * G_1 = sigma (phi(z) + z Phi(z)) and
  * G_2 = sigma^2 ((z^2 + 1) Phi(z) + z phi(z)). The lognormal uses its
  * partial moments m_k(t), see primarycensored_lognormal_tilt_series().
  *
  * @param t Point
  * @param dist_id Distribution identifier: 1 (Lognormal), 2 (Gamma),
  *   4 (Exponential) or 18 (Normal), see check_for_tilt_transform()
  * @param params Array of distribution parameters, as for dist_lcdf()
  *
  * @return Vector [log G_1(t), log G_2(t)]
  */
vector primarycensored_tilt_moments(real t, int dist_id,
                                    array[] real params) {
  if (dist_id == 2) {
    if (t <= 0) return rep_vector(negative_infinity(), 2);
    return primarycensored_gamma_tilt_moments(t, params[1], params[2]);
  } else if (dist_id == 4) {
    if (t <= 0) return rep_vector(negative_infinity(), 2);
    return primarycensored_gamma_tilt_moments(t, 1, params[1]);
  } else if (dist_id == 18) {
    real mu = params[1];
    real sigma = params[2];
    real z = (t - mu) / sigma;
    real log_phi = std_normal_lpdf(z);
    real log_Phi = primarycensored_log_std_normal_cdf(z);
    real log_g1;
    real log_g2;
    if (z > -1) {
      // The terms do not cancel badly here, and the direct form has a
      // derivative at z = 0, where the log form does not
      log_g1 = log(exp(log_phi) + z * exp(log_Phi));
      log_g2 = log((square(z) + 1) * exp(log_Phi) + z * exp(log_phi));
    } else {
      // For z < 0 the terms of each sum have opposite signs, the first
      // larger, and the log form keeps the tail precise
      log_g1 = primarycensored_log_diff_exp(log_phi, log(-z) + log_Phi);
      log_g2 = primarycensored_log_diff_exp(
        log1p(square(z)) + log_Phi, log(-z) + log_phi
      );
    }
    return [log(sigma) + log_g1, 2 * log(sigma) + log_g2]';
  } else if (dist_id == 1) {
    if (t <= 0) return rep_vector(negative_infinity(), 2);
    real mu = params[1];
    real sigma = params[2];
    real log_t = log(t);
    real z = (log_t - mu) / sigma;
    real log_m0 = primarycensored_log_std_normal_cdf(z);
    real log_m1 = mu + 0.5 * square(sigma)
                  + primarycensored_log_std_normal_cdf(z - sigma);
    real log_m2 = 2 * mu + 2 * square(sigma)
                  + primarycensored_log_std_normal_cdf(z - 2 * sigma);
    real log_g1 = primarycensored_log_diff_exp(log_t + log_m0, log_m1);
    real log_h = primarycensored_log_diff_exp(log_t + log_m1, log_m2);
    real log_g2 = primarycensored_log_diff_exp(log_t + log_g1, log_h);
    return [log_g1, log_g2]';
  }
  reject("Invalid distribution identifier: ", dist_id);
}
