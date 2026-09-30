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
 * branch in check_for_tilt_transform(), log_tilt_transform_pair() and, for
 * the small tilt forms of the exponentially tilted window,
 * primarycensored_tilt_moments().
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
  * whose derivative is exact. It underflows below -37.5, where
  * std_normal_lcdf() is used. For z at least 0 log(Phi(z)) is accurate and
  * saturates at 0 above 8.25.
  *
  * @param z Point
  *
  * @return log(Phi(z))
  */
real primarycensored_log_std_normal_cdf(real z) {
  if (z >= 0) return log(Phi(z));
  if (z > -37) return log(0.5 * erfc(-z / sqrt(2)));
  return std_normal_lcdf(z);
}

/**
  * Log of the regularised lower incomplete gamma function
  * @ingroup tilt_transforms
  *
  * The shape derivative of gamma_lcdf() is inaccurate well below the
  * shape, with a relative error of 1.7e-2 for shape 20 at 2 and of 0.5 for
  * shape 100 at 30, where the value is correct. Its accuracy recovers above
  * about 0.4 of the shape, to 1e-8 or better. Below x = 0.5 shape this uses
  * the lower series
  * P(shape, x) = x^shape exp(-x) / Gamma(shape + 1) *
  *   sum_k x^k / ((shape + 1) ... (shape + k)),
  * built from elementary operations so that autodiff is exact. Each term is
  * less than half the previous one, so it converges to double precision in
  * at most about 55 terms, and in fewer for x well below the shape.
  *
  * @param x Point, positive
  * @param shape Shape, positive
  *
  * @return log P(shape, x)
  */
real primarycensored_log_gamma_p(real x, real shape) {
  if (x >= 0.5 * shape) return gamma_lcdf(x | shape, 1);
  real term = 1;
  real total = 1;
  for (k in 1:200) {
    term *= x / (shape + k);
    total += term;
    if (term < 1e-17 * total) break;
  }
  return shape * log(x) - x - lgamma(shape + 1) + log(total);
}

/**
  * Test whether the gamma lower tail underflows at these arguments
  * @ingroup tilt_transforms
  *
  * Underflow of `gamma_p` makes `gamma_lcdf` `-inf` with an autodiff partial
  * that is not finite, which Stan's reverse pass chains into the parameters.
  * Callers test this before calling `gamma_lcdf`, rather than checking the
  * result. For x < shape + 1 it uses the bound
  * P(shape, x) <= x^shape exp(-x) / (Gamma(shape + 1) (1 - x / (shape + 1)))
  * with a margin below the smallest double, so a term dropped on this test is
  * below 1e-300 of a probability. The transforms of a gamma with a large total
  * can be representable below that, and are then dropped with it.
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
  * As for gamma_lcdf_underflows(), for `gamma_lccdf`. It uses the asymptotic
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
  * if rate - xi > 0. The normal (18) form has no restriction. Callers use the
  * numerical path when this is 0. The ODE is less accurate there for the
  * lower tail of a gamma with a shape below 1 (a relative error of about
  * 4e-2 for shape 0.3 at 1e-3 and tilt -1 with rate 1).
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
  * The lower transform is T_f(xi; t). For xi = 0 it is the log CDF of the
  * delay. For delays on the non-negative reals it is `-inf` for t <= 0. The
  * upper transform is T_f(xi; Inf) - T_f(xi; t), the log of the total for
  * t <= 0 for those delays. Evaluating the tail directly keeps precision
  * where the lower transform is close to its total, see
  * primarycensored_tail_diff(). Only defined where check_for_tilt_transform()
  * is 1.
  *
  * For the gamma the tilted density is a gamma density with the rate lowered
  * by xi, times the total (rate / (rate - xi))^shape. One tail is evaluated
  * and the other follows from it, which halves the cost of the incomplete
  * gamma function and of its derivative in the shape. The lower tail is
  * primarycensored_log_gamma_p(), which is gamma_lcdf() except well below
  * the shape, where the shape derivative of gamma_lcdf() is inaccurate. The
  * upper tail is log1m_exp() of it,
  * which is accurate in value and derivative for an upper tail down to
  * about 1e-9. The shape derivative of gamma_lccdf() is worse, with
  * relative errors of 1e-3 to 1e-2 over much of the bulk (for example shape
  * 2.5 at 7 or shape 20 at 24) and above 1e-4 into the far tail for large
  * shapes. It is used only where the upper tail is below 1e-8, where
  * log1m_exp() of the lower tail loses the tail or is `-inf` with a
  * derivative that is not finite. For the exponential
  * T_f = rate / (rate - xi) (1 - exp(-(rate - xi) t)). For the normal,
  * completing the square gives a normal density with mean mu + xi sigma^2.
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
    real tilted_rate = rate - xi;
    real log_total = shape * (log(rate) - log(tilted_rate));
    if (t <= 0) return [negative_infinity(), log_total]';
    real x = t * tilted_rate;
    if (gamma_lcdf_underflows(x, shape)) {
      return [negative_infinity(), log_total]';
    }
    // Beyond the point where the upper tail underflows the CDF is 1
    if (gamma_lccdf_underflows(x, shape)) {
      return [log_total, negative_infinity()]';
    }
    real log_lower = primarycensored_log_gamma_p(x, shape);
    if (log_lower > -1e-8) {
      real log_upper = gamma_lccdf(x | shape, 1);
      return [log_total + log1m_exp(log_upper), log_total + log_upper]';
    }
    return [log_total + log_lower, log_total + log1m_exp(log_lower)]';
  } else if (dist_id == 4) {
    real rate = params[1];
    real tilted_rate = rate - xi;
    real log_total = log(rate) - log(tilted_rate);
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
  }
  reject("Invalid distribution identifier: ", dist_id);
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
  * @param dist_id Distribution identifier: 2 (Gamma), 4 (Exponential) or 18
  *   (Normal), see check_for_tilt_transform()
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
  * @param dist_id Distribution identifier: 2 (Gamma), 4 (Exponential) or 18
  *   (Normal), see check_for_tilt_transform()
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
  * G_2 = sigma^2 ((z^2 + 1) Phi(z) + z phi(z)).
  *
  * @param t Point
  * @param dist_id Distribution identifier: 2 (Gamma), 4 (Exponential) or 18
  *   (Normal), see check_for_tilt_transform()
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
  }
  reject("Invalid distribution identifier: ", dist_id);
}
