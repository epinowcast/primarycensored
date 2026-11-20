#ifndef PRIMARYCENSORED_EXPTILT_HPP
#define PRIMARYCENSORED_EXPTILT_HPP

/**
 * C++ version, with closed-form gradients, of the primarycensored Stan
 * function primarycensored_tilt_pair(): the log tilt transforms
 * log T_f(xi; t) and log(T_f(xi; Inf) - T_f(xi; t)) of gamma, Weibull,
 * exponential, generalised gamma and normal delays, used by the
 * exponentially tilted primary.
 *
 * Its Stan body, in inst/stan/functions/primarycensored_exptilt.stan of the
 * primarycensored package, is the reference. Values match the Stan code.
 *
 * Exponential (rate lambda), with r = lambda - xi:
 *   log T = log(lambda / r) + log(1 - exp(-r t)),
 *   log U = log(lambda / r) - r t.
 * Gamma (k, beta), with r = beta - xi and y = r t:
 *   log T = k log(beta / r) + log P(k, y),
 *   log U = k log(beta / r) + log Q(k, y).
 * The log y partials of log P and log Q are y f(y) / P and -y f(y) / Q,
 * with f the unit rate gamma density. The k partials differentiate the
 * same series and continued fraction as the Stan code, so they match
 * Stan's autodiff through it.
 * Normal (mu, sigma), with z = (t - mu - xi sigma^2) / sigma:
 *   log T = xi mu + (xi sigma)^2 / 2 + log Phi(z),
 *   log U = xi mu + (xi sigma)^2 / 2 + log Phi(-z).
 * Weibull and generalised gamma (shape, scale, k), x = (t / scale)^shape:
 * for xi = 0 the CDF P(k, x) and its complement, as for the gamma. For
 * xi other than 0 the series of primarycensored_stacy_series(), whose
 * terms are differentiated one by one, see stacy_tilt_pair().
 * Lognormal (mu, sigma), for xi >= 0 only: the CDF for xi = 0 and the
 * partial moment series for xi > 0, see lognormal_tilt_pair(). The
 * quadrature for xi < 0 stays in Stan.
 * The third element, the log loss of the series, is 0 for the others.
 * Each call is one autodiff node.
 *
 * The function is put in stan::math, as in uniform_lcdf.hpp.
 */

#include "uniform_lcdf.hpp"
#include <stan/math.hpp>
#include <array>
#include <cmath>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace primarycensored {
namespace internal {

/**
 * Log CDF of a unit rate gamma from log(x), as gamma_lcdf_logx() in Stan,
 * with its derivative in the shape as Stan's autodiff finds it.
 *
 * The series is differentiated term by term. Elsewhere the derivative is
 * that of Stan's gamma_lcdf(), by nested autodiff.
 *
 * @param log_x Log of the argument.
 * @param a Shape, positive and finite.
 * @param dlog_p Set to d log P / da if not null.
 * @return log P(a, exp(log_x)).
 */
inline double gamma_lcdf_logx_da(double log_x, double a, double* dlog_p) {
  if (!dlog_p) {
    return gamma_lcdf_logx(log_x, a);
  }
  if (!(a > 0) || std::isinf(a)) {
    throw std::domain_error(
        "gamma_lcdf_logx: shape must be positive finite");
  }
  const double x = std::exp(log_x);
  if (!(x < 0.5 * (a + 1))) {
    stan::math::nested_rev_autodiff nested;
    stan::math::var a_v = a;
    stan::math::var log_p = stan::math::gamma_lcdf(x, a_v, 1.0);
    log_p.grad();
    *dlog_p = a_v.adj();
    return log_p.val();
  }
  double term = 1;
  double total = 1;
  double dterm = 0;
  double dtotal = 0;
  double n = a;
  while (term > 1e-16 * total) {
    n += 1;
    const double ratio = x / n;
    dterm = dterm * ratio - term * (ratio / n);
    term *= ratio;
    total += term;
    dtotal += dterm;
  }
  *dlog_p = log_x - stan::math::digamma(a + 1) + dtotal / total;
  return a * log_x - x - stan::math::lgamma(a + 1) + std::log(total);
}

/**
 * Log CDF and log CCDF of a unit rate gamma from log(x), as
 * gamma_lcdf_lccdf_logx() in Stan, with their derivatives in the shape as
 * Stan's autodiff finds them. The continued fraction is differentiated
 * step by step.
 *
 * @param log_x Log of the argument.
 * @param a Shape, positive and finite.
 * @param da If not null, set to {d log P / da, d log Q / da}.
 * @return {log P(a, x), log Q(a, x)}.
 */
inline std::array<double, 2> gamma_lcdf_lccdf_logx(double log_x, double a,
                                                   std::array<double, 2>* da) {
  using stan::math::log1m_exp;
  if (!(a > 0) || std::isinf(a)) {
    throw std::domain_error(
        "gamma_lcdf_lccdf_logx: shape must be positive finite");
  }
  const double x = std::exp(log_x);
  if (!(x >= a + 1)) {
    double dlog_p = 0;
    const double log_p
        = gamma_lcdf_logx_da(log_x, a, da ? &dlog_p : nullptr);
    if (da) {
      // d log1m_exp(u) / du = -1 / expm1(-u)
      (*da)[0] = dlog_p;
      (*da)[1] = -dlog_p / std::expm1(-log_p);
    }
    return {log_p, log1m_exp(log_p)};
  }
  if (std::isinf(x)) {
    if (da) {
      *da = {0, 0};
    }
    return {0, neg_inf};
  }
  // Each quantity and its derivative in a
  double b = x + 1 - a, db = -1;
  double c = 1e300, dc = 0;
  double d = 1 / b, dd = -db / (b * b);
  double h = d, dh = dd;
  double delta = 1, ddelta = 0;
  double err = 1;
  int i = 0;
  while (err >= 1e-15) {
    i += 1;
    const double an = -i * (i - a);
    const double dan = i;
    const double d_prev = d;
    b += 2;
    dd = dan * d + an * dd + db;
    d = an * d + b;
    if (std::fabs(d) < 1e-300) {
      d = 1e-300;
      dd = 0;
    }
    dc = db + (dan * c - an * dc) / (c * c);
    c = b + an / c;
    if (std::fabs(c) < 1e-300) {
      c = 1e-300;
      dc = 0;
    }
    d = 1 / d;
    dd = -dd * (d * d);
    err *= i * std::fmax(std::fabs(i - a), 1) * std::fabs(d * d_prev / delta);
    ddelta = dd * c + d * dc;
    delta = d * c;
    dh = dh * delta + h * ddelta;
    h *= delta;
  }
  const double log_q = a * log_x - x - stan::math::lgamma(a) + std::log(h);
  if (da) {
    (*da)[1] = log_x - stan::math::digamma(a) + dh / h;
    (*da)[0] = -(*da)[1] / std::expm1(-log_q);
  }
  return {log1m_exp(log_q), log_q};
}

/**
 * Log of the standard normal CDF, as primarycensored_log_std_normal_cdf()
 * in Stan.
 *
 * @param z Point.
 * @return log Phi(z).
 */
inline double log_std_normal_cdf(double z) {
  if (z >= 0) {
    return std::log(stan::math::Phi(z));
  }
  if (z > -37) {
    return std::log(0.5 * stan::math::erfc(-z / std::sqrt(2.0)));
  }
  const double inv_z2 = 1 / (z * z);
  // Nested as in the Stan code, so that the value matches
  double series = 1 - 13 * inv_z2;
  for (int j = 11; j >= 1; j -= 2) {
    series = 1 - j * inv_z2 * series;
  }
  return stan::math::std_normal_lpdf(z) - std::log(-z) + std::log(series);
}

/**
 * d log Phi(z) / dz = phi(z) / Phi(z), as Stan's autodiff through
 * primarycensored_log_std_normal_cdf() finds it in each of its regimes.
 *
 * Below -37 it differentiates the asymptotic series term by term, as
 * exp(log phi - log Phi) there loses about 13 digits to the size of the
 * logs, which the normal tilt partials then amplify.
 *
 * @param z Point.
 * @param log_Phi log Phi(z), unused.
 * @return The derivative.
 */
inline double dlog_std_normal_cdf(double z, double /* log_Phi */) {
  if (z >= 0) {
    return stan::math::INV_SQRT_TWO_PI * std::exp(-0.5 * z * z)
           / stan::math::Phi(z);
  }
  if (z > -37) {
    // d log erfc(x) / dz with x = -z / sqrt(2)
    const double x = -z / std::sqrt(2.0);
    return stan::math::SQRT_TWO_OVER_SQRT_PI * std::exp(-x * x)
           / stan::math::erfc(x);
  }
  // log Phi = log phi(z) - log(-z) + log(series(u)), u = 1 / z^2, with
  // series(u) = 1 - u (1 - 3 u (1 - 5 u (...)))
  const double u = 1 / (z * z);
  double series = 1 - 13 * u;
  double dseries = -13;
  for (int j = 11; j >= 1; j -= 2) {
    dseries = -j * (series + u * dseries);
    series = 1 - j * u * series;
  }
  // du / dz = -2 / z^3
  return -z - 1 / z + dseries * (-2 * u / z) / series;
}

/**
 * Values of the three outputs and their partials with respect to xi, up to
 * three delay parameters and the point t.
 */
struct TiltPair {
  std::array<double, 3> val{{neg_inf, neg_inf, 0}};
  // grad[i][0] is d val[i] / d xi, grad[i][j] is d val[i] / d params[j - 1]
  // for j = 1, 2, 3, and grad[i][4] is d val[i] / d t
  std::array<std::array<double, 5>, 3> grad{};
};

/// Column of the partials in t, see TiltPair.
constexpr int t_col = 4;

/**
 * Exponential tilt pair, as the dist_id 4 branch in Stan.
 *
 * @param t Point.
 * @param xi Tilt.
 * @param rate Rate.
 * @param grad Whether to compute the partials.
 * @return The pair.
 */
inline TiltPair exponential_tilt_pair(double t, double xi, double rate,
                                      bool grad) {
  TiltPair out;
  const double tilted_rate = rate - xi;
  const double log_total = -stan::math::log1m(xi / rate);
  // d log_total / d rate and d xi
  const double tot_rate = 1 / rate - 1 / tilted_rate;
  const double tot_xi = 1 / tilted_rate;
  if (t <= 0) {
    out.val[1] = log_total;
    if (grad) {
      out.grad[1][0] = tot_xi;
      out.grad[1][1] = tot_rate;
    }
    return out;
  }
  out.val[0] = log_total + stan::math::log1m_exp(-tilted_rate * t);
  out.val[1] = log_total - tilted_rate * t;
  if (grad) {
    // d log1m_exp(-r t) / dr
    const double dlow = t / std::expm1(tilted_rate * t);
    out.grad[0][0] = tot_xi - dlow;
    out.grad[0][1] = tot_rate + dlow;
    out.grad[1][0] = tot_xi + t;
    out.grad[1][1] = tot_rate - t;
    out.grad[0][t_col] = tilted_rate / std::expm1(tilted_rate * t);
    out.grad[1][t_col] = -tilted_rate;
  }
  return out;
}

/**
 * Gamma tilt pair, as the dist_id 2 branch in Stan.
 *
 * @param t Point.
 * @param xi Tilt.
 * @param shape Shape.
 * @param rate Rate.
 * @param grad Whether to compute the partials.
 * @return The pair.
 */
inline TiltPair gamma_tilt_pair(double t, double xi, double shape,
                                double rate, bool grad) {
  TiltPair out;
  const double log1m_ratio = stan::math::log1m(xi / rate);
  const double log_total = -shape * log1m_ratio;
  const double tilted_rate = rate - xi;
  if (t <= 0) {
    out.val[1] = log_total;
    if (grad) {
      out.grad[1][0] = shape / tilted_rate;
      out.grad[1][1] = -log1m_ratio;
      out.grad[1][2] = shape * (1 / rate - 1 / tilted_rate);
    }
    return out;
  }
  const double log_y = std::log(t) + std::log(tilted_rate);
  std::array<double, 2> pq_da;
  const std::array<double, 2> pq
      = gamma_lcdf_lccdf_logx(log_y, shape, grad ? &pq_da : nullptr);
  out.val[0] = log_total + pq[0];
  out.val[1] = log_total + pq[1];
  if (!grad) {
    return out;
  }
  // log of y f(y), f the unit rate gamma density
  const double log_yf
      = shape * log_y - std::exp(log_y) - stan::math::lgamma(shape);
  // d log P / d log y and d log Q / d log y
  const double dlog_y[2]
      = {std::exp(log_yf - pq[0]), -std::exp(log_yf - pq[1])};
  const double tot_xi = shape / tilted_rate;
  const double tot_rate = shape * (1 / rate - 1 / tilted_rate);
  for (int i = 0; i < 2; ++i) {
    out.grad[i][0] = tot_xi - dlog_y[i] / tilted_rate;
    out.grad[i][1] = -log1m_ratio + pq_da[i];
    out.grad[i][2] = tot_rate + dlog_y[i] / tilted_rate;
    out.grad[i][t_col] = dlog_y[i] / t;
  }
  return out;
}

/**
 * Normal tilt pair, as the dist_id 18 branch in Stan.
 *
 * @param t Point.
 * @param xi Tilt.
 * @param mu Mean.
 * @param sigma Standard deviation.
 * @param grad Whether to compute the partials.
 * @return The pair.
 */
inline TiltPair normal_tilt_pair(double t, double xi, double mu, double sigma,
                                 bool grad) {
  TiltPair out;
  const double z = (t - mu - xi * (sigma * sigma)) / sigma;
  const double xs = xi * sigma;
  const double log_total = xi * mu + 0.5 * (xs * xs);
  const double log_lo = log_std_normal_cdf(z);
  const double log_up = log_std_normal_cdf(-z);
  out.val[0] = log_total + log_lo;
  out.val[1] = log_total + log_up;
  if (grad) {
    const double r_lo = dlog_std_normal_cdf(z, log_lo);
    const double r_up = -dlog_std_normal_cdf(-z, log_up);
    // Partials of z and log_total with respect to xi, mu and sigma
    const double z_p[3]
        = {-sigma, -1 / sigma, -(t - mu) / (sigma * sigma) - xi};
    const double tot_p[3] = {mu + xi * sigma * sigma, xi, xi * xs};
    for (int j = 0; j < 3; ++j) {
      out.grad[0][j] = tot_p[j] + r_lo * z_p[j];
      out.grad[1][j] = tot_p[j] + r_up * z_p[j];
    }
    out.grad[0][t_col] = r_lo / sigma;
    out.grad[1][t_col] = r_up / sigma;
  }
  return out;
}

/**
 * Stacy family tilt pair, as the dist_id 3 and 5 branch in Stan, with
 * primarycensored_stacy_series() for xi other than 0.
 *
 * With x = (t / scale)^shape, the CDF is P(k, x) and each series term is
 *   log a_n = n log(|xi| scale) - log n! + log Gamma(s) - log Gamma(k)
 *             + log P(s, x),  s = k + n / shape.
 * The partials of the sum and of the sum of absolute values are the
 * weighted sums of the term partials.
 *
 * @param t Point.
 * @param xi Tilt.
 * @param shape Power of the Stacy family.
 * @param scale Scale.
 * @param k Shape of the underlying gamma.
 * @param has_k Whether k is a parameter (the generalised gamma).
 * @param grad Whether to compute the partials.
 * @return The pair. Columns 1 to 3 of the partials are shape, scale and k.
 */
inline TiltPair stacy_tilt_pair(double t, double xi, double shape,
                                double scale, double k, bool has_k,
                                bool grad) {
  TiltPair out;
  const double nan = std::numeric_limits<double>::quiet_NaN();
  if (xi == 0) {
    out.val[1] = 0;
    if (t <= 0) {
      return out;
    }
    const double log_ts = std::log(t) - std::log(scale);
    const double log_x = shape * log_ts;
    std::array<double, 2> pq_da;
    const std::array<double, 2> pq
        = gamma_lcdf_lccdf_logx(log_x, k, grad ? &pq_da : nullptr);
    out.val[0] = pq[0];
    out.val[1] = pq[1];
    if (grad) {
      const double log_xf = k * log_x - std::exp(log_x) - stan::math::lgamma(k);
      const double dlog_x[2]
          = {std::exp(log_xf - pq[0]), -std::exp(log_xf - pq[1])};
      for (int i = 0; i < 2; ++i) {
        out.grad[i][1] = dlog_x[i] * log_ts;
        out.grad[i][2] = -dlog_x[i] * shape / scale;
        out.grad[i][3] = has_k ? pq_da[i] : 0;
        out.grad[i][t_col] = dlog_x[i] * shape / t;
      }
    }
    return out;
  }
  out.val[1] = nan;
  if (t <= 0) {
    return out;
  }
  const double log_ts = std::log(t) - std::log(scale);
  const double log_x = shape * log_ts;
  const double x = std::exp(log_x);
  const bool alternating = xi < 0;
  const double log_z = std::log(std::fabs(xi)) + std::log(scale);
  const double log_w = std::log(std::fabs(xi)) + std::log(t);
  const double w = std::exp(log_w);
  auto unreliable = [&]() {
    out.val[0] = nan;
    out.val[2] = nan;
    out.grad = {};
    return out;
  };
  if (w >= 400) {
    return unreliable();
  }
  const double log_k = stan::math::lgamma(k);
  const double digamma_k = grad ? stan::math::digamma(k) : 0;
  double log_pos = 0;
  double log_neg = neg_inf;
  double log_first = 0;
  // Log terms, whether each is negative, and their partials in xi, shape,
  // scale and k
  std::vector<double> log_terms;
  std::vector<bool> negative;
  std::vector<std::array<double, 5>> partials;
  for (int n = 0; n <= 400; ++n) {
    const double s = k + n / shape;
    double dlog_p = 0;
    const double log_p = gamma_lcdf_logx_da(log_x, s, grad ? &dlog_p : nullptr);
    const double log_term = n * log_z - stan::math::lgamma(n + 1)
                            + stan::math::lgamma(s) - log_k + log_p;
    const bool is_neg = alternating && n % 2 == 1;
    if (n == 0) {
      log_first = log_term;
      log_pos = log_term;
    } else if (is_neg) {
      log_neg = stan::math::log_sum_exp(log_neg, log_term);
    } else {
      log_pos = stan::math::log_sum_exp(log_pos, log_term);
    }
    if (grad) {
      // d log P(s, x) / d log x
      const double dlog_x
          = std::exp(s * log_x - x - stan::math::lgamma(s) - log_p);
      const double ds = stan::math::digamma(s) + dlog_p;
      log_terms.push_back(log_term);
      negative.push_back(is_neg);
      partials.push_back({n / xi,
                          -ds * n / (shape * shape) + dlog_x * log_ts,
                          n / scale - dlog_x * shape / scale,
                          ds - digamma_k, dlog_x * shape / t});
    }
    const double log_abs = stan::math::log_sum_exp(log_pos, log_neg);
    if (alternating && log_abs - log_first > std::log(1e5)) {
      return unreliable();
    }
    if (n + 1 > w
        && log_term + log_w - std::log(n + 1 - w)
               < log_abs + std::log(1e-17)) {
      if (log_pos <= log_neg) {
        return unreliable();
      }
      const double log_sum = log_neg == neg_inf
                                 ? log_pos
                                 : stan::math::log_diff_exp(log_pos, log_neg);
      if (log_abs - log_sum > std::log(1e5)) {
        return unreliable();
      }
      out.val[0] = log_sum;
      out.val[2] = log_abs - log_sum;
      if (grad) {
        // d log_sum and d log_abs
        std::array<double, 5> g_sum{}, g_abs{};
        for (std::size_t m = 0; m < log_terms.size(); ++m) {
          const double w_sum = std::exp(log_terms[m] - log_sum);
          const double w_abs = std::exp(log_terms[m] - log_abs);
          for (int j = 0; j < 5; ++j) {
            g_sum[j] += (negative[m] ? -w_sum : w_sum) * partials[m][j];
            g_abs[j] += w_abs * partials[m][j];
          }
        }
        for (int j = 0; j < 5; ++j) {
          out.grad[0][j] = g_sum[j];
          out.grad[2][j] = g_abs[j] - g_sum[j];
        }
        if (!has_k) {
          out.grad[0][3] = 0;
          out.grad[2][3] = 0;
        }
      }
      return out;
    }
  }
  return unreliable();
}

/**
 * Lognormal tilt pair for xi >= 0, as primarycensored_lognormal_tilt_pair()
 * in Stan.
 *
 * With z = (log t - mu) / sigma, the CDF is Phi(z). For xi > 0 the lower
 * transform is the series of the terms
 *   log a_k = k (log xi + mu) + (k sigma)^2 / 2 - log k!
 *             + log Phi(z - k sigma),
 * whose partials are summed with weights a_k / sum as the terms are added.
 * The upper transform is inf for xi > 0, as the total diverges.
 *
 * @param t Point.
 * @param xi Tilt, at least 0.
 * @param mu Mean of the log of the delay.
 * @param sigma Standard deviation of the log of the delay.
 * @param grad Whether to compute the partials.
 * @return The pair.
 */
inline TiltPair lognormal_tilt_pair(double t, double xi, double mu,
                                    double sigma, bool grad) {
  TiltPair out;
  const double inf = std::numeric_limits<double>::infinity();
  if (xi < 0) {
    throw std::domain_error(
        "primarycensored_tilt_pair: the lognormal needs xi >= 0");
  }
  out.val[1] = xi == 0 ? 0 : inf;
  if (t <= 0) {
    return out;
  }
  const double z = (std::log(t) - mu) / sigma;
  const double log_lo = log_std_normal_cdf(z);
  if (xi == 0) {
    const double log_up = log_std_normal_cdf(-z);
    out.val[0] = log_lo;
    out.val[1] = log_up;
    if (grad) {
      const double r_lo = dlog_std_normal_cdf(z, log_lo);
      const double r_up = -dlog_std_normal_cdf(-z, log_up);
      out.grad[0][1] = -r_lo / sigma;
      out.grad[0][2] = -r_lo * z / sigma;
      out.grad[1][1] = -r_up / sigma;
      out.grad[1][2] = -r_up * z / sigma;
      out.grad[0][t_col] = r_lo / (sigma * t);
      out.grad[1][t_col] = r_up / (sigma * t);
    }
    return out;
  }
  const double log_xi_mu = std::log(xi) + mu;
  double log_factorial = 0;
  double total = log_lo;
  // Partials of the log sum in xi, mu, sigma and t
  double g[4] = {0, 0, 0, 0};
  if (grad) {
    const double r = dlog_std_normal_cdf(z, log_lo);
    g[1] = -r / sigma;
    g[2] = -r * z / sigma;
    g[3] = r / (sigma * t);
  }
  for (int k = 1; k <= 20000; ++k) {
    log_factorial += std::log(k);
    const double ks = k * sigma;
    const double z_k = z - ks;
    const double log_phi_k = log_std_normal_cdf(z_k);
    const double log_term
        = k * log_xi_mu + 0.5 * (ks * ks) - log_factorial + log_phi_k;
    const double new_total = stan::math::log_sum_exp(total, log_term);
    if (grad) {
      const double keep = std::exp(total - new_total);
      const double weight = std::exp(log_term - new_total);
      const double r = dlog_std_normal_cdf(z_k, log_phi_k);
      g[0] = g[0] * keep + weight * k / xi;
      g[1] = g[1] * keep + weight * (k - r / sigma);
      g[2] = g[2] * keep + weight * (k * ks - r * (z / sigma + k));
      g[3] = g[3] * keep + weight * r / (sigma * t);
    }
    total = new_total;
    if (k > xi * t && log_term < total - 41) {
      out.val[0] = total;
      for (int j = 0; j < 3; ++j) {
        out.grad[0][j] = g[j];
      }
      out.grad[0][t_col] = g[3];
      return out;
    }
  }
  throw std::domain_error(
      "The lognormal tilt transform needs more than 20000 terms. Use the "
      "numerical path.");
}

/**
 * Tilt pair for the given delay distribution.
 *
 * @param t Point.
 * @param dist_id 1 (lognormal, for xi >= 0), 2 (gamma), 3 (Weibull),
 *   4 (exponential), 5 (generalised gamma) or 18 (normal).
 * @param xi Tilt.
 * @param p Delay parameters.
 * @param grad Whether to compute the partials.
 * @return The pair.
 */
inline TiltPair tilt_pair(double t, int dist_id, double xi,
                          const std::vector<double>& p, bool grad) {
  if (dist_id == 1) {
    return lognormal_tilt_pair(t, xi, p[0], p[1], grad);
  } else if (dist_id == 2) {
    return gamma_tilt_pair(t, xi, p[0], p[1], grad);
  } else if (dist_id == 3) {
    return stacy_tilt_pair(t, xi, p[0], p[1], 1.0, false, grad);
  } else if (dist_id == 4) {
    return exponential_tilt_pair(t, xi, p[0], grad);
  } else if (dist_id == 5) {
    return stacy_tilt_pair(t, xi, p[0], p[1], p[2], true, grad);
  }
  return normal_tilt_pair(t, xi, p[0], p[1], grad);
}

}  // namespace internal

/**
 * Log of the tilt transforms over the lower and upper support.
 *
 * The point t may be a var, as for the split point of the truncated
 * logistic series. The last argument is the output stream of the Stan
 * calling convention, which is not used.
 *
 * @param t Point (var or double).
 * @param dist_id 1 (lognormal, for xi >= 0), 2 (gamma), 3 (Weibull),
 *   4 (exponential), 5 (generalised gamma) or 18 (normal).
 * @param xi Tilt (var or double).
 * @param params Delay parameters (var or double).
 * @return Vector [log T_f(xi; t), log(T_f(xi; Inf) - T_f(xi; t)), log loss],
 *   as log_tilt_transform_pair_shared() in Stan.
 */
template <typename T_t, typename T_xi, typename T_p,
          stan::require_all_stan_scalar_t<T_t, T_xi, T_p>* = nullptr>
inline Eigen::Matrix<stan::return_type_t<T_t, T_xi, T_p>, -1, 1>
primarycensored_tilt_pair(const T_t& t, const int& dist_id, const T_xi& xi,
                          const std::vector<T_p>& params,
                          std::ostream* /* pstream__ */) {
  constexpr const char* fn = "primarycensored_tilt_pair";
  if (dist_id < 1 || (dist_id > 5 && dist_id != 18)) {
    throw std::domain_error(std::string(fn)
                            + ": dist_id must be 1, 2, 3, 4, 5 or 18");
  }
  const std::size_t n_params = dist_id == 4 ? 1 : (dist_id == 5 ? 3 : 2);
  if (params.size() < n_params) {
    throw std::domain_error(std::string(fn) + ": too few params");
  }
  constexpr bool t_var = stan::is_var<T_t>::value;
  constexpr bool xi_var = stan::is_var<T_xi>::value;
  constexpr bool p_var = stan::is_var<T_p>::value;
  std::vector<double> p(n_params);
  for (std::size_t j = 0; j < n_params; ++j) {
    p[j] = stan::math::value_of(params[j]);
  }
  const internal::TiltPair pair = internal::tilt_pair(
      stan::math::value_of(t), dist_id, stan::math::value_of(xi), p,
      t_var || xi_var || p_var);
  Eigen::VectorXd val(3);
  for (int i = 0; i < 3; ++i) {
    val(i) = pair.val[i];
  }
  using ret_t = stan::return_type_t<T_t, T_xi, T_p>;
  if constexpr (!t_var && !xi_var && !p_var) {
    return val.template cast<ret_t>();
  } else {
    using stan::math::var;
    // Columns of the partials that are inputs: xi, the parameters used and t
    std::vector<int> cols;
    stan::arena_t<std::vector<var>> ins;
    if constexpr (xi_var) {
      cols.push_back(0);
      ins.push_back(xi);
    }
    if constexpr (p_var) {
      for (std::size_t j = 0; j < n_params; ++j) {
        cols.push_back(static_cast<int>(j) + 1);
        ins.push_back(params[j]);
      }
    }
    if constexpr (t_var) {
      cols.push_back(internal::t_col);
      ins.push_back(t);
    }
    const int n_in = static_cast<int>(cols.size());
    stan::arena_t<Eigen::MatrixXd> jac(3, n_in);
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < n_in; ++j) {
        jac(i, j) = pair.grad[i][cols[j]];
      }
    }
    stan::arena_t<Eigen::Matrix<var, -1, 1>> res = val.template cast<var>();
    stan::math::reverse_pass_callback([res, jac, ins, n_in]() mutable {
      for (int j = 0; j < n_in; ++j) {
        double adj = 0;
        for (int i = 0; i < 3; ++i) {
          // Elements that do not depend on the inputs take no adjoint, as
          // in Stan
          if (jac(i, j) != 0) {
            adj += res(i).adj() * jac(i, j);
          }
        }
        ins[j].adj() += adj;
      }
    });
    return res;
  }
}

}  // namespace primarycensored

namespace stan {
namespace math {
using ::primarycensored::primarycensored_tilt_pair;
}  // namespace math
}  // namespace stan

#endif
