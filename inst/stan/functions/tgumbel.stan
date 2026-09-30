/*
 * Truncated Gumbel primary event distribution
 *
 * A Gumbel distribution with location mu and scale beta, truncated to
 * [xmin, xmax]. With G(z) = exp(-s(z)) and s(z) = exp(-(z - mu) / beta) the
 * density is G'(x) / (G(xmax) - G(xmin)) and the CDF is
 * (G(x) - G(xmin)) / (G(xmax) - G(xmin)). Differences of s are written as
 * s(x) (exp((x - xmin) / beta) - 1) so that the terms do not cancel when mu
 * is far from the window. The R equivalents are dtgumbel(), ptgumbel() and
 * rtgumbel(). It is primary distribution 4 with primary_params = [mu, beta].
 */

/**
  * Log of the difference of s across the window of a truncated Gumbel
  * @ingroup truncated_gumbel_distributions
  *
  * The log of s(xmin) - s(xmax), which is s(xmax) (exp((xmax - xmin) / beta)
  * - 1), with s(z) = exp(-(z - mu) / beta).
  *
  * @param xmin Lower bound of the distribution
  * @param xmax Upper bound of the distribution
  * @param mu Location
  * @param beta Scale, positive
  * @return log(s(xmin) - s(xmax))
  */
real tgumbel_log_delta_window(real xmin, real xmax, real mu, real beta) {
  return -(xmax - mu) / beta + log_diff_exp((xmax - xmin) / beta, 0);
}

/**
  * Log of the normalisation of a truncated Gumbel
  * @ingroup truncated_gumbel_distributions
  *
  * The log of G(xmax) - G(xmin), which is -s(xmax) plus the log of
  * 1 - exp(-(s(xmin) - s(xmax))).
  *
  * @param xmin Lower bound of the distribution
  * @param xmax Upper bound of the distribution
  * @param mu Location
  * @param beta Scale, positive
  * @return log(G(xmax) - G(xmin))
  */
real tgumbel_log_norm(real xmin, real xmax, real mu, real beta) {
  real log_delta = tgumbel_log_delta_window(xmin, xmax, mu, beta);
  // Beyond this the log term is 0 to double precision, and exp() overflows
  if (log_delta > 700) {
    return -exp(-(xmax - mu) / beta);
  }
  return -exp(-(xmax - mu) / beta) + log1m_exp(-exp(log_delta));
}

/**
  * Truncated Gumbel log probability density function (log PDF)
  * @ingroup truncated_gumbel_distributions
  *
  * @param x Value at which to evaluate the log PDF
  * @param xmin Lower bound of the distribution
  * @param xmax Upper bound of the distribution
  * @param mu Location of the Gumbel before truncation
  * @param beta Scale of the Gumbel before truncation, positive
  * @return The log PDF evaluated at x, `-inf` outside [xmin, xmax]
  */
real tgumbel_lpdf(real x, real xmin, real xmax, real mu, real beta) {
  if (x < xmin || x > xmax) {
    return negative_infinity();
  }
  real log_s = -(x - mu) / beta;
  return -log(beta) + log_s - exp(log_s)
         - tgumbel_log_norm(xmin, xmax, mu, beta);
}

/**
  * Truncated Gumbel log cumulative distribution function (log CDF)
  * @ingroup truncated_gumbel_distributions
  *
  * With delta_lower = s(xmin) - s(x), delta_upper = s(x) - s(xmax) and
  * delta_window = s(xmin) - s(xmax), the CDF is
  * (exp(delta_lower) - 1) / (exp(delta_window) - 1), which on the log scale
  * is log(1 - exp(-delta_lower)) - delta_upper - log(1 - exp(-delta_window)).
  *
  * @param x Value at which to evaluate the log CDF
  * @param xmin Lower bound of the distribution
  * @param xmax Upper bound of the distribution
  * @param mu Location of the Gumbel before truncation
  * @param beta Scale of the Gumbel before truncation, positive
  * @return The log CDF evaluated at x
  */
real tgumbel_lcdf(real x, real xmin, real xmax, real mu, real beta) {
  if (x <= xmin) {
    return negative_infinity();
  }
  if (x >= xmax) {
    return 0;
  }
  real log_delta_window = tgumbel_log_delta_window(xmin, xmax, mu, beta);
  real log_delta_lower = -(x - mu) / beta
                         + log_diff_exp((x - xmin) / beta, 0);
  real log_delta_upper = -(xmax - mu) / beta
                         + log_diff_exp((xmax - x) / beta, 0);
  real log_norm = log_delta_window > 700
                  ? 0 : log1m_exp(-exp(log_delta_window));
  return log1m_exp(-exp(log_delta_lower)) - exp(log_delta_upper) - log_norm;
}

/**
  * Truncated Gumbel random number generator
  * @ingroup truncated_gumbel_distributions
  *
  * Inverts the upper tail, 1 - F(x) = (1 - exp(-delta_upper)) /
  * (1 - exp(-delta_window)), which adds delta_upper = s(x) - s(xmax) to
  * s(xmax) and so does not cancel.
  *
  * @param xmin Lower bound of the distribution
  * @param xmax Upper bound of the distribution
  * @param mu Location of the Gumbel before truncation
  * @param beta Scale of the Gumbel before truncation, positive
  * @return A random draw from the truncated Gumbel distribution
  */
real tgumbel_rng(real xmin, real xmax, real mu, real beta) {
  real u = uniform_rng(0, 1);
  real log_delta_window = tgumbel_log_delta_window(xmin, xmax, mu, beta);
  real log_norm = log_delta_window > 700
                  ? 0 : log1m_exp(-exp(log_delta_window));
  real delta_upper = -log1m_exp(log1m(u) + log_norm);
  real s = exp(-(xmax - mu) / beta) + delta_upper;
  return fmin(fmax(mu - beta * log(s), xmin), xmax);
}
