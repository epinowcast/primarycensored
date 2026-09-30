/*
 * Truncated Gumbel primary event distribution
 *
 * A Gumbel distribution with location mu and scale beta, truncated to
 * [xmin, xmax]. With G(z) = exp(-s(z)) and s(z) = exp(-(z - mu) / beta) the
 * density is G'(x) / (G(xmax) - G(xmin)). Differences of s are written as
 * s(x) (exp((x - xmin) / beta) - 1) so that the terms do not cancel when mu
 * is far from the window. The R equivalents are dtgumbel(), ptgumbel() and
 * rtgumbel(). It is primary distribution 4 with primary_params = [mu, beta].
 */

/**
  * Log of the difference of s across the window of a truncated Gumbel
  * @ingroup truncated_gumbel_distributions
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
  * Log of 1 - exp(-exp(x))
  * @ingroup truncated_gumbel_distributions
  *
  * Below -37 the value is x, where exp(x) would underflow for a window far
  * above the location. Above 700 it is 0 to double precision.
  *
  * @param x Real number
  * @return log(1 - exp(-exp(x)))
  */
real tgumbel_log1m_exp_neg_exp(real x) {
  if (x < -37) {
    return x;
  }
  if (x > 700) {
    return 0;
  }
  return log1m_exp(-exp(x));
}

/**
  * Log of the normalisation of a truncated Gumbel
  * @ingroup truncated_gumbel_distributions
  *
  * The log of (G(xmax) - G(xmin)) / G(xmax), as G(xmax) cancels against the
  * density.
  *
  * @param xmin Lower bound of the distribution
  * @param xmax Upper bound of the distribution
  * @param mu Location
  * @param beta Scale, positive
  * @return log(1 - exp(-(s(xmin) - s(xmax))))
  */
real tgumbel_log_norm(real xmin, real xmax, real mu, real beta) {
  return tgumbel_log1m_exp_neg_exp(
    tgumbel_log_delta_window(xmin, xmax, mu, beta)
  );
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
  real exponent = 0;
  if (x < xmax) {
    exponent = exp(-(xmax - mu) / beta + log_diff_exp((xmax - x) / beta, 0));
  }
  return -log(beta) + log_s - exponent - tgumbel_log_norm(xmin, xmax, mu, beta);
}

/**
  * Truncated Gumbel log cumulative distribution function (log CDF)
  * @ingroup truncated_gumbel_distributions
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
  real log_norm = tgumbel_log1m_exp_neg_exp(log_delta_window);
  return tgumbel_log1m_exp_neg_exp(log_delta_lower) - exp(log_delta_upper)
         - log_norm;
}

/**
  * Truncated Gumbel random number generator
  * @ingroup truncated_gumbel_distributions
  *
  * Inverts the upper tail, which does not cancel.
  *
  * @param xmin Lower bound of the distribution
  * @param xmax Upper bound of the distribution
  * @param mu Location of the Gumbel before truncation
  * @param beta Scale of the Gumbel before truncation, positive
  * @return A random draw from the truncated Gumbel distribution
  */
real tgumbel_rng(real xmin, real xmax, real mu, real beta) {
  real u = uniform_rng(0, 1);
  real log_norm = tgumbel_log_norm(xmin, xmax, mu, beta);
  real log_delta_upper = log(-log1m_exp(log1m(u) + log_norm));
  // x = xmax - beta log(1 + delta_upper / s(xmax)), which does not cancel
  real log_ratio = log_delta_upper + (xmax - mu) / beta;
  return fmin(fmax(xmax - beta * log1p_exp(log_ratio), xmin), xmax);
}
