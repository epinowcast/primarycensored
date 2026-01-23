/**
  * Exponential growth probability density function (PDF)
  * @ingroup exponential_growth_distributions
  *
  * @param x Value at which to evaluate the PDF
  * @param xmin Lower bound of the distribution
  * @param xmax Upper bound of the distribution
  * @param r Rate parameter for exponential growth
  * @return The PDF evaluated at x
  */
real expgrowth_pdf(real x, real xmin, real xmax, real r) {
  if (x < xmin || x > xmax) {
    return 0;
  }
  if (abs(r) < 1e-10) {
    return 1 / (xmax - xmin);
  }
  return r * exp(r * (x - xmin)) / (exp(r * xmax) - exp(r * xmin));
}

/**
  * Exponential growth log probability density function (log PDF)
  * @ingroup exponential_growth_distributions
  *
  * @param x Value at which to evaluate the log PDF
  * @param xmin Lower bound of the distribution
  * @param xmax Upper bound of the distribution
  * @param r Rate parameter for exponential growth
  * @return The log PDF evaluated at x
  */
real expgrowth_lpdf(real x, real xmin, real xmax, real r) {
  if (x < xmin || x > xmax) {
    return negative_infinity();
  }
  if (abs(r) < 1e-10) {
    return -log(xmax - xmin);
  }
  return log(r) + r * (x - xmin) - log(exp(r * xmax) - exp(r * xmin));
}

/**
  * Exponential growth cumulative distribution function (CDF)
  * @ingroup exponential_growth_distributions
  *
  * @param x Value at which to evaluate the CDF
  * @param xmin Lower bound of the distribution
  * @param xmax Upper bound of the distribution
  * @param r Rate parameter for exponential growth
  * @return The CDF evaluated at x
  */
real expgrowth_cdf(real x, real xmin, real xmax, real r) {
  if (x < xmin) {
    return 0;
  }
  if (x > xmax) {
    return 1;
  }
  if (abs(r) < 1e-10) {
    return (x - xmin) / (xmax - xmin);
  }
  return (exp(r * (x - xmin)) - exp(r * xmin)) / (exp(r * xmax) - exp(r * xmin));
}

/**
  * Exponential growth log cumulative distribution function (log CDF)
  * @ingroup exponential_growth_distributions
  *
  * @param x Value at which to evaluate the log CDF
  * @param xmin Lower bound of the distribution
  * @param xmax Upper bound of the distribution
  * @param r Rate parameter for exponential growth
  * @return The log CDF evaluated at x
  */
real expgrowth_lcdf(real x, real xmin, real xmax, real r) {
  if (x < xmin) {
    return negative_infinity();
  }
  if (x > xmax) {
    return 0;
  }
  return log(expgrowth_cdf(x | xmin, xmax, r));
}

/**
  * Exponential growth random number generator
  * @ingroup exponential_growth_distributions
  *
  * @param xmin Lower bound of the distribution
  * @param xmax Upper bound of the distribution
  * @param r Rate parameter for exponential growth
  * @return A random draw from the exponential growth distribution
  */
real expgrowth_rng(real xmin, real xmax, real r) {
  real u = uniform_rng(0, 1);
  if (abs(r) < 1e-10) {
    return xmin + u * (xmax - xmin);
  }
  return xmin + log(u * (exp(r * xmax) - exp(r * xmin)) + exp(r * xmin)) / r;
}
