/**
  * Log of the mass of the logistic distribution on an interval
  * @ingroup truncated_logistic_distributions
  *
  * Evaluates log(L(hi) - L(lo)) for lo <= hi from the lower tails of the
  * logistic CDF L where lo is below the location and from the upper tails
  * otherwise. The difference is then not lost to rounding when the interval is
  * far from the location, where L(hi) and L(lo) both round to 0 or 1. It is
  * `-inf` where the mass is zero to rounding.
  *
  * @param lo Lower end of the interval
  * @param hi Upper end of the interval
  * @param location Location of the logistic distribution
  * @param scale Scale of the logistic distribution
  *
  * @return log(L(hi) - L(lo))
  */
real tlogis_log_mass(real lo, real hi, real location, real scale) {
  if (lo >= location) {
    return primarycensored_log_diff_exp(
      log_inv_logit(-(lo - location) / scale),
      log_inv_logit(-(hi - location) / scale)
    );
  }
  return primarycensored_log_diff_exp(
    log_inv_logit((hi - location) / scale),
    log_inv_logit((lo - location) / scale)
  );
}

/**
  * Truncated logistic log probability density function (log PDF)
  * @ingroup truncated_logistic_distributions
  *
  * The density of a logistic distribution truncated to [xmin, xmax],
  * L'(x) / (L(xmax) - L(xmin)) with L the logistic CDF. This is the primary
  * event distribution with primary_id 3.
  *
  * @param x Value at which to evaluate the log PDF
  * @param xmin Lower bound of the distribution
  * @param xmax Upper bound of the distribution
  * @param location Location of the logistic distribution before truncation
  * @param scale Scale of the logistic distribution before truncation, positive
  *
  * @return The log PDF evaluated at x
  */
real tlogis_lpdf(real x, real xmin, real xmax, real location, real scale) {
  if (x < xmin || x > xmax) {
    return negative_infinity();
  }
  return logistic_lpdf(x | location, scale)
         - tlogis_log_mass(xmin, xmax, location, scale);
}

/**
  * Truncated logistic log cumulative distribution function (log CDF)
  * @ingroup truncated_logistic_distributions
  *
  * @param x Value at which to evaluate the log CDF
  * @param xmin Lower bound of the distribution
  * @param xmax Upper bound of the distribution
  * @param location Location of the logistic distribution before truncation
  * @param scale Scale of the logistic distribution before truncation, positive
  *
  * @return The log CDF evaluated at x
  */
real tlogis_lcdf(real x, real xmin, real xmax, real location, real scale) {
  if (x <= xmin) {
    return negative_infinity();
  }
  if (x >= xmax) {
    return 0;
  }
  return tlogis_log_mass(xmin, x, location, scale)
         - tlogis_log_mass(xmin, xmax, location, scale);
}

/**
  * Truncated logistic cumulative distribution function (CDF)
  * @ingroup truncated_logistic_distributions
  *
  * @param x Value at which to evaluate the CDF
  * @param xmin Lower bound of the distribution
  * @param xmax Upper bound of the distribution
  * @param location Location of the logistic distribution before truncation
  * @param scale Scale of the logistic distribution before truncation, positive
  *
  * @return The CDF evaluated at x
  */
real tlogis_cdf(real x, real xmin, real xmax, real location, real scale) {
  return exp(tlogis_lcdf(x | xmin, xmax, location, scale));
}

/**
  * Truncated logistic random number generator
  * @ingroup truncated_logistic_distributions
  *
  * Inverse transform sampling on the log scale from the tails of the logistic
  * CDF that are nearer the interval, see tlogis_log_mass().
  *
  * @param xmin Lower bound of the distribution
  * @param xmax Upper bound of the distribution
  * @param location Location of the logistic distribution before truncation
  * @param scale Scale of the logistic distribution before truncation, positive
  *
  * @return A random draw from the truncated logistic distribution
  */
real tlogis_rng(real xmin, real xmax, real location, real scale) {
  real u = uniform_rng(0, 1);
  real z;
  if (xmin >= location) {
    // The upper tail S at the draw is (1 - u) S(xmin) + u S(xmax)
    real log_tail = log_mix(
      u, log_inv_logit(-(xmax - location) / scale),
      log_inv_logit(-(xmin - location) / scale)
    );
    z = log1m_exp(log_tail) - log_tail;
  } else {
    real log_tail = log_mix(
      u, log_inv_logit((xmax - location) / scale),
      log_inv_logit((xmin - location) / scale)
    );
    z = log_tail - log1m_exp(log_tail);
  }
  return fmin(fmax(location + scale * z, xmin), xmax);
}

/**
  * Times at which to restart the ODE for a narrow truncated logistic primary
  * @ingroup truncated_logistic_distributions
  *
  * The integrand of the primary event censored CDF is the delay CDF times the
  * primary density at d - t. A narrow logistic density has its mass within a
  * few scales of its location, or of the edge of the window nearest to it
  * when the location is outside. A solver that takes large steps across the
  * flat region steps over that spike, and returns a CDF that is too small
  * without an error. The integral is split at multiples of the scale from
  * that centre, so each part is solved on its own. The density falls by
  * about e^-x at x scales, so the mass beyond 30 scales is below 1e-13.
  *
  * The integral does not depend on where it is split, so the gradients with
  * respect to the times cancel.
  *
  * @param d Delay, the end of the integral
  * @param start Start of the integral
  * @param pwindow Primary event window
  * @param location Location of the logistic distribution before truncation
  * @param scale Scale of the logistic distribution before truncation, positive
  *
  * @return Strictly increasing times in (start, d), possibly empty
  */
vector tlogis_spike_times(data real d, real start, data real pwindow,
                          real location, real scale) {
  // Multiples of the scale in the primary event time, so decreasing to
  // give increasing times
  array[11] real multiples = {30, 12, 6, 3, 1.5, 0, -1.5, -3, -6, -12, -30};
  real centre = fmin(fmax(location, 0), pwindow);
  vector[11] times;
  int n = 0;
  for (k in 1:11) {
    real time = d - (centre + scale * multiples[k]);
    if (time > start && time < d && (n == 0 || time > times[n])) {
      n += 1;
      times[n] = time;
    }
  }
  return head(times, n);
}
