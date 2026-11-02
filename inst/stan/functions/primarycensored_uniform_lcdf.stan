/**
  * Compute the uniform primary censored log CDF at integer delays
  * @ingroup primary_event_analytical_distributions
  *
  * For lognormal (1), gamma (2) and Weibull (3) delays. This Stan body is
  * the reference implementation. `pcd_load_stan_functions(cpp = TRUE)` in R
  * replaces it with a declaration, and the C++ version in the header given
  * by `pcd_stan_header()` is used instead. That version returns the same
  * values and computes the gradient in closed form as one autodiff node.
  *
  * @param start First delay to compute
  * @param n Last delay to compute, and the length of the result
  * @param dist_id Distribution identifier (1, 2 or 3)
  * @param params Array of distribution parameters
  * @param pwindow Primary event window, a positive integer
  *
  * @return Vector whose element d is the log CDF at d, for d in start:n.
  * Elements before start are not computed.
  */
vector primarycensored_uniform_lcdf_vectorized(data int start, data int n,
                                               data int dist_id,
                                               array[] real params,
                                               data real pwindow) {
  return primarycensored_uniform_terms_lcdf_vectorized(
    start, n, dist_id, params, pwindow
  );
}
