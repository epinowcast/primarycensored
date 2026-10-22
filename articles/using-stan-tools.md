# How to use primarycensored with Stan

## 1 Introduction

### 1.1 What are we going to do in this vignette

In this vignette, we’ll explore how to use the `primarycensored` package
in your Stan modelling workflow. We’ll cover the following key points:

1.  Introduction to Stan and its relevance for our analysis
2.  Overview of the packages we’ll be using
3.  How to access and use Stan functions provided by `primarycensored`
4.  Methods for integrating these Stan functions into your workflow

If you are instead interested in fitting a delay distribution using
`primarycensored` in R see the
[`vignette("fitting-dists-with-stan")`](https://primarycensored.epinowcast.org/articles/fitting-dists-with-stan.md)
vignette or [`epidist`](https://epidist.epinowcast.org) package (which
uses `primarycensored` under the hood).

### 1.2 What is Stan and why are we using it

Stan is a probabilistic programming language for statistical inference.
It provides a flexible and efficient platform for Bayesian modeling and
is widely used in various fields of data science and statistics. In this
vignette, we’ll use Stan in conjunction with `primarycensored` to
perform Bayesian inference on censored data.

For more information on Stan:

- [Stan’s official website](https://mc-stan.org/)
- [Stan documentation](https://mc-stan.org/users/documentation/)
- [Stan forums](https://discourse.mc-stan.org/) for community support
  and discussions

### 1.3 Packages used in this vignette

Alongside the `primarycensored` package we will use the `cmdstanr`
package for interfacing with cmdstan.

``` r

library(primarycensored)
library(cmdstanr)
```

## 2 Using Stan code in primarycensored

`primarycensored` includes a set of Stan functions that mirror the R
functions in `primarycensored`. Documentation for these functions can be
found [here](https://primarycensored.epinowcast.org/stan/). We support a
range of approaches to integrate this Stan code into your workflow.

### 2.1 Checking available Stan functions using `pcd_stan_functions()`

Aside from reading the documentation it is also possible to list the
available Stan functions using a helper function directly in R.

``` r

pcd_stan_functions()
```

    ##  [1] "expgrowth_pdf"                                       
    ##  [2] "expgrowth_lpdf"                                      
    ##  [3] "expgrowth_cdf"                                       
    ##  [4] "expgrowth_lcdf"                                      
    ##  [5] "expgrowth_rng"                                       
    ##  [6] "pstep_lcdf"                                          
    ##  [7] "hazards_to_pmf"                                      
    ##  [8] "phazard_lcdf"                                        
    ##  [9] "primary_lcdf_vec"                                    
    ## [10] "discretestep_lcdf"                                   
    ## [11] "discretehazard_lcdf"                                 
    ## [12] "check_for_uniform_terms"                             
    ## [13] "check_for_analytical"                                
    ## [14] "check_for_analytical_params"                         
    ## [15] "primarycensored_uniform_lcdf_from_terms"             
    ## [16] "primarycensored_gamma_uniform_terms"                 
    ## [17] "primarycensored_lognormal_uniform_terms"             
    ## [18] "log_weibull_g"                                       
    ## [19] "primarycensored_weibull_uniform_terms"               
    ## [20] "primarycensored_gengamma_uniform_terms"              
    ## [21] "primarycensored_uniform_terms"                       
    ## [22] "primarycensored_gamma_uniform_lcdf"                  
    ## [23] "primarycensored_lognormal_uniform_lcdf"              
    ## [24] "primarycensored_weibull_uniform_lcdf"                
    ## [25] "primarycensored_gengamma_uniform_lcdf"               
    ## [26] "primarycensored_analytical_lcdf_raw"                 
    ## [27] "primarycensored_analytical_lcdf"                     
    ## [28] "primarycensored_analytical_cdf"                      
    ## [29] "check_for_analytical_vectorized"                     
    ## [30] "primarycensored_analytical_lcdf_vectorized"          
    ## [31] "primarycensored_log_diff_exp"                        
    ## [32] "primarycensored_log_std_normal_cdf"                  
    ## [33] "check_for_tilt_transform"                            
    ## [34] "log_tilt_transform_pair"                             
    ## [35] "primarycensored_moments_from_partial"                
    ## [36] "primarycensored_gamma_tilt_moments"                  
    ## [37] "primarycensored_tilt_moments"                        
    ## [38] "check_for_exptilt"                                   
    ## [39] "primarycensored_exptilt_form"                        
    ## [40] "primarycensored_tail_diff"                           
    ## [41] "primarycensored_exptilt_lcdf_from_terms"             
    ## [42] "primarycensored_exptilt_small_window_lcdf_from_terms"
    ## [43] "primarycensored_exptilt_small_delay_lcdf_from_terms" 
    ## [44] "primarycensored_exptilt_lcdf"                        
    ## [45] "gamma_log_lead_logx"                                 
    ## [46] "gamma_log_pq_logx"                                   
    ## [47] "gamma_lcdf_logx"                                     
    ## [48] "gamma_lcdf_logx_pair"                                
    ## [49] "gamma_log_rate"                                      
    ## [50] "gengamma_lcdf"                                       
    ## [51] "dist_has_positive_support"                           
    ## [52] "lognormal_lcdf_underflows"                           
    ## [53] "dist_lcdf"                                           
    ## [54] "primary_lcdf"                                        
    ## [55] "primary_lpdf"                                        
    ## [56] "primarycensored_ode"                                 
    ## [57] "primarycensored_log_normalizer"                      
    ## [58] "primarycensored_apply_truncation"                    
    ## [59] "primarycensored_truncation_bounds"                   
    ## [60] "primarycensored_numeric_cdf"                         
    ## [61] "primarycensored_cdf"                                 
    ## [62] "primarycensored_lcdf"                                
    ## [63] "primarycensored_lpmf"                                
    ## [64] "primarycensored_pmf"                                 
    ## [65] "primarycensored_lcdf_vectorized"                     
    ## [66] "primarycensored_sone_lpmf_vectorized"                
    ## [67] "primarycensored_sone_pmf_vectorized"

### 2.2 Accessing Stan functions

Stan functions are accessed using the
[`pcd_load_stan_functions()`](https://primarycensored.epinowcast.org/reference/pcd_load_stan_functions.md)
function. This function takes the name of the function as an argument
and returns the function as a string. It can additionally write the
functions to a file and wrap them in a `functions{}` block.

``` r

pcd_load_stan_functions("primarycensored_lpmf")
```

    ## [1] "// Stan functions from primarycensored version 1.6.0\nreal primarycensored_lpmf(data int d, data int dist_id, array[] real params,\n                                data real pwindow, data real d_upper,\n                                data real L, data real D, data int primary_id,\n                                array[] real primary_params) {\n  if (d_upper > D) {\n    reject(\"Upper truncation point is greater than D. It is \", d_upper,\n           \" and D is \", D, \". Resolve this by increasing D to be greater or equal to d + swindow or decreasing swindow.\");\n  }\n  if (d_upper <= d) {\n    reject(\"Upper truncation point is less than or equal to d. It is \", d_upper,\n           \" and d is \", d, \". Resolve this by increasing d to be less than d_upper.\");\n  }\n  if (d < L) {\n    return negative_infinity();\n  }\n  real log_cdf_upper = primarycensored_lcdf(\n    d_upper | dist_id, params, pwindow,\n    dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),\n    positive_infinity(), primary_id, primary_params\n  );\n  real log_cdf_lower = primarycensored_lcdf(\n    d | dist_id, params, pwindow,\n    dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),\n    positive_infinity(), primary_id, primary_params\n  );\n\n  // Apply truncation normalization: log((F(d_upper) - F(d)) / (F(D) - F(L))).\n  // Skip when F(L) = 0 makes it a no-op (positive support, L <= 0).\n  if (!is_inf(D) || L > 0 ||\n      (!is_inf(L) && !dist_has_positive_support(dist_id))) {\n    real log_cdf_D;\n    real log_cdf_L;\n\n    // Get CDF at lower truncation point L\n    if (is_inf(L)) {\n      // No left truncation (L = -inf sentinel)\n      log_cdf_L = negative_infinity();\n    } else if (d == L) {\n      // Reuse already computed CDF at d\n      log_cdf_L = log_cdf_lower;\n    } else {\n      // Compute CDF at L directly\n      log_cdf_L = primarycensored_lcdf(\n        L | dist_id, params, pwindow,\n        dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),\n        positive_infinity(), primary_id, primary_params\n      );\n    }\n\n    // Get CDF at upper truncation point D\n    if (d_upper == D) {\n      log_cdf_D = log_cdf_upper;\n    } else if (is_inf(D)) {\n      log_cdf_D = 0;\n    } else {\n      log_cdf_D = primarycensored_lcdf(\n        D | dist_id, params, pwindow,\n        dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),\n        positive_infinity(), primary_id, primary_params\n      );\n    }\n\n    real log_normalizer = primarycensored_log_normalizer(log_cdf_D, log_cdf_L, L);\n    return log_diff_exp(log_cdf_upper, log_cdf_lower) - log_normalizer;\n  } else {\n    return log_diff_exp(log_cdf_upper, log_cdf_lower);\n  }\n}"

#### 2.2.1 Including dependencies automatically

Many Stan functions in `primarycensored` depend on other functions. For
example, `primarycensored_lpmf` calls `primarycensored_lcdf`, which in
turn may call analytical or ODE-based implementations. When using these
functions in your own Stan models, you need all the dependencies to be
available.

The `dependencies` argument automatically resolves and includes all
required functions in the correct order (dependencies before the
functions that use them):

``` r

pcd_load_stan_functions("primarycensored_lpmf", dependencies = TRUE)
```

    ## [1] "// Stan functions from primarycensored version 1.6.0\nreal expgrowth_cdf(real x, real xmin, real xmax, real r) {\n  if (x < xmin) {\n    return 0;\n  }\n  if (x > xmax) {\n    return 1;\n  }\n  if (abs(r) < 1e-10) {\n    return (x - xmin) / (xmax - xmin);\n  }\n  return (exp(r * x) - exp(r * xmin)) / (exp(r * xmax) - exp(r * xmin));\n}\nreal expgrowth_lcdf(real x, real xmin, real xmax, real r) {\n  if (x < xmin) {\n    return negative_infinity();\n  }\n  if (x > xmax) {\n    return 0;\n  }\n  return log(expgrowth_cdf(x | xmin, xmax, r));\n}\nvector primary_lcdf_vec(vector p, int primary_id,\n                        array[] real primary_params, data real pwindow) {\n  int N = num_elements(p);\n  vector[N] out;\n  for (i in 1:N) {\n    out[i] = primary_lcdf(p[i] | primary_id, primary_params, pwindow);\n  }\n  return out;\n}\nreal discretestep_lcdf(\n  data real d, vector boundaries, vector pmf,\n  int primary_id, array[] real primary_params, data real pwindow\n) {\n  int K = num_elements(pmf);\n  // Integration support in u = d - p for p in [0, pwindow]. It is not\n  // clipped at 0 so boundaries that start below zero (delays with negative\n  // support) are handled; for non-negative boundaries the per-bin clip to\n  // [boundaries[k], boundaries[k + 1]] below gives the same result.\n  real u_min = d - pwindow;\n  real u_max = d;\n\n  // Structural-zero short-circuit. Below `boundaries[2]` F_step is zero\n  // and the bin-1 contribution carries `cum_before = 0`, so the integral\n  // collapses to 0. Returning `negative_infinity()` directly keeps\n  // `log(0)` off the autodiff tape so downstream `log_diff_exp(a, -inf)`\n  // evaluates cleanly with a zero gradient w.r.t. `pmf`.\n  if (u_max <= boundaries[2]) return negative_infinity();\n\n  // Sub-interval endpoints in u-space, clipped to [u_min, u_max].\n  vector[K] lo = fmax(u_min, head(boundaries, K));\n  vector[K] hi = fmin(u_max, tail(boundaries, K));\n\n  // F_step is right-continuous and on [b_k, b_{k+1}) takes the value\n  // sum_{j < k} pmf[j] (mass before bin k). cumulative_sum(pmf) gives\n  // the mass through and including bin k, so we shift right by one.\n  vector[K] cum_before;\n  cum_before[1] = 0;\n  if (K > 1) cum_before[2:K] = head(cumulative_sum(pmf), K - 1);\n\n  // 0/1 mask drops bins with `hi <= lo` from the reduction without a\n  // branch in the inner expression. Built on `data`-level inputs.\n  vector[K] active;\n  for (k in 1:K) active[k] = hi[k] > lo[k] ? 1 : 0;\n\n  // F_primary at lo/hi via two vectorised calls; one masked subtraction\n  // gives the per-bin difference for the dot product.\n  vector[K] f_lo = primary_lcdf_vec(d - lo, primary_id, primary_params,\n                                    pwindow);\n  vector[K] f_hi = primary_lcdf_vec(d - hi, primary_id, primary_params,\n                                    pwindow);\n  vector[K] f_diff = (exp(f_lo) - exp(f_hi)) .* active;\n\n  real integral = dot_product(cum_before, f_diff);\n\n  // Tail region [boundaries[K+1], u_max]: F_step = 1, contributing\n  // F_primary(d - tail_start) - F_primary(d - u_max).\n  real tail_start = fmax(boundaries[K + 1], u_min);\n  if (tail_start < u_max) {\n    real fp_tail = exp(primary_lcdf(d - tail_start | primary_id,\n                                    primary_params, pwindow));\n    real fp_end = exp(primary_lcdf(d - u_max | primary_id,\n                                   primary_params, pwindow));\n    integral += fp_tail - fp_end;\n  }\n\n  return log(integral);\n}\nvector hazards_to_pmf(vector hazards) {\n  int K = num_elements(hazards);\n  vector[K] log_surv;\n  log_surv[1] = 0;\n  if (K > 1) {\n    log_surv[2:K] = cumulative_sum(log1m(hazards[1:(K - 1)]));\n  }\n  return hazards .* exp(log_surv);\n}\nreal discretehazard_lcdf(\n  data real d, vector boundaries, vector hazards,\n  int primary_id, array[] real primary_params, data real pwindow\n) {\n  return discretestep_lcdf(\n    d | boundaries, hazards_to_pmf(hazards), primary_id, primary_params,\n    pwindow\n  );\n}\nreal pstep_lcdf(real t, vector boundaries, vector pmf) {\n  int K = num_elements(pmf);\n  if (t < boundaries[2]) return negative_infinity();\n  if (t >= boundaries[K + 1]) return 0;\n  // Right-continuous CDF with jumps at the right edges\n  // boundaries[2], ..., boundaries[K + 1]. F(t) = cum_pmf[k] for\n  // t in [boundaries[k + 1], boundaries[k + 2]); equivalently the\n  // largest k with boundaries[k + 1] <= t. Boundary-on-jump cases\n  // (t == boundaries[k + 1]) advance k, matching R's\n  // `findInterval(left.open = FALSE)`.\n  int k = 1;\n  while (k < K && boundaries[k + 2] <= t) k += 1;\n  return log(cumulative_sum(pmf)[k]);\n}\nreal phazard_lcdf(real t, vector boundaries, vector hazards) {\n  return pstep_lcdf(t | boundaries, hazards_to_pmf(hazards));\n}\nint check_for_uniform_terms(int dist_id, int primary_id) {\n  if (primary_id != 1) return 0;\n  return dist_id == 2 || dist_id == 1 || dist_id == 3 || dist_id == 5;\n}\nint check_for_analytical(int dist_id, int primary_id) {\n  // Gamma, Lognormal, Weibull and generalised gamma with a Uniform primary\n  if (check_for_uniform_terms(dist_id, primary_id)) return 1;\n  // Keep this primary list in sync with `primary_lcdf`; see the note above.\n  if (dist_id == 26 || dist_id == 27 || dist_id == 28) {\n    return primary_id == 1 || primary_id == 2;\n  }\n  return 0; // No analytical solution for other combinations\n}\nint check_for_tilt_transform(int dist_id, real xi, array[] real params) {\n  if (dist_id == 4) return params[1] - xi > 0;\n  if (dist_id == 2) return params[2] - xi > 0;\n  if (dist_id == 18) return 1;\n  return 0;\n}\nint check_for_exptilt(int dist_id, int primary_id) {\n  return primary_id == 2 && (dist_id == 2 || dist_id == 4 || dist_id == 18);\n}\nint check_for_analytical_params(int dist_id, array[] real params,\n                                int primary_id, array[] real primary_params,\n                                data real pwindow) {\n  if (check_for_exptilt(dist_id, primary_id)) {\n    return check_for_tilt_transform(dist_id, -primary_params[1], params);\n  }\n  return check_for_analytical(dist_id, primary_id);\n}\nreal primarycensored_uniform_lcdf_from_terms(vector terms_d, vector terms_q,\n                                             data real pwindow) {\n  real log_A = log_sum_exp(terms_d[1], terms_q[2]);\n  real log_B = log_sum_exp(terms_q[1], terms_d[2]);\n  // Deep enough into the lower tail every term underflows together. Both\n  // are then `-inf` and `log_diff_exp` would give NaN, so return the limit\n  // directly.\n  if (log_A == negative_infinity() && log_B == negative_infinity()) {\n    return negative_infinity();\n  }\n  return log_diff_exp(log_A, log_B) - log(pwindow);\n}\nvector primarycensored_gamma_uniform_terms(real t,\n                                           array[] real params) {\n  if (t <= 0) {\n    return rep_vector(negative_infinity(), 2);\n  }\n  real shape = params[1];\n  real rate = params[2];\n  real log_rate = gamma_log_rate(rate);\n  // log E where E = k * theta = shape / rate is the mean of the delay\n  real log_E = log(shape) - log_rate;\n  // Both CDFs share one series or fraction evaluation at log(rate * t)\n  real log_t = log(t);\n  vector[2] log_F_T = gamma_lcdf_logx_pair(log_t + log_rate, shape);\n  return [log_t + log_F_T[1], log_E + log_F_T[2]]';\n}\nreal primarycensored_gamma_uniform_lcdf(data real d, real q,\n                                        array[] real params,\n                                        data real pwindow) {\n  return primarycensored_uniform_lcdf_from_terms(\n    primarycensored_gamma_uniform_terms(d, params),\n    primarycensored_gamma_uniform_terms(q, params), pwindow\n  );\n}\nvector primarycensored_lognormal_uniform_terms(real t,\n                                               array[] real params) {\n  real mu = params[1];\n  real sigma = params[2];\n  real mu_sigma2 = mu + square(sigma);\n  // log E where E = exp(mu + sigma^2/2) is the mean of the delay\n  real log_E = mu + 0.5 * square(sigma);\n  real log_t_F_T = lognormal_lcdf_underflows(t, mu, sigma)\n                   ? negative_infinity()\n                   : log(t) + lognormal_lcdf(t | mu, sigma);\n  real log_E_tF_T = lognormal_lcdf_underflows(t, mu_sigma2, sigma)\n                    ? negative_infinity()\n                    : log_E + lognormal_lcdf(t | mu_sigma2, sigma);\n  return [log_t_F_T, log_E_tF_T]';\n}\nreal primarycensored_lognormal_uniform_lcdf(data real d, real q,\n                                            array[] real params,\n                                            data real pwindow) {\n  return primarycensored_uniform_lcdf_from_terms(\n    primarycensored_lognormal_uniform_terms(d, params),\n    primarycensored_lognormal_uniform_terms(q, params), pwindow\n  );\n}\nreal log_weibull_g(real t, real shape, real scale) {\n  real x = pow(t * inv(scale), shape);\n  real a = 1 + inv(shape);\n  // gamma_lcdf(x | a, 1) is log(gamma_p(a, x)), but reverse-mode gamma_p()\n  // returns zero gradients for x / a > 10\n  // (https://github.com/stan-dev/math/issues/2006). gamma_lcdf() has the\n  // same value and computes its own gradients without that cutoff.\n  return gamma_lcdf(x | a, 1) + lgamma(a);\n}\nvector primarycensored_weibull_uniform_terms(real t,\n                                             array[] real params) {\n  if (t <= 0) {\n    return rep_vector(negative_infinity(), 2);\n  }\n  real shape = params[1];\n  real scale = params[2];\n  return [\n    log(t) + weibull_lcdf(t | shape, scale),\n    log(scale) + log_weibull_g(t, shape, scale)\n  ]';\n}\nreal primarycensored_weibull_uniform_lcdf(data real d, real q,\n                                          array[] real params,\n                                          data real pwindow) {\n  return primarycensored_uniform_lcdf_from_terms(\n    primarycensored_weibull_uniform_terms(d, params),\n    primarycensored_weibull_uniform_terms(q, params), pwindow\n  );\n}\nvector primarycensored_gengamma_uniform_terms(real t,\n                                              array[] real params) {\n  if (t <= 0) {\n    return rep_vector(negative_infinity(), 2);\n  }\n  real shape = params[1];\n  real scale = params[2];\n  real k = params[3];\n  real k_shift = k + inv(shape);\n  real log_E = log(scale) + lgamma(k_shift) - lgamma(k);\n  real log_t = log(t);\n  real log_x = shape * (log_t - log(scale));\n  return [\n    log_t + gamma_lcdf_logx(log_x, k),\n    log_E + gamma_lcdf_logx(log_x, k_shift)\n  ]';\n}\nreal primarycensored_gengamma_uniform_lcdf(data real d, real q,\n                                           array[] real params,\n                                           data real pwindow) {\n  return primarycensored_uniform_lcdf_from_terms(\n    primarycensored_gengamma_uniform_terms(d, params),\n    primarycensored_gengamma_uniform_terms(q, params), pwindow\n  );\n}\nreal primarycensored_log_std_normal_cdf(real z) {\n  if (z >= 0) return log(Phi(z));\n  if (z > -37) return log(0.5 * erfc(-z / sqrt(2)));\n  real inv_z2 = 1 / square(z);\n  real series = 1 - inv_z2 * (1 - 3 * inv_z2 * (1 - 5 * inv_z2 * (\n    1 - 7 * inv_z2 * (1 - 9 * inv_z2 * (1 - 11 * inv_z2 * (\n      1 - 13 * inv_z2)))\n  )));\n  return std_normal_lpdf(z) - log(-z) + log(series);\n}\nvector log_tilt_transform_pair(real t, int dist_id, real xi,\n                               array[] real params) {\n  if (dist_id == 2) {\n    real shape = params[1];\n    real rate = params[2];\n    real log_total = -shape * log1m(xi / rate);\n    if (t <= 0) return [negative_infinity(), log_total]';\n    return log_total + gamma_log_pq_logx(log(t) + log(rate - xi), shape);\n  } else if (dist_id == 4) {\n    real rate = params[1];\n    real tilted_rate = rate - xi;\n    real log_total = -log1m(xi / rate);\n    if (t <= 0) return [negative_infinity(), log_total]';\n    return [\n      log_total + log1m_exp(-tilted_rate * t), log_total - tilted_rate * t\n    ]';\n  } else if (dist_id == 18) {\n    // The upper tail is the lower tail of the reflected normal\n    real mu = params[1];\n    real sigma = params[2];\n    real z = (t - mu - xi * square(sigma)) / sigma;\n    real log_total = xi * mu + 0.5 * square(xi * sigma);\n    return [\n      log_total + primarycensored_log_std_normal_cdf(z),\n      log_total + primarycensored_log_std_normal_cdf(-z)\n    ]';\n  }\n  reject(\"Invalid distribution identifier: \", dist_id);\n}\nreal primarycensored_log_diff_exp(real a, real b) {\n  if (b == negative_infinity()) return a;\n  if (a <= b) return negative_infinity();\n  return log_diff_exp(a, b);\n}\nvector primarycensored_moments_from_partial(real log_t, vector log_m) {\n  real log_g1 = primarycensored_log_diff_exp(log_t + log_m[1], log_m[2]);\n  real log_h = primarycensored_log_diff_exp(log_t + log_m[2], log_m[3]);\n  real log_g2 = primarycensored_log_diff_exp(log_t + log_g1, log_h);\n  real log_a = primarycensored_log_diff_exp(log_t + log_m[3], log_m[4]);\n  real log_b = primarycensored_log_diff_exp(log_t + log_h, log_a);\n  real log_g3 = primarycensored_log_diff_exp(log_t + log_g2, log_b);\n  return [log_g1, log_g2, log_g3]';\n}\nvector primarycensored_gamma_tilt_moments(real t, real shape, real rate) {\n  real log_t = log(t);\n  real log_x = log_t + log(rate);\n  // log_mk is the log of int_0^t u^k f(u) du\n  vector[2] log_p01 = gamma_lcdf_logx_pair(log_x, shape);\n  vector[2] log_p23 = gamma_lcdf_logx_pair(log_x, shape + 2);\n  real log_m0 = log_p01[1];\n  real log_m1 = log(shape) - log(rate) + log_p01[2];\n  real log_m2 = log(shape) + log(shape + 1) - 2 * log(rate) + log_p23[1];\n  real log_m3 = log(shape) + log(shape + 1) + log(shape + 2) - 3 * log(rate)\n                + log_p23[2];\n  return primarycensored_moments_from_partial(\n    log_t, [log_m0, log_m1, log_m2, log_m3]'\n  );\n}\nvector primarycensored_tilt_moments(real t, int dist_id,\n                                    array[] real params) {\n  if (dist_id == 2) {\n    if (t <= 0) return rep_vector(negative_infinity(), 3);\n    return primarycensored_gamma_tilt_moments(t, params[1], params[2]);\n  } else if (dist_id == 4) {\n    if (t <= 0) return rep_vector(negative_infinity(), 3);\n    return primarycensored_gamma_tilt_moments(t, 1, params[1]);\n  } else if (dist_id == 18) {\n    real mu = params[1];\n    real sigma = params[2];\n    real z = (t - mu) / sigma;\n    real log_phi = std_normal_lpdf(z);\n    real log_Phi = primarycensored_log_std_normal_cdf(z);\n    real log_g1;\n    real log_g2;\n    real log_g3;\n    if (z > -1) {\n      // The direct form has a derivative at z = 0, unlike the log form\n      log_g1 = log(exp(log_phi) + z * exp(log_Phi));\n      log_g2 = log((square(z) + 1) * exp(log_Phi) + z * exp(log_phi));\n      log_g3 = log(\n        (square(z) + 2) * exp(log_phi) + z * (square(z) + 3) * exp(log_Phi)\n      );\n    } else {\n      // The terms have opposite signs, so the log form keeps the tail\n      log_g1 = primarycensored_log_diff_exp(log_phi, log(-z) + log_Phi);\n      log_g2 = primarycensored_log_diff_exp(\n        log1p(square(z)) + log_Phi, log(-z) + log_phi\n      );\n      log_g3 = primarycensored_log_diff_exp(\n        log(square(z) + 2) + log_phi, log(-z) + log(square(z) + 3) + log_Phi\n      );\n    }\n    return [\n      log(sigma) + log_g1, 2 * log(sigma) + log_g2, 3 * log(sigma) + log_g3\n    ]';\n  }\n  reject(\"Invalid distribution identifier: \", dist_id);\n}\nint primarycensored_exptilt_form(data real d, int dist_id, real rho,\n                                 data real pwindow) {\n  int positive = dist_has_positive_support(dist_id);\n  if (abs(rho) * pwindow < (positive ? 1e-2 : 1e-3)\n      && abs(rho) * (abs(d) + pwindow) < 0.1) {\n    return 1;\n  }\n  if (positive && d < pwindow && abs(rho) * d < 5e-3) return 2;\n  return 0;\n}\nreal primarycensored_tail_diff(real lower_d, real lower_q, real upper_d,\n                               real upper_q) {\n  // NaN from terms that both underflow is a zero difference\n  if (is_nan(lower_q - lower_d) || lower_q - lower_d <= upper_d - upper_q) {\n    return primarycensored_log_diff_exp(lower_d, lower_q);\n  }\n  return primarycensored_log_diff_exp(upper_q, upper_d);\n}\nreal primarycensored_exptilt_lcdf_from_terms(vector terms_d, vector terms_q,\n                                             data real d, real rho,\n                                             data real pwindow) {\n  real log_diff_f = primarycensored_tail_diff(\n    terms_d[1], terms_q[1], terms_d[2], terms_q[2]\n  );\n  real log_diff_j = primarycensored_tail_diff(\n    terms_d[3], terms_q[3], terms_d[4], terms_q[4]\n  );\n  real log_num;\n  real log_den;\n  if (rho > 0) {\n    log_num = primarycensored_log_diff_exp(rho * d + log_diff_j, log_diff_f);\n    log_den = log_diff_exp(rho * pwindow, 0);\n  } else {\n    log_num = primarycensored_log_diff_exp(log_diff_f, rho * d + log_diff_j);\n    log_den = log1m_exp(rho * pwindow);\n  }\n  // log_sum_exp of two `-inf` has a NaN derivative\n  if (terms_q[1] == negative_infinity() && log_num == negative_infinity()) {\n    return negative_infinity();\n  }\n  // Rounding can put the log CDF above 0\n  return fmin(log_sum_exp(terms_q[1], log_num - log_den), 0);\n}\nreal primarycensored_exptilt_small_window_lcdf_from_terms(\n  vector moments_d, vector moments_q, real rho, data real pwindow\n) {\n  real scale = moments_d[1];\n  if (scale == negative_infinity()) return negative_infinity();\n  real g1_q = exp(moments_q[1] - scale);\n  real g2_d = exp(moments_d[2] - scale);\n  real g2_q = exp(moments_q[2] - scale);\n  real g3_d = exp(moments_d[3] - scale);\n  real g3_q = exp(moments_q[3] - scale);\n  real relative = (1 - g1_q)\n                  + 0.5 * rho * (g2_d - pwindow - g2_q - pwindow * g1_q)\n                  + square(rho) * ((g3_d - g3_q) / 6\n                                   - pwindow * (g2_d + g2_q) / 4\n                                   + square(pwindow) * (1 - g1_q) / 12);\n  if (relative <= 0) return negative_infinity();\n  return fmin(scale + log(relative) - log(pwindow), 0);\n}\nreal primarycensored_exptilt_small_delay_lcdf_from_terms(\n  vector moments_d, real rho, data real pwindow\n) {\n  if (moments_d[1] == negative_infinity()) return negative_infinity();\n  real log_den = rho > 0 ? log_diff_exp(rho * pwindow, 0)\n                         : log1m_exp(rho * pwindow);\n  return fmin(\n    moments_d[1] + log(abs(rho)) - log_den\n    + log1p(0.5 * rho * exp(moments_d[2] - moments_d[1])\n            + square(rho) * exp(moments_d[3] - moments_d[1]) / 6),\n    0\n  );\n}\nreal primarycensored_exptilt_lcdf(data real d, int dist_id,\n                                  array[] real params, data real pwindow,\n                                  real rho) {\n  int positive = dist_has_positive_support(dist_id);\n  if (positive && d <= 0) {\n    return negative_infinity();\n  }\n  if (pwindow == 0) return dist_lcdf(d | params, dist_id);\n  real q = d - pwindow;\n  int form = primarycensored_exptilt_form(d, dist_id, rho, pwindow);\n  if (form == 1) {\n    return primarycensored_exptilt_small_window_lcdf_from_terms(\n      primarycensored_tilt_moments(d, dist_id, params),\n      primarycensored_tilt_moments(q, dist_id, params), rho, pwindow\n    );\n  }\n  if (form == 2) {\n    return primarycensored_exptilt_small_delay_lcdf_from_terms(\n      primarycensored_tilt_moments(d, dist_id, params), rho, pwindow\n    );\n  }\n  return primarycensored_exptilt_lcdf_from_terms(\n    append_row(\n      log_tilt_transform_pair(d, dist_id, 0, params),\n      log_tilt_transform_pair(d, dist_id, -rho, params)\n    ),\n    append_row(\n      log_tilt_transform_pair(q, dist_id, 0, params),\n      log_tilt_transform_pair(q, dist_id, -rho, params)\n    ),\n    d, rho, pwindow\n  );\n}\nreal primarycensored_analytical_lcdf_raw(data real d, int dist_id,\n                                         array[] real params,\n                                         data real pwindow,\n                                         int primary_id,\n                                         array[] real primary_params) {\n  real q = max({d - pwindow, 0});\n\n  if (check_for_exptilt(dist_id, primary_id)) {\n    if (!check_for_tilt_transform(dist_id, -primary_params[1], params)) {\n      reject(\n        \"The tilted delay distribution does not exist for tilt \",\n        primary_params[1], \". Use the numerical path, see \",\n        \"check_for_analytical_params().\"\n      );\n    }\n    return primarycensored_exptilt_lcdf(\n      d | dist_id, params, pwindow, primary_params[1]\n    );\n  }\n  if (dist_id == 2 && primary_id == 1) {\n    return primarycensored_gamma_uniform_lcdf(d | q, params, pwindow);\n  } else if (dist_id == 1 && primary_id == 1) {\n    return primarycensored_lognormal_uniform_lcdf(d | q, params, pwindow);\n  } else if (dist_id == 3 && primary_id == 1) {\n    return primarycensored_weibull_uniform_lcdf(d | q, params, pwindow);\n  } else if (dist_id == 5 && primary_id == 1) {\n    return primarycensored_gengamma_uniform_lcdf(d | q, params, pwindow);\n  } else if (dist_id == 26) {\n    // params = [boundaries (K+1), pmf (K)]; length 2*K + 1.\n    int K = (size(params) - 1) %/% 2;\n    return discretestep_lcdf(\n      d | to_vector(segment(params, 1, K + 1)),\n          to_vector(segment(params, K + 2, K)),\n          primary_id, primary_params, pwindow\n    );\n  } else if (dist_id == 27 || dist_id == 28) {\n    // params = [boundaries (K+1), hazards (K)]; length 2*K + 1. The last\n    // hazard must equal 1. RW (27) and RE (28) only differ in their\n    // prior so they share this likelihood dispatch.\n    int K = (size(params) - 1) %/% 2;\n    return discretehazard_lcdf(\n      d | to_vector(segment(params, 1, K + 1)),\n          to_vector(segment(params, K + 2, K)),\n          primary_id, primary_params, pwindow\n    );\n  }\n  return negative_infinity();\n}\nreal primarycensored_analytical_lcdf(data real d, int dist_id,\n                                           array[] real params,\n                                           data real pwindow, data real L,\n                                           data real D, int primary_id,\n                                           array[] real primary_params) {\n  if (d <= L) return negative_infinity();\n  if (d >= D) return 0;\n\n  real result = primarycensored_analytical_lcdf_raw(\n    d, dist_id, params, pwindow, primary_id, primary_params\n  );\n\n  // Apply truncation normalization, also for a finite negative L for delays\n  // on the reals\n  if (!is_inf(D) || L > 0\n      || (!is_inf(L) && !dist_has_positive_support(dist_id))) {\n    vector[2] bounds = primarycensored_truncation_bounds(\n      L, D, dist_id, params, pwindow, primary_id, primary_params\n    );\n    real log_cdf_L = bounds[1];\n    real log_cdf_D = bounds[2];\n\n    real log_normalizer = primarycensored_log_normalizer(log_cdf_D, log_cdf_L, L);\n    result = primarycensored_apply_truncation(result, log_cdf_L, log_normalizer, L);\n  }\n\n  return result;\n}\nreal primarycensored_analytical_cdf(data real d, int dist_id,\n                                          array[] real params,\n                                          data real pwindow, data real L,\n                                          data real D, int primary_id,\n                                          array[] real primary_params) {\n  return exp(primarycensored_analytical_lcdf(d | dist_id, params, pwindow, L, D, primary_id, primary_params));\n}\nint dist_has_positive_support(data int dist_id) {\n  if (dist_id == 1) return 1;   // Lognormal\n  if (dist_id == 2) return 1;   // Gamma\n  if (dist_id == 3) return 1;   // Weibull\n  if (dist_id == 4) return 1;   // Exponential\n  if (dist_id == 5) return 1;   // Generalised gamma\n  if (dist_id == 9) return 1;   // Beta (support on [0, 1])\n  if (dist_id == 13) return 1;  // Chi-square\n  if (dist_id == 16) return 1;  // Inverse Gamma\n  if (dist_id == 19) return 1;  // Inverse Chi-square\n  if (dist_id == 21) return 1;  // Pareto\n  if (dist_id == 22) return 1;  // Scaled inverse Chi-square\n  return 0;\n}\nreal primary_lcdf(real p, int primary_id, array[] real primary_params,\n                  data real pwindow) {\n  if (primary_id == 1) {\n    // Uniform on [0, pwindow]: built-in uniform_lcdf matches the package\n    // primary semantics over [0, pwindow].\n    if (p <= 0) return negative_infinity();\n    if (p >= pwindow) return 0;\n    return uniform_lcdf(p | 0, pwindow);\n  } else if (primary_id == 2) {\n    return expgrowth_lcdf(p | 0, pwindow, primary_params[1]);\n  }\n  reject(\"primary_lcdf: unsupported primary_id \", primary_id);\n}\nreal gamma_log_lead_logx(real log_x, real a) {\n  real x = exp(log_x);\n  if (a < 100) {\n    return a * log_x - x - lgamma(a + 1);\n  }\n  real inv_a = inv(a);\n  real inv_a2 = square(inv_a);\n  real delta = inv_a * (1.0 / 12\n                        - inv_a2 * (1.0 / 360\n                                    - inv_a2 * (1.0 / 1260\n                                                - inv_a2 / 1680)));\n  real r = x / a;\n  // r - 1 is exact near 1, so log1p avoids cancellation\n  real log_r = r > 0.5 && r < 2 ? log1p(r - 1) : log_x - log(a);\n  return a * (log_r - (r - 1))\n         - 0.5 * (log(a) + log(2 * pi())) - delta;\n}\nvector gamma_log_pq_logx(real log_x, real a) {\n  if (!(a > 0) || is_inf(a)) {\n    reject(\"gamma_log_pq_logx: shape must be positive finite, found a = \", a);\n  }\n  real x = exp(log_x);\n  if (log_x == negative_infinity()) {\n    return [negative_infinity(), 0]';\n  }\n  if (is_inf(x)) {\n    return [0, negative_infinity()]';\n  }\n  real log_lead = gamma_log_lead_logx(log_x, a);\n  // Both need about 9 sqrt(a) terms at worst, near x = a + 1\n  real max_terms = 20 * sqrt(a) + 200;\n  if (x < a + 1) {\n    real term = 1;\n    real total = 1;\n    int n = 0;\n    while (term >= 1e-17 * total) {\n      n += 1;\n      if (n > max_terms) {\n        reject(\"gamma_log_pq_logx: series did not converge for log_x = \",\n               log_x, \" and a = \", a);\n      }\n      term *= x / (a + n);\n      total += term;\n    }\n    real log_p = log_lead + log(total);\n    return [log_p, log1m_exp(log_p)]';\n  }\n  real b = x + 1 - a;\n  real c = 1e300;\n  real d = 1 / b;\n  real h = d;\n  real delta = 1;\n  // err tracks |delta - 1|, which each step multiplies by\n  // |an * d * d_prev / delta_prev|. Flooring |i - a| in an at 1 keeps it\n  // from reaching 0 at an integer shape, where the fraction terminates but\n  // its shape gradient does not.\n  real err = 1;\n  int i = 0;\n  while (err >= 1e-15) {\n    i += 1;\n    if (i > max_terms) {\n      reject(\"gamma_log_pq_logx: continued fraction did not converge for \",\n             \"log_x = \", log_x, \" and a = \", a);\n    }\n    real an = -i * (i - a);\n    real d_prev = d;\n    b += 2;\n    d = an * d + b;\n    if (abs(d) < 1e-300) d = 1e-300;\n    c = b + an / c;\n    if (abs(c) < 1e-300) c = 1e-300;\n    d = 1 / d;\n    err *= i * fmax(abs(i - a), 1) * abs(d * d_prev / delta);\n    delta = d * c;\n    h *= delta;\n  }\n  real log_q = log_lead + log(a) + log(h);\n  return [log1m_exp(log_q), log_q]';\n}\nreal gamma_lcdf_logx(real log_x, real a) {\n  return gamma_log_pq_logx(log_x, a)[1];\n}\nvector gamma_lcdf_logx_pair(real log_x, real a) {\n  if (is_inf(log_x)) {\n    real log_p = gamma_lcdf_logx(log_x, a);\n    return [log_p, log_p]';\n  }\n  real log_lead = gamma_log_lead_logx(log_x, a);\n  if (exp(log_x) < a + 1) {\n    real log_p_ap1 = gamma_lcdf_logx(log_x, a + 1);\n    return [log_sum_exp(log_p_ap1, log_lead), log_p_ap1]';\n  }\n  real log_q = gamma_log_pq_logx(log_x, a)[2];\n  return [log1m_exp(log_q), log1m_exp(log_sum_exp(log_q, log_lead))]';\n}\nreal gamma_log_rate(real rate) {\n  if (!(rate > 0) || is_inf(rate)) {\n    reject(\"Gamma rate must be positive finite, found \", rate);\n  }\n  return log(rate);\n}\nint lognormal_lcdf_underflows(real y, real mu, real sigma) {\n  if (y <= 0) {\n    return 1;\n  }\n  return (log(y) - mu) / sigma < -38 ? 1 : 0;\n}\nreal gengamma_lcdf(real y, real shape, real scale, real k) {\n  if (y < 0) {\n    reject(\"gengamma_lcdf: y must be non-negative, found y = \", y);\n  }\n  if (y == 0) {\n    return negative_infinity();\n  }\n  return gamma_lcdf_logx(shape * (log(y) - log(scale)), k);\n}\nreal dist_lcdf(real delay, array[] real params, int dist_id) {\n  if (dist_has_positive_support(dist_id) && delay <= 0) {\n    return negative_infinity();\n  }\n\n  // IDs match pcd_distributions$stan_id in R\n  // Guarded so a lower-tail underflow cannot put a NaN partial on the tape.\n  // The downstream `exp(-inf)` differentiates to 0.\n  if (dist_id == 1) {\n    return lognormal_lcdf_underflows(delay, params[1], params[2])\n           ? negative_infinity()\n           : lognormal_lcdf(delay | params[1], params[2]);\n  }\n  else if (dist_id == 2) {\n    return gamma_lcdf_logx(log(delay) + gamma_log_rate(params[2]), params[1]);\n  }\n  else if (dist_id == 3) return weibull_lcdf(delay | params[1], params[2]);\n  else if (dist_id == 4) return exponential_lcdf(delay | params[1]);\n  else if (dist_id == 5) return gengamma_lcdf(delay | params[1], params[2], params[3]);\n  else if (dist_id == 9) return beta_lcdf(delay | params[1], params[2]);\n  else if (dist_id == 12) return cauchy_lcdf(delay | params[1], params[2]);\n  else if (dist_id == 13) return chi_square_lcdf(delay | params[1]);\n  else if (dist_id == 15) return gumbel_lcdf(delay | params[1], params[2]);\n  else if (dist_id == 16) return inv_gamma_lcdf(delay | params[1], params[2]);\n  else if (dist_id == 17) return logistic_lcdf(delay | params[1], params[2]);\n  else if (dist_id == 18) return normal_lcdf(delay | params[1], params[2]);\n  else if (dist_id == 19) return inv_chi_square_lcdf(delay | params[1]);\n  else if (dist_id == 20) return double_exponential_lcdf(delay | params[1], params[2]);\n  else if (dist_id == 21) return pareto_lcdf(delay | params[1], params[2]);\n  else if (dist_id == 22) return scaled_inv_chi_square_lcdf(delay | params[1], params[2]);\n  else if (dist_id == 23) return student_t_lcdf(delay | params[1], params[2], params[3]);\n  else if (dist_id == 24) return uniform_lcdf(delay | params[1], params[2]);\n  else if (dist_id == 25) return von_mises_lcdf(delay | params[1], params[2]);\n  else if (dist_id == 26) {\n    // Non-parametric step: params = [boundaries (K+1), pmf (K)].\n    int K = (size(params) - 1) %/% 2;\n    return pstep_lcdf(\n      delay | to_vector(segment(params, 1, K + 1)),\n              to_vector(segment(params, K + 2, K))\n    );\n  }\n  else if (dist_id == 27 || dist_id == 28) {\n    // Non-parametric discrete hazard: params = [boundaries (K+1),\n    // hazards (K)] with hazards[K] = 1. RW (27) and RE (28) share the\n    // same likelihood; they only differ in the prior.\n    int K = (size(params) - 1) %/% 2;\n    return phazard_lcdf(\n      delay | to_vector(segment(params, 1, K + 1)),\n              to_vector(segment(params, K + 2, K))\n    );\n  }\n  else reject(\"Invalid distribution identifier: \", dist_id);\n}\nreal primarycensored_log_normalizer(real log_cdf_D, real log_cdf_L, real L) {\n  if (!is_inf(L)) {\n    return log_diff_exp(log_cdf_D, log_cdf_L);\n  } else {\n    return log_cdf_D;\n  }\n}\nreal primarycensored_apply_truncation(real log_cdf, real log_cdf_L,\n                                      real log_normalizer, real L) {\n  if (!is_inf(L)) {\n    return log_diff_exp(log_cdf, log_cdf_L) - log_normalizer;\n  } else {\n    return log_cdf - log_normalizer;\n  }\n}\nvector primarycensored_truncation_bounds(\n  data real L, data real D,\n  data int dist_id, array[] real params, data real pwindow,\n  data int primary_id, array[] real primary_params\n) {\n  vector[2] result;\n  // Internal lower bound for the un-truncated distribution: 0 lets the\n  // `d <= L` early-exit in primarycensored_lcdf return -inf for delays below\n  // the natural support of positive-support distributions; -inf disables that\n  // short-circuit so distributions with support on the reals are integrated.\n  // Expression is inlined (rather than bound to a local) so Stan's data-flow\n  // checker recognises it as data-only.\n\n  // Get CDF at lower truncation point L\n  if (is_inf(L)) {\n    result[1] = negative_infinity();\n  } else {\n    result[1] = primarycensored_lcdf(\n      L | dist_id, params, pwindow,\n      dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),\n      positive_infinity(), primary_id, primary_params\n    );\n  }\n\n  // Get CDF at upper truncation point D\n  if (is_inf(D)) {\n    result[2] = 0;\n  } else {\n    result[2] = primarycensored_lcdf(\n      D | dist_id, params, pwindow,\n      dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),\n      positive_infinity(), primary_id, primary_params\n    );\n  }\n\n  return result;\n}\nreal primarycensored_numeric_cdf(data real d, data int dist_id,\n                                 array[] real params, data real pwindow,\n                                 data int primary_id,\n                                 array[] real primary_params) {\n  real lower_bound = dist_has_positive_support(dist_id)\n                     ? fmax(d - pwindow, 0) : d - pwindow;\n  int n_params = num_elements(params);\n  int n_primary_params = num_elements(primary_params);\n  array[n_params + n_primary_params] real theta = append_array(\n    params, primary_params\n  );\n  array[4] int ids = {dist_id, primary_id, n_params, n_primary_params};\n\n  vector[1] y0 = rep_vector(0.0, 1);\n  return ode_rk45(\n    primarycensored_ode, y0, lower_bound, {d}, theta, {d, pwindow}, ids\n  )[1, 1];\n}\nreal primarycensored_cdf(data real d, data int dist_id, array[] real params,\n                               data real pwindow, data real L, data real D,\n                               data int primary_id,\n                               array[] real primary_params) {\n  real result;\n  if (d <= L) {\n    return 0;\n  }\n\n  if (d >= D) {\n    return 1;\n  }\n\n  // Check if an analytical solution exists and applies for these parameters\n  if (check_for_analytical_params(dist_id, params, primary_id,\n                                  primary_params, pwindow)) {\n    // Use analytical solution\n    result = primarycensored_analytical_cdf(\n      d | dist_id, params, pwindow, L, D, primary_id, primary_params\n    );\n  } else {\n    // Use numerical integration for other cases. For positive-support\n    // delays the integral is 0 for d <= 0.\n    if (dist_has_positive_support(dist_id) && d <= 0) {\n      return 0;\n    }\n    result = primarycensored_numeric_cdf(\n      d | dist_id, params, pwindow, primary_id, primary_params\n    );\n\n    // Apply truncation normalization on log scale for numerical stability.\n    // Skip when F(L) = 0 makes it a no-op (positive support, L <= 0).\n    if (!is_inf(D) || L > 0 ||\n        (!is_inf(L) && !dist_has_positive_support(dist_id))) {\n      real log_result = log(result);\n      vector[2] bounds = primarycensored_truncation_bounds(\n        L, D, dist_id, params, pwindow, primary_id, primary_params\n      );\n      real log_cdf_L = bounds[1];\n      real log_cdf_D = bounds[2];\n\n      real log_normalizer = primarycensored_log_normalizer(log_cdf_D, log_cdf_L, L);\n      log_result = primarycensored_apply_truncation(\n        log_result, log_cdf_L, log_normalizer, L\n      );\n      result = exp(log_result);\n    }\n  }\n\n  return result;\n}\nreal primarycensored_lcdf(data real d, data int dist_id, array[] real params,\n                                data real pwindow, data real L, data real D,\n                                data int primary_id,\n                                array[] real primary_params) {\n  real result;\n\n  if (d <= L) {\n    return negative_infinity();\n  }\n\n  if (d >= D) {\n    return 0;\n  }\n\n  // Check if an analytical solution exists. The internal lower bound is 0 for\n  // positive-support delays (lets the d <= L early-exit return -inf for d <= 0)\n  // and -inf for distributions with support on the reals.\n  if (check_for_analytical_params(dist_id, params, primary_id,\n                                  primary_params, pwindow)) {\n    result = primarycensored_analytical_lcdf(\n      d | dist_id, params, pwindow,\n      dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),\n      positive_infinity(), primary_id, primary_params\n    );\n  } else {\n    // Use numerical integration\n    result = log(primarycensored_cdf(\n      d | dist_id, params, pwindow,\n      dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),\n      positive_infinity(), primary_id, primary_params\n    ));\n  }\n\n  // Handle truncation normalization. Skip when F(L) = 0 makes it a no-op\n  // (positive support, L <= 0) to avoid the cancelling log_diff_exp.\n  if (!is_inf(D) || L > 0 ||\n      (!is_inf(L) && !dist_has_positive_support(dist_id))) {\n    vector[2] bounds = primarycensored_truncation_bounds(\n      L, D, dist_id, params, pwindow, primary_id, primary_params\n    );\n    real log_cdf_L = bounds[1];\n    real log_cdf_D = bounds[2];\n\n    real log_normalizer = primarycensored_log_normalizer(log_cdf_D, log_cdf_L, L);\n    result = primarycensored_apply_truncation(result, log_cdf_L, log_normalizer, L);\n  }\n\n  return result;\n}\nreal primarycensored_lpmf(data int d, data int dist_id, array[] real params,\n                                data real pwindow, data real d_upper,\n                                data real L, data real D, data int primary_id,\n                                array[] real primary_params) {\n  if (d_upper > D) {\n    reject(\"Upper truncation point is greater than D. It is \", d_upper,\n           \" and D is \", D, \". Resolve this by increasing D to be greater or equal to d + swindow or decreasing swindow.\");\n  }\n  if (d_upper <= d) {\n    reject(\"Upper truncation point is less than or equal to d. It is \", d_upper,\n           \" and d is \", d, \". Resolve this by increasing d to be less than d_upper.\");\n  }\n  if (d < L) {\n    return negative_infinity();\n  }\n  real log_cdf_upper = primarycensored_lcdf(\n    d_upper | dist_id, params, pwindow,\n    dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),\n    positive_infinity(), primary_id, primary_params\n  );\n  real log_cdf_lower = primarycensored_lcdf(\n    d | dist_id, params, pwindow,\n    dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),\n    positive_infinity(), primary_id, primary_params\n  );\n\n  // Apply truncation normalization: log((F(d_upper) - F(d)) / (F(D) - F(L))).\n  // Skip when F(L) = 0 makes it a no-op (positive support, L <= 0).\n  if (!is_inf(D) || L > 0 ||\n      (!is_inf(L) && !dist_has_positive_support(dist_id))) {\n    real log_cdf_D;\n    real log_cdf_L;\n\n    // Get CDF at lower truncation point L\n    if (is_inf(L)) {\n      // No left truncation (L = -inf sentinel)\n      log_cdf_L = negative_infinity();\n    } else if (d == L) {\n      // Reuse already computed CDF at d\n      log_cdf_L = log_cdf_lower;\n    } else {\n      // Compute CDF at L directly\n      log_cdf_L = primarycensored_lcdf(\n        L | dist_id, params, pwindow,\n        dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),\n        positive_infinity(), primary_id, primary_params\n      );\n    }\n\n    // Get CDF at upper truncation point D\n    if (d_upper == D) {\n      log_cdf_D = log_cdf_upper;\n    } else if (is_inf(D)) {\n      log_cdf_D = 0;\n    } else {\n      log_cdf_D = primarycensored_lcdf(\n        D | dist_id, params, pwindow,\n        dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),\n        positive_infinity(), primary_id, primary_params\n      );\n    }\n\n    real log_normalizer = primarycensored_log_normalizer(log_cdf_D, log_cdf_L, L);\n    return log_diff_exp(log_cdf_upper, log_cdf_lower) - log_normalizer;\n  } else {\n    return log_diff_exp(log_cdf_upper, log_cdf_lower);\n  }\n}"

This is useful when you want to extract a self-contained set of Stan
functions for use in another project.

#### 2.2.2 Exploring function dependencies

To understand which functions a particular Stan function depends on, use
[`pcd_stan_function_deps()`](https://primarycensored.epinowcast.org/reference/pcd_stan_function_deps.md):

``` r

pcd_stan_function_deps("primarycensored_lpmf")
```

    ##  [1] "dist_has_positive_support"                           
    ##  [2] "primarycensored_log_normalizer"                      
    ##  [3] "check_for_uniform_terms"                             
    ##  [4] "check_for_analytical"                                
    ##  [5] "check_for_tilt_transform"                            
    ##  [6] "check_for_exptilt"                                   
    ##  [7] "check_for_analytical_params"                         
    ##  [8] "expgrowth_cdf"                                       
    ##  [9] "expgrowth_lcdf"                                      
    ## [10] "primary_lcdf"                                        
    ## [11] "primary_lcdf_vec"                                    
    ## [12] "discretestep_lcdf"                                   
    ## [13] "hazards_to_pmf"                                      
    ## [14] "discretehazard_lcdf"                                 
    ## [15] "primarycensored_uniform_lcdf_from_terms"             
    ## [16] "gamma_log_lead_logx"                                 
    ## [17] "gamma_log_pq_logx"                                   
    ## [18] "gamma_lcdf_logx"                                     
    ## [19] "gamma_lcdf_logx_pair"                                
    ## [20] "gamma_log_rate"                                      
    ## [21] "primarycensored_gamma_uniform_terms"                 
    ## [22] "primarycensored_gamma_uniform_lcdf"                  
    ## [23] "lognormal_lcdf_underflows"                           
    ## [24] "primarycensored_lognormal_uniform_terms"             
    ## [25] "primarycensored_lognormal_uniform_lcdf"              
    ## [26] "log_weibull_g"                                       
    ## [27] "primarycensored_weibull_uniform_terms"               
    ## [28] "primarycensored_weibull_uniform_lcdf"                
    ## [29] "primarycensored_gengamma_uniform_terms"              
    ## [30] "primarycensored_gengamma_uniform_lcdf"               
    ## [31] "primarycensored_log_std_normal_cdf"                  
    ## [32] "log_tilt_transform_pair"                             
    ## [33] "primarycensored_log_diff_exp"                        
    ## [34] "primarycensored_moments_from_partial"                
    ## [35] "primarycensored_gamma_tilt_moments"                  
    ## [36] "primarycensored_tilt_moments"                        
    ## [37] "primarycensored_exptilt_form"                        
    ## [38] "primarycensored_tail_diff"                           
    ## [39] "primarycensored_exptilt_lcdf_from_terms"             
    ## [40] "primarycensored_exptilt_small_window_lcdf_from_terms"
    ## [41] "primarycensored_exptilt_small_delay_lcdf_from_terms" 
    ## [42] "pstep_lcdf"                                          
    ## [43] "phazard_lcdf"                                        
    ## [44] "gengamma_lcdf"                                       
    ## [45] "dist_lcdf"                                           
    ## [46] "primarycensored_exptilt_lcdf"                        
    ## [47] "primarycensored_analytical_lcdf_raw"                 
    ## [48] "primarycensored_apply_truncation"                    
    ## [49] "primarycensored_truncation_bounds"                   
    ## [50] "primarycensored_analytical_lcdf"                     
    ## [51] "primarycensored_analytical_cdf"                      
    ## [52] "primarycensored_numeric_cdf"                         
    ## [53] "primarycensored_cdf"                                 
    ## [54] "primarycensored_lcdf"                                
    ## [55] "primarycensored_lpmf"

The result is ordered so that dependencies come before the functions
that use them, with the requested function last. This can help you
understand the structure of the package’s Stan code or decide which
functions you need to include.

### 2.3 Linking the Stan functions to your workflow

#### 2.3.1 Writing functions to a file

One option for using Stan functions is to write them to a file and then
compile them using `cmdstanr`. This is a good approach as it means that
once the functions are written they can be used in the same way as any
other stan functions you might use. The downside is that it may mean
more work keeping up to date with changes to the functions. We can do
this using the
[`pcd_load_stan_functions()`](https://primarycensored.epinowcast.org/reference/pcd_load_stan_functions.md)
function.

``` r

expgrowth_rng_file <- file.path(tempdir(), "expgrowth_rng.stan")
exp_model <- pcd_load_stan_functions(
  "expgrowth_rng",
  write_to_file = TRUE,
  output_file = expgrowth_rng_file,
  wrap_in_block = TRUE
)
```

    ## Stan functions written to: /tmp/Rtmpm0hyXG/expgrowth_rng.stan

This can now be compiled and used in the same way as any other
`cmdstanr` model.

``` r

model <- cmdstan_model(expgrowth_rng_file)
model
```

    ## functions {
    ## // Stan functions from primarycensored version 1.6.0
    ## real expgrowth_rng(real xmin, real xmax, real r) {
    ##   real u = uniform_rng(0, 1);
    ##   if (abs(r) < 1e-10) {
    ##     return xmin + u * (xmax - xmin);
    ##   }
    ##   return log(u * exp(r * xmax) + (1 - u) * exp(r * xmin)) / r;
    ## }
    ## }

Alternatively, you could use `#include expgrowth_rng.stan` in a stan
file functions block to include the function along with the path to that
file as with any other stan file (see
[here](https://mc-stan.org/cmdstanr/reference/model-method-compile.html)).

#### 2.3.2 Including the functions directly via `include_paths`

Rather than writing the functions to a file it is also possible to
include them directly in the stan file using the `include_paths`
argument to
[`cmdstan_model()`](https://mc-stan.org/cmdstanr/reference/cmdstan_model.html).
This is useful if you don’t to clutter your model with the stan code
from `primarycensored` and want automatic updating of the functions. To
demonstrate we will first write a small model has has `expgrowth.stan`
in its include paths (rather than writing it to a file and then
including it). The first step is find the file and path for the
`expgrowth_rng` function.

``` r

pcd_stan_files("expgrowth_rng")
```

    ## [1] "expgrowth.stan"

With that done we now write stan wrapper model.

``` r

expgrowth_stan_file <- file.path(tempdir(), "expgrowth.stan")
writeLines(
  text = c(
    "functions {",
    "#include expgrowth.stan",
    "}",
    "generated quantities {",
    "  real y = expgrowth_rng(0, 1, 0.4);",
    "}"
  ),
  con = expgrowth_stan_file
)
```

We can now use this file to compile a model. **Note** that we need to
include the path to the `primarycensored` Stan functions using the
`include_paths` argument to
[`cmdstan_model()`](https://mc-stan.org/cmdstanr/reference/cmdstan_model.html).

``` r

model <- cmdstan_model(expgrowth_stan_file, include_paths = pcd_stan_path())
model
```

    ## functions {
    ## #include expgrowth.stan
    ## }
    ## generated quantities {
    ##   real y = expgrowth_rng(0, 1, 0.4);
    ## }

We can then sample from the model (we set `fixed_param = TRUE` here as
our toy example doesn’t require MCMC sampling).

``` r

samples <- model$sample(chains = 1, fixed_param = TRUE)
```

    ## Running MCMC with 1 chain...
    ## 
    ## Chain 1 Iteration:   1 / 1000 [  0%]  (Sampling) 
    ## Chain 1 Iteration: 100 / 1000 [ 10%]  (Sampling) 
    ## Chain 1 Iteration: 200 / 1000 [ 20%]  (Sampling) 
    ## Chain 1 Iteration: 300 / 1000 [ 30%]  (Sampling) 
    ## Chain 1 Iteration: 400 / 1000 [ 40%]  (Sampling) 
    ## Chain 1 Iteration: 500 / 1000 [ 50%]  (Sampling) 
    ## Chain 1 Iteration: 600 / 1000 [ 60%]  (Sampling) 
    ## Chain 1 Iteration: 700 / 1000 [ 70%]  (Sampling) 
    ## Chain 1 Iteration: 800 / 1000 [ 80%]  (Sampling) 
    ## Chain 1 Iteration: 900 / 1000 [ 90%]  (Sampling) 
    ## Chain 1 Iteration: 1000 / 1000 [100%]  (Sampling) 
    ## Chain 1 finished in 0.0 seconds.

``` r

samples
```

    ##  variable mean median   sd  mad   q5  q95 rhat ess_bulk ess_tail
    ##         y 0.52   0.54 0.29 0.37 0.06 0.95 1.00      979      802

### 2.4 Using Stan functions directly in R

Whilst it is possible to use Stan functions directly in R it is not
recommended for most use cases (use the R functions in `primarycensored`
instead). However, it can be useful to understand what is going on under
the hood or for exploration (indeed we use this internally in
`primarycensored` to check our functions against the R implementations).
To do this we use the `expose_functions()` method on our already
compiled model. **This can take some time (~30 seconds) to compile all
of the functions.**

``` r

model$expose_functions(global = TRUE)
```

We can now use the function in R. Note that this may get slightly more
complicated if our stan function depends on other stan functions
(i.e. you need to have those included in your compiled model as well).

``` r

expgrowth_rng(0, 1, 0.4)
```

    ## [1] 0.1772395

### 2.5 Summary

In this vignette we have shown approaches for using the Stan functions
provided by `primarycensored` in your Stan modelling workflow. We have
also shown how to use the `expose_functions()` method to access the Stan
functions directly in R for exploration and testing.
