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
[`vignette("fitting-dists-with-stan")`](https://primarycensored.epinowcast.org/dev/articles/fitting-dists-with-stan.md)
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
    ## [12] "check_for_analytical"                  
    ## [13] "primarycensored_gamma_uniform_lcdf"    
    ## [14] "primarycensored_lognormal_uniform_lcdf"
    ## [15] "log_weibull_g"                         
    ## [16] "primarycensored_weibull_uniform_lcdf"  
    ## [17] "primarycensored_gengamma_uniform_lcdf" 
    ## [18] "primarycensored_analytical_lcdf_raw"   
    ## [19] "primarycensored_analytical_lcdf"       
    ## [20] "primarycensored_analytical_cdf"        
    ## [21] "gengamma_lcdf"                         
    ## [22] "dist_has_positive_support"             
    ## [23] "dist_lcdf"                             
    ## [24] "primary_lcdf"                          
    ## [25] "primary_lpdf"                          
    ## [26] "primarycensored_ode"                   
    ## [27] "primarycensored_log_normalizer"        
    ## [28] "primarycensored_apply_truncation"      
    ## [29] "primarycensored_truncation_bounds"     
    ## [30] "primarycensored_cdf"                   
    ## [31] "primarycensored_lcdf"                  
    ## [32] "primarycensored_lpmf"                  
    ## [33] "primarycensored_pmf"                   
    ## [34] "primarycensored_sone_lpmf_vectorized"  
    ## [35] "primarycensored_sone_pmf_vectorized"

### 2.2 Accessing Stan functions

Stan functions are accessed using the
[`pcd_load_stan_functions()`](https://primarycensored.epinowcast.org/dev/reference/pcd_load_stan_functions.md)
function. This function takes the name of the function as an argument
and returns the function as a string. It can additionally write the
functions to a file and wrap them in a `functions{}` block.

``` r

pcd_load_stan_functions("primarycensored_lpmf")
```

    ## [1] "// Stan functions from primarycensored version 1.5.2.1000\nreal primarycensored_lpmf(data int d, data int dist_id, array[] real params,\n                                data real pwindow, data real d_upper,\n                                data real L, data real D, data int primary_id,\n                                array[] real primary_params) {\n  if (d_upper > D) {\n    reject(\"Upper truncation point is greater than D. It is \", d_upper,\n           \" and D is \", D, \". Resolve this by increasing D to be greater or equal to d + swindow or decreasing swindow.\");\n  }\n  if (d_upper <= d) {\n    reject(\"Upper truncation point is less than or equal to d. It is \", d_upper,\n           \" and d is \", d, \". Resolve this by increasing d to be less than d_upper.\");\n  }\n  if (d < L) {\n    return negative_infinity();\n  }\n  real log_cdf_upper = primarycensored_lcdf(\n    d_upper | dist_id, params, pwindow,\n    dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),\n    positive_infinity(), primary_id, primary_params\n  );\n  real log_cdf_lower = primarycensored_lcdf(\n    d | dist_id, params, pwindow,\n    dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),\n    positive_infinity(), primary_id, primary_params\n  );\n\n  // Apply truncation normalization: log((F(d_upper) - F(d)) / (F(D) - F(L))).\n  // Skip when F(L) = 0 makes it a no-op (positive support, L <= 0).\n  if (!is_inf(D) || L > 0 ||\n      (!is_inf(L) && !dist_has_positive_support(dist_id))) {\n    real log_cdf_D;\n    real log_cdf_L;\n\n    // Get CDF at lower truncation point L\n    if (is_inf(L)) {\n      // No left truncation (L = -inf sentinel)\n      log_cdf_L = negative_infinity();\n    } else if (d == L) {\n      // Reuse already computed CDF at d\n      log_cdf_L = log_cdf_lower;\n    } else {\n      // Compute CDF at L directly\n      log_cdf_L = primarycensored_lcdf(\n        L | dist_id, params, pwindow,\n        dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),\n        positive_infinity(), primary_id, primary_params\n      );\n    }\n\n    // Get CDF at upper truncation point D\n    if (d_upper == D) {\n      log_cdf_D = log_cdf_upper;\n    } else if (is_inf(D)) {\n      log_cdf_D = 0;\n    } else {\n      log_cdf_D = primarycensored_lcdf(\n        D | dist_id, params, pwindow,\n        dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),\n        positive_infinity(), primary_id, primary_params\n      );\n    }\n\n    real log_normalizer = primarycensored_log_normalizer(log_cdf_D, log_cdf_L, L);\n    return log_diff_exp(log_cdf_upper, log_cdf_lower) - log_normalizer;\n  } else {\n    return log_diff_exp(log_cdf_upper, log_cdf_lower);\n  }\n}"

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

    ## [1] "// Stan functions from primarycensored version 1.5.2.1000\nreal expgrowth_cdf(real x, real xmin, real xmax, real r) {\n  if (x < xmin) {\n    return 0;\n  }\n  if (x > xmax) {\n    return 1;\n  }\n  if (abs(r) < 1e-10) {\n    return (x - xmin) / (xmax - xmin);\n  }\n  return (exp(r * x) - exp(r * xmin)) / (exp(r * xmax) - exp(r * xmin));\n}\nreal expgrowth_lcdf(real x, real xmin, real xmax, real r) {\n  if (x < xmin) {\n    return negative_infinity();\n  }\n  if (x > xmax) {\n    return 0;\n  }\n  return log(expgrowth_cdf(x | xmin, xmax, r));\n}\nvector primary_lcdf_vec(vector p, int primary_id,\n                        array[] real primary_params, data real pwindow) {\n  int N = num_elements(p);\n  vector[N] out;\n  for (i in 1:N) {\n    out[i] = primary_lcdf(p[i] | primary_id, primary_params, pwindow);\n  }\n  return out;\n}\nreal discretestep_lcdf(\n  data real d, vector boundaries, vector pmf,\n  int primary_id, array[] real primary_params, data real pwindow\n) {\n  int K = num_elements(pmf);\n  // Integration support in u = d - p for p in [0, pwindow]. It is not\n  // clipped at 0 so boundaries that start below zero (delays with negative\n  // support) are handled; for non-negative boundaries the per-bin clip to\n  // [boundaries[k], boundaries[k + 1]] below gives the same result.\n  real u_min = d - pwindow;\n  real u_max = d;\n\n  // Structural-zero short-circuit. Below `boundaries[2]` F_step is zero\n  // and the bin-1 contribution carries `cum_before = 0`, so the integral\n  // collapses to 0. Returning `negative_infinity()` directly keeps\n  // `log(0)` off the autodiff tape so downstream `log_diff_exp(a, -inf)`\n  // evaluates cleanly with a zero gradient w.r.t. `pmf`.\n  if (u_max <= boundaries[2]) return negative_infinity();\n\n  // Sub-interval endpoints in u-space, clipped to [u_min, u_max].\n  vector[K] lo = fmax(u_min, head(boundaries, K));\n  vector[K] hi = fmin(u_max, tail(boundaries, K));\n\n  // F_step is right-continuous and on [b_k, b_{k+1}) takes the value\n  // sum_{j < k} pmf[j] (mass before bin k). cumulative_sum(pmf) gives\n  // the mass through and including bin k, so we shift right by one.\n  vector[K] cum_before;\n  cum_before[1] = 0;\n  if (K > 1) cum_before[2:K] = head(cumulative_sum(pmf), K - 1);\n\n  // 0/1 mask drops bins with `hi <= lo` from the reduction without a\n  // branch in the inner expression. Built on `data`-level inputs.\n  vector[K] active;\n  for (k in 1:K) active[k] = hi[k] > lo[k] ? 1 : 0;\n\n  // F_primary at lo/hi via two vectorised calls; one masked subtraction\n  // gives the per-bin difference for the dot product.\n  vector[K] f_lo = primary_lcdf_vec(d - lo, primary_id, primary_params,\n                                    pwindow);\n  vector[K] f_hi = primary_lcdf_vec(d - hi, primary_id, primary_params,\n                                    pwindow);\n  vector[K] f_diff = (exp(f_lo) - exp(f_hi)) .* active;\n\n  real integral = dot_product(cum_before, f_diff);\n\n  // Tail region [boundaries[K+1], u_max]: F_step = 1, contributing\n  // F_primary(d - tail_start) - F_primary(d - u_max).\n  real tail_start = fmax(boundaries[K + 1], u_min);\n  if (tail_start < u_max) {\n    real fp_tail = exp(primary_lcdf(d - tail_start | primary_id,\n                                    primary_params, pwindow));\n    real fp_end = exp(primary_lcdf(d - u_max | primary_id,\n                                   primary_params, pwindow));\n    integral += fp_tail - fp_end;\n  }\n\n  return log(integral);\n}\nvector hazards_to_pmf(vector hazards) {\n  int K = num_elements(hazards);\n  vector[K] log_surv;\n  log_surv[1] = 0;\n  if (K > 1) {\n    log_surv[2:K] = cumulative_sum(log1m(hazards[1:(K - 1)]));\n  }\n  return hazards .* exp(log_surv);\n}\nreal discretehazard_lcdf(\n  data real d, vector boundaries, vector hazards,\n  int primary_id, array[] real primary_params, data real pwindow\n) {\n  return discretestep_lcdf(\n    d | boundaries, hazards_to_pmf(hazards), primary_id, primary_params,\n    pwindow\n  );\n}\nint check_for_analytical(int dist_id, int primary_id) {\n  if (dist_id == 2 && primary_id == 1) return 1; // Gamma, Uniform\n  if (dist_id == 1 && primary_id == 1) return 1; // Lognormal, Uniform\n  if (dist_id == 3 && primary_id == 1) return 1; // Weibull, Uniform\n  if (dist_id == 5 && primary_id == 1) return 1; // Generalised gamma, Uniform\n  // Keep this primary list in sync with `primary_lcdf`; see the note above.\n  if (dist_id == 26 || dist_id == 27 || dist_id == 28) {\n    return primary_id == 1 || primary_id == 2;\n  }\n  return 0; // No analytical solution for other combinations\n}\nreal primarycensored_gamma_uniform_lcdf(data real d, real q, array[] real params, data real pwindow) {\n  real shape = params[1];\n  real rate = params[2];\n  real log_window = log(pwindow);\n  // log E where E = k * theta = shape / rate is the mean of the delay\n  real log_E = log(shape) - log(rate);\n\n  // F_T(d; k) and the recursion to F_T(d; k+1):\n  // P(k+1, y) = P(k, y) - y^k e^{-y} / Gamma(k+1), with y = rate * d\n  real log_F_T_d_k = gamma_lcdf(d | shape, rate);\n  real gamma_kp1_pdf_log_d\n    = shape * log(rate * d) - rate * d - lgamma(shape + 1);\n  real log_F_T_d_kp1 = log_diff_exp(log_F_T_d_k, gamma_kp1_pdf_log_d);\n\n  // q-dependent terms. Final algebra is unified; only a guard to avoid\n  // log_diff_exp(-inf, -inf) and log(0) when q == 0 (q is data, so autodiff\n  // is unaffected by this branch).\n  real log_q_F_T_q;    // log(q * F_T(q; k))\n  real log_E_tF_T_q;   // log(E * F_T(q; k+1))\n  if (q > 0) {\n    real log_F_T_q_k = gamma_lcdf(q | shape, rate);\n    real gamma_kp1_pdf_log_q\n      = shape * log(rate * q) - rate * q - lgamma(shape + 1);\n    real log_F_T_q_kp1 = log_diff_exp(log_F_T_q_k, gamma_kp1_pdf_log_q);\n    log_q_F_T_q = log(q) + log_F_T_q_k;\n    log_E_tF_T_q = log_E + log_F_T_q_kp1;\n  } else {\n    log_q_F_T_q = negative_infinity();\n    log_E_tF_T_q = negative_infinity();\n  }\n\n  // Unified form: F_{S+}(d) = (A - B) / w_P with A, B sums of positives:\n  //   A = d * F_T(d; k)   + E * F_T(q; k+1)\n  //   B = q * F_T(q; k)   + E * F_T(d; k+1)\n  // Ordering A >= B is guaranteed by F_{S+}(d) >= 0.\n  real log_A = log_sum_exp(log(d) + log_F_T_d_k, log_E_tF_T_q);\n  real log_B = log_sum_exp(log_q_F_T_q, log_E + log_F_T_d_kp1);\n\n  return log_diff_exp(log_A, log_B) - log_window;\n}\nreal primarycensored_lognormal_uniform_lcdf(data real d, real q, array[] real params, data real pwindow) {\n  real mu = params[1];\n  real sigma = params[2];\n  real mu_sigma2 = mu + square(sigma);\n  real log_window = log(pwindow);\n  // log E where E = exp(mu + sigma^2/2) is the mean of the delay\n  real log_E = mu + 0.5 * square(sigma);\n\n  real log_F_T_d = lognormal_lcdf(d | mu, sigma);\n  real log_tF_T_d = lognormal_lcdf(d | mu_sigma2, sigma);\n\n  // q-dependent terms (guard only to avoid log(0); final algebra is unified).\n  real log_q_F_T_q;    // log(q * F_T(q))\n  real log_E_tF_T_q;   // log(E * tilde F_T(q))\n  if (q > 0) {\n    real log_F_T_q = lognormal_lcdf(q | mu, sigma);\n    real log_tF_T_q = lognormal_lcdf(q | mu_sigma2, sigma);\n    log_q_F_T_q = log(q) + log_F_T_q;\n    log_E_tF_T_q = log_E + log_tF_T_q;\n  } else {\n    log_q_F_T_q = negative_infinity();\n    log_E_tF_T_q = negative_infinity();\n  }\n\n  // Unified form: F_{S+}(d) = (A - B) / w_P with\n  //   A = d * F_T(d) + E * tilde F_T(q)\n  //   B = q * F_T(q) + E * tilde F_T(d)\n  // Ordering A >= B is guaranteed by F_{S+}(d) >= 0.\n  real log_A = log_sum_exp(log(d) + log_F_T_d, log_E_tF_T_q);\n  real log_B = log_sum_exp(log_q_F_T_q, log_E + log_tF_T_d);\n\n  return log_diff_exp(log_A, log_B) - log_window;\n}\nreal log_weibull_g(real t, real shape, real scale) {\n  real x = pow(t * inv(scale), shape);\n  real a = 1 + inv(shape);\n  return log(gamma_p(a, x)) + lgamma(a);\n}\nreal primarycensored_weibull_uniform_lcdf(data real d, real q, array[] real params, data real pwindow) {\n  real shape = params[1];\n  real scale = params[2];\n  real log_window = log(pwindow);\n  real log_scale = log(scale);\n\n  // For Weibull: E = scale (lambda) and tilde F_T(t) = g(t; lambda, k), so\n  // log(E * tilde F_T(t)) = log(scale) + log_weibull_g(t, shape, scale).\n  real log_F_T_d = weibull_lcdf(d | shape, scale);\n  real log_E_tF_T_d = log_scale + log_weibull_g(d, shape, scale);\n\n  // q-dependent terms (guard only to avoid log(0); final algebra is unified).\n  real log_q_F_T_q;    // log(q * F_T(q))\n  real log_E_tF_T_q;   // log(E * tilde F_T(q)) = log(scale * g(q; lambda, k))\n  if (q > 0) {\n    log_q_F_T_q = log(q) + weibull_lcdf(q | shape, scale);\n    log_E_tF_T_q = log_scale + log_weibull_g(q, shape, scale);\n  } else {\n    log_q_F_T_q = negative_infinity();\n    log_E_tF_T_q = negative_infinity();\n  }\n\n  // Unified form: F_{S+}(d) = (A - B) / w_P with\n  //   A = d * F_T(d)    + scale * g(q; lambda, k)\n  //   B = q * F_T(q)    + scale * g(d; lambda, k)\n  // Ordering A >= B is guaranteed by F_{S+}(d) >= 0.\n  real log_A = log_sum_exp(log(d) + log_F_T_d, log_E_tF_T_q);\n  real log_B = log_sum_exp(log_q_F_T_q, log_E_tF_T_d);\n\n  return log_diff_exp(log_A, log_B) - log_window;\n}\nreal primarycensored_gengamma_uniform_lcdf(data real d, real q, array[] real params, data real pwindow) {\n  real shape = params[1];\n  real scale = params[2];\n  real k = params[3];\n  real k_shift = k + inv(shape);\n  real log_window = log(pwindow);\n  // log E where E = scale * Gamma(k + 1/shape) / Gamma(k) is the mean of the\n  // delay\n  real log_E = log(scale) + lgamma(k_shift) - lgamma(k);\n\n  real log_F_T_d = gengamma_lcdf(d | shape, scale, k);\n  real log_tF_T_d = gengamma_lcdf(d | shape, scale, k_shift);\n\n  // q-dependent terms (guard only to avoid log(0); final algebra is unified).\n  real log_q_F_T_q;    // log(q * F_T(q))\n  real log_E_tF_T_q;   // log(E * tilde F_T(q))\n  if (q > 0) {\n    log_q_F_T_q = log(q) + gengamma_lcdf(q | shape, scale, k);\n    log_E_tF_T_q = log_E + gengamma_lcdf(q | shape, scale, k_shift);\n  } else {\n    log_q_F_T_q = negative_infinity();\n    log_E_tF_T_q = negative_infinity();\n  }\n\n  // Unified form: F_{S+}(d) = (A - B) / w_P with\n  //   A = d * F_T(d) + E * tilde F_T(q)\n  //   B = q * F_T(q) + E * tilde F_T(d)\n  // Ordering A >= B is guaranteed by F_{S+}(d) >= 0.\n  real log_A = log_sum_exp(log(d) + log_F_T_d, log_E_tF_T_q);\n  real log_B = log_sum_exp(log_q_F_T_q, log_E + log_tF_T_d);\n\n  return log_diff_exp(log_A, log_B) - log_window;\n}\nreal primarycensored_analytical_lcdf_raw(data real d, int dist_id,\n                                         array[] real params,\n                                         data real pwindow,\n                                         int primary_id,\n                                         array[] real primary_params) {\n  real q = max({d - pwindow, 0});\n\n  if (dist_id == 2 && primary_id == 1) {\n    return primarycensored_gamma_uniform_lcdf(d | q, params, pwindow);\n  } else if (dist_id == 1 && primary_id == 1) {\n    return primarycensored_lognormal_uniform_lcdf(d | q, params, pwindow);\n  } else if (dist_id == 3 && primary_id == 1) {\n    return primarycensored_weibull_uniform_lcdf(d | q, params, pwindow);\n  } else if (dist_id == 5 && primary_id == 1) {\n    return primarycensored_gengamma_uniform_lcdf(d | q, params, pwindow);\n  } else if (dist_id == 26) {\n    // params = [boundaries (K+1), pmf (K)]; length 2*K + 1.\n    int K = (size(params) - 1) %/% 2;\n    return discretestep_lcdf(\n      d | to_vector(segment(params, 1, K + 1)),\n          to_vector(segment(params, K + 2, K)),\n          primary_id, primary_params, pwindow\n    );\n  } else if (dist_id == 27 || dist_id == 28) {\n    // params = [boundaries (K+1), hazards (K)]; length 2*K + 1. The last\n    // hazard must equal 1. RW (27) and RE (28) only differ in their\n    // prior so they share this likelihood dispatch.\n    int K = (size(params) - 1) %/% 2;\n    return discretehazard_lcdf(\n      d | to_vector(segment(params, 1, K + 1)),\n          to_vector(segment(params, K + 2, K)),\n          primary_id, primary_params, pwindow\n    );\n  }\n  return negative_infinity();\n}\nreal primarycensored_analytical_lcdf(data real d, int dist_id,\n                                           array[] real params,\n                                           data real pwindow, data real L,\n                                           data real D, int primary_id,\n                                           array[] real primary_params) {\n  if (d <= L) return negative_infinity();\n  if (d >= D) return 0;\n\n  real result = primarycensored_analytical_lcdf_raw(\n    d, dist_id, params, pwindow, primary_id, primary_params\n  );\n\n  // Apply truncation normalization\n  if (!is_inf(D) || L > 0) {\n    vector[2] bounds = primarycensored_truncation_bounds(\n      L, D, dist_id, params, pwindow, primary_id, primary_params\n    );\n    real log_cdf_L = bounds[1];\n    real log_cdf_D = bounds[2];\n\n    real log_normalizer = primarycensored_log_normalizer(log_cdf_D, log_cdf_L, L);\n    result = primarycensored_apply_truncation(result, log_cdf_L, log_normalizer, L);\n  }\n\n  return result;\n}\nreal primarycensored_analytical_cdf(data real d, int dist_id,\n                                          array[] real params,\n                                          data real pwindow, data real L,\n                                          data real D, int primary_id,\n                                          array[] real primary_params) {\n  return exp(primarycensored_analytical_lcdf(d | dist_id, params, pwindow, L, D, primary_id, primary_params));\n}\nint dist_has_positive_support(data int dist_id) {\n  if (dist_id == 1) return 1;   // Lognormal\n  if (dist_id == 2) return 1;   // Gamma\n  if (dist_id == 3) return 1;   // Weibull\n  if (dist_id == 4) return 1;   // Exponential\n  if (dist_id == 5) return 1;   // Generalised gamma\n  if (dist_id == 9) return 1;   // Beta (support on [0, 1])\n  if (dist_id == 13) return 1;  // Chi-square\n  if (dist_id == 16) return 1;  // Inverse Gamma\n  if (dist_id == 19) return 1;  // Inverse Chi-square\n  if (dist_id == 21) return 1;  // Pareto\n  if (dist_id == 22) return 1;  // Scaled inverse Chi-square\n  return 0;\n}\nreal primary_lcdf(real p, int primary_id, array[] real primary_params,\n                  data real pwindow) {\n  if (primary_id == 1) {\n    // Uniform on [0, pwindow]: built-in uniform_lcdf matches the package\n    // primary semantics over [0, pwindow].\n    if (p <= 0) return negative_infinity();\n    if (p >= pwindow) return 0;\n    return uniform_lcdf(p | 0, pwindow);\n  } else if (primary_id == 2) {\n    return expgrowth_lcdf(p | 0, pwindow, primary_params[1]);\n  }\n  reject(\"primary_lcdf: unsupported primary_id \", primary_id);\n}\nreal gengamma_lcdf(real y, real shape, real scale, real k) {\n  return gamma_lcdf(pow(y / scale, shape) | k, 1);\n}\nreal primarycensored_log_normalizer(real log_cdf_D, real log_cdf_L, real L) {\n  if (!is_inf(L)) {\n    return log_diff_exp(log_cdf_D, log_cdf_L);\n  } else {\n    return log_cdf_D;\n  }\n}\nreal primarycensored_apply_truncation(real log_cdf, real log_cdf_L,\n                                      real log_normalizer, real L) {\n  if (!is_inf(L)) {\n    return log_diff_exp(log_cdf, log_cdf_L) - log_normalizer;\n  } else {\n    return log_cdf - log_normalizer;\n  }\n}\nvector primarycensored_truncation_bounds(\n  data real L, data real D,\n  data int dist_id, array[] real params, data real pwindow,\n  data int primary_id, array[] real primary_params\n) {\n  vector[2] result;\n  // Internal lower bound for the un-truncated distribution: 0 lets the\n  // `d <= L` early-exit in primarycensored_lcdf return -inf for delays below\n  // the natural support of positive-support distributions; -inf disables that\n  // short-circuit so distributions with support on the reals are integrated.\n  // Expression is inlined (rather than bound to a local) so Stan's data-flow\n  // checker recognises it as data-only.\n\n  // Get CDF at lower truncation point L\n  if (is_inf(L)) {\n    result[1] = negative_infinity();\n  } else {\n    result[1] = primarycensored_lcdf(\n      L | dist_id, params, pwindow,\n      dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),\n      positive_infinity(), primary_id, primary_params\n    );\n  }\n\n  // Get CDF at upper truncation point D\n  if (is_inf(D)) {\n    result[2] = 0;\n  } else {\n    result[2] = primarycensored_lcdf(\n      D | dist_id, params, pwindow,\n      dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),\n      positive_infinity(), primary_id, primary_params\n    );\n  }\n\n  return result;\n}\nreal primarycensored_cdf(data real d, data int dist_id, array[] real params,\n                               data real pwindow, data real L, data real D,\n                               data int primary_id,\n                               array[] real primary_params) {\n  real result;\n  if (d <= L) {\n    return 0;\n  }\n\n  if (d >= D) {\n    return 1;\n  }\n\n  // Check if an analytical solution exists\n  if (check_for_analytical(dist_id, primary_id)) {\n    // Use analytical solution\n    result = primarycensored_analytical_cdf(\n      d | dist_id, params, pwindow, L, D, primary_id, primary_params\n    );\n  } else {\n    // Use numerical integration for other cases. The integration variable\n    // ranges over the primary-event time, so the natural lower bound is\n    // d - pwindow. For positive-support delays the integrand `F_delay(t)` is\n    // 0 for t <= 0, so an unclipped lower bound just adds a flat zero region\n    // for negative t. Distributions with support on the reals also accept the\n    // unclipped lower bound directly.\n    real lower_bound = d - pwindow;\n    int n_params = num_elements(params);\n    int n_primary_params = num_elements(primary_params);\n    array[n_params + n_primary_params] real theta = append_array(params, primary_params);\n    array[4] int ids = {dist_id, primary_id, n_params, n_primary_params};\n\n    vector[1] y0 = rep_vector(0.0, 1);\n    result = ode_rk45(primarycensored_ode, y0, lower_bound, {d}, theta, {d, pwindow}, ids)[1, 1];\n\n    // Apply truncation normalization on log scale for numerical stability.\n    // Skip when F(L) = 0 makes it a no-op (positive support, L <= 0).\n    if (!is_inf(D) || L > 0 ||\n        (!is_inf(L) && !dist_has_positive_support(dist_id))) {\n      real log_result = log(result);\n      vector[2] bounds = primarycensored_truncation_bounds(\n        L, D, dist_id, params, pwindow, primary_id, primary_params\n      );\n      real log_cdf_L = bounds[1];\n      real log_cdf_D = bounds[2];\n\n      real log_normalizer = primarycensored_log_normalizer(log_cdf_D, log_cdf_L, L);\n      log_result = primarycensored_apply_truncation(\n        log_result, log_cdf_L, log_normalizer, L\n      );\n      result = exp(log_result);\n    }\n  }\n\n  return result;\n}\nreal primarycensored_lcdf(data real d, data int dist_id, array[] real params,\n                                data real pwindow, data real L, data real D,\n                                data int primary_id,\n                                array[] real primary_params) {\n  real result;\n\n  if (d <= L) {\n    return negative_infinity();\n  }\n\n  if (d >= D) {\n    return 0;\n  }\n\n  // Check if an analytical solution exists. The internal lower bound is 0 for\n  // positive-support delays (lets the d <= L early-exit return -inf for d <= 0)\n  // and -inf for distributions with support on the reals.\n  if (check_for_analytical(dist_id, primary_id)) {\n    result = primarycensored_analytical_lcdf(\n      d | dist_id, params, pwindow,\n      dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),\n      positive_infinity(), primary_id, primary_params\n    );\n  } else {\n    // Use numerical integration\n    result = log(primarycensored_cdf(\n      d | dist_id, params, pwindow,\n      dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),\n      positive_infinity(), primary_id, primary_params\n    ));\n  }\n\n  // Handle truncation normalization. Skip when F(L) = 0 makes it a no-op\n  // (positive support, L <= 0) to avoid the cancelling log_diff_exp.\n  if (!is_inf(D) || L > 0 ||\n      (!is_inf(L) && !dist_has_positive_support(dist_id))) {\n    vector[2] bounds = primarycensored_truncation_bounds(\n      L, D, dist_id, params, pwindow, primary_id, primary_params\n    );\n    real log_cdf_L = bounds[1];\n    real log_cdf_D = bounds[2];\n\n    real log_normalizer = primarycensored_log_normalizer(log_cdf_D, log_cdf_L, L);\n    result = primarycensored_apply_truncation(result, log_cdf_L, log_normalizer, L);\n  }\n\n  return result;\n}\nreal primarycensored_lpmf(data int d, data int dist_id, array[] real params,\n                                data real pwindow, data real d_upper,\n                                data real L, data real D, data int primary_id,\n                                array[] real primary_params) {\n  if (d_upper > D) {\n    reject(\"Upper truncation point is greater than D. It is \", d_upper,\n           \" and D is \", D, \". Resolve this by increasing D to be greater or equal to d + swindow or decreasing swindow.\");\n  }\n  if (d_upper <= d) {\n    reject(\"Upper truncation point is less than or equal to d. It is \", d_upper,\n           \" and d is \", d, \". Resolve this by increasing d to be less than d_upper.\");\n  }\n  if (d < L) {\n    return negative_infinity();\n  }\n  real log_cdf_upper = primarycensored_lcdf(\n    d_upper | dist_id, params, pwindow,\n    dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),\n    positive_infinity(), primary_id, primary_params\n  );\n  real log_cdf_lower = primarycensored_lcdf(\n    d | dist_id, params, pwindow,\n    dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),\n    positive_infinity(), primary_id, primary_params\n  );\n\n  // Apply truncation normalization: log((F(d_upper) - F(d)) / (F(D) - F(L))).\n  // Skip when F(L) = 0 makes it a no-op (positive support, L <= 0).\n  if (!is_inf(D) || L > 0 ||\n      (!is_inf(L) && !dist_has_positive_support(dist_id))) {\n    real log_cdf_D;\n    real log_cdf_L;\n\n    // Get CDF at lower truncation point L\n    if (is_inf(L)) {\n      // No left truncation (L = -inf sentinel)\n      log_cdf_L = negative_infinity();\n    } else if (d == L) {\n      // Reuse already computed CDF at d\n      log_cdf_L = log_cdf_lower;\n    } else {\n      // Compute CDF at L directly\n      log_cdf_L = primarycensored_lcdf(\n        L | dist_id, params, pwindow,\n        dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),\n        positive_infinity(), primary_id, primary_params\n      );\n    }\n\n    // Get CDF at upper truncation point D\n    if (d_upper == D) {\n      log_cdf_D = log_cdf_upper;\n    } else if (is_inf(D)) {\n      log_cdf_D = 0;\n    } else {\n      log_cdf_D = primarycensored_lcdf(\n        D | dist_id, params, pwindow,\n        dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),\n        positive_infinity(), primary_id, primary_params\n      );\n    }\n\n    real log_normalizer = primarycensored_log_normalizer(log_cdf_D, log_cdf_L, L);\n    return log_diff_exp(log_cdf_upper, log_cdf_lower) - log_normalizer;\n  } else {\n    return log_diff_exp(log_cdf_upper, log_cdf_lower);\n  }\n}"

This is useful when you want to extract a self-contained set of Stan
functions for use in another project.

#### 2.2.2 Exploring function dependencies

To understand which functions a particular Stan function depends on, use
[`pcd_stan_function_deps()`](https://primarycensored.epinowcast.org/dev/reference/pcd_stan_function_deps.md):

``` r

pcd_stan_function_deps("primarycensored_lpmf")
```

    ##  [1] "dist_has_positive_support"             
    ##  [2] "primarycensored_log_normalizer"        
    ##  [3] "check_for_analytical"                  
    ##  [4] "expgrowth_cdf"                         
    ##  [5] "expgrowth_lcdf"                        
    ##  [6] "primary_lcdf"                          
    ##  [7] "primary_lcdf_vec"                      
    ##  [8] "discretestep_lcdf"                     
    ##  [9] "hazards_to_pmf"                        
    ## [10] "discretehazard_lcdf"                   
    ## [11] "primarycensored_gamma_uniform_lcdf"    
    ## [12] "primarycensored_lognormal_uniform_lcdf"
    ## [13] "log_weibull_g"                         
    ## [14] "primarycensored_weibull_uniform_lcdf"  
    ## [15] "gengamma_lcdf"                         
    ## [16] "primarycensored_gengamma_uniform_lcdf" 
    ## [17] "primarycensored_analytical_lcdf_raw"   
    ## [18] "primarycensored_apply_truncation"      
    ## [19] "primarycensored_truncation_bounds"     
    ## [20] "primarycensored_analytical_lcdf"       
    ## [21] "primarycensored_analytical_cdf"        
    ## [22] "primarycensored_cdf"                   
    ## [23] "primarycensored_lcdf"                  
    ## [24] "primarycensored_lpmf"

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
[`pcd_load_stan_functions()`](https://primarycensored.epinowcast.org/dev/reference/pcd_load_stan_functions.md)
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

    ## Stan functions written to: /tmp/RtmpXmrHF5/expgrowth_rng.stan

This can now be compiled and used in the same way as any other
`cmdstanr` model.

``` r

model <- cmdstan_model(expgrowth_rng_file)
model
```

    ## functions {
    ## // Stan functions from primarycensored version 1.5.2.1000
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
