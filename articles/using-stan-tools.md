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
    ## [13] "check_uniform_terms_params"                
    ## [14] "check_for_analytical"                      
    ## [15] "primarycensored_uniform_lcdf_from_terms"   
    ## [16] "primarycensored_gamma_uniform_terms"       
    ## [17] "primarycensored_lognormal_uniform_terms"   
    ## [18] "log_weibull_g"                             
    ## [19] "primarycensored_weibull_uniform_terms"     
    ## [20] "primarycensored_gengamma_uniform_terms"    
    ## [21] "primarycensored_exponential_uniform_terms" 
    ## [22] "primarycensored_normal_uniform_terms"      
    ## [23] "primarycensored_beta_uniform_terms"        
    ## [24] "primarycensored_inv_gamma_lcdf"            
    ## [25] "primarycensored_invgamma_uniform_terms"    
    ## [26] "primarycensored_pareto_uniform_terms"      
    ## [27] "primarycensored_uniform_terms"             
    ## [28] "primarycensored_gamma_uniform_lcdf"        
    ## [29] "primarycensored_lognormal_uniform_lcdf"    
    ## [30] "primarycensored_weibull_uniform_lcdf"      
    ## [31] "primarycensored_gengamma_uniform_lcdf"     
    ## [32] "primarycensored_analytical_lcdf_raw"       
    ## [33] "primarycensored_analytical_lcdf"           
    ## [34] "primarycensored_analytical_cdf"            
    ## [35] "check_for_analytical_vectorized"           
    ## [36] "primarycensored_analytical_lcdf_vectorized"
    ## [37] "gamma_lcdf_logx"                           
    ## [38] "gamma_log_rate"                            
    ## [39] "gengamma_lcdf"                             
    ## [40] "dist_has_positive_support"                 
    ## [41] "lognormal_lcdf_underflows"                 
    ## [42] "dist_lcdf"                                 
    ## [43] "primary_lcdf"                              
    ## [44] "primary_lpdf"                              
    ## [45] "primarycensored_ode"                       
    ## [46] "primarycensored_log_cdf_diff"              
    ## [47] "primarycensored_log_normalizer"            
    ## [48] "primarycensored_apply_truncation"          
    ## [49] "primarycensored_truncation_bounds"         
    ## [50] "primarycensored_cdf"                       
    ## [51] "primarycensored_lcdf"                      
    ## [52] "primarycensored_lpmf"                      
    ## [53] "primarycensored_pmf"                       
    ## [54] "primarycensored_lcdf_vectorized"           
    ## [55] "primarycensored_sone_lpmf_vectorized"      
    ## [56] "primarycensored_sone_pmf_vectorized"

### 2.2 Accessing Stan functions

Stan functions are accessed using the
[`pcd_load_stan_functions()`](https://primarycensored.epinowcast.org/reference/pcd_load_stan_functions.md)
function. This function takes the name of the function as an argument
and returns the function as a string. It can additionally write the
functions to a file and wrap them in a `functions{}` block.

``` r

pcd_load_stan_functions("primarycensored_lpmf")
```

    ## [1] "// Stan functions from primarycensored version 1.6.0\nreal primarycensored_lpmf(data int d, data int dist_id, array[] real params,\n                                data real pwindow, data real d_upper,\n                                data real L, data real D, data int primary_id,\n                                array[] real primary_params) {\n  if (d_upper > D) {\n    reject(\"Upper truncation point is greater than D. It is \", d_upper,\n           \" and D is \", D, \". Resolve this by increasing D to be greater or equal to d + swindow or decreasing swindow.\");\n  }\n  if (d_upper <= d) {\n    reject(\"Upper truncation point is less than or equal to d. It is \", d_upper,\n           \" and d is \", d, \". Resolve this by increasing d to be less than d_upper.\");\n  }\n  if (d < L) {\n    return negative_infinity();\n  }\n  real log_cdf_upper = primarycensored_lcdf(\n    d_upper | dist_id, params, pwindow,\n    dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),\n    positive_infinity(), primary_id, primary_params\n  );\n  real log_cdf_lower = primarycensored_lcdf(\n    d | dist_id, params, pwindow,\n    dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),\n    positive_infinity(), primary_id, primary_params\n  );\n\n  // Apply truncation normalization: log((F(d_upper) - F(d)) / (F(D) - F(L))).\n  // Skip when F(L) = 0 makes it a no-op (positive support, L <= 0).\n  if (!is_inf(D) || L > 0 ||\n      (!is_inf(L) && !dist_has_positive_support(dist_id))) {\n    real log_cdf_D;\n    real log_cdf_L;\n\n    // Get CDF at lower truncation point L\n    if (is_inf(L)) {\n      // No left truncation (L = -inf sentinel)\n      log_cdf_L = negative_infinity();\n    } else if (d == L) {\n      // Reuse already computed CDF at d\n      log_cdf_L = log_cdf_lower;\n    } else {\n      // Compute CDF at L directly\n      log_cdf_L = primarycensored_lcdf(\n        L | dist_id, params, pwindow,\n        dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),\n        positive_infinity(), primary_id, primary_params\n      );\n    }\n\n    // Get CDF at upper truncation point D\n    if (d_upper == D) {\n      log_cdf_D = log_cdf_upper;\n    } else if (is_inf(D)) {\n      log_cdf_D = 0;\n    } else {\n      log_cdf_D = primarycensored_lcdf(\n        D | dist_id, params, pwindow,\n        dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),\n        positive_infinity(), primary_id, primary_params\n      );\n    }\n\n    real log_normalizer = primarycensored_log_normalizer(log_cdf_D, log_cdf_L, L);\n    return primarycensored_log_cdf_diff(log_cdf_upper, log_cdf_lower)\n           - log_normalizer;\n  } else {\n    return primarycensored_log_cdf_diff(log_cdf_upper, log_cdf_lower);\n  }\n}"

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

    ## [1] "// Stan functions from primarycensored version 1.6.0\nreal expgrowth_cdf(real x, real xmin, real xmax, real r) {\n  if (x < xmin) {\n    return 0;\n  }\n  if (x > xmax) {\n    return 1;\n  }\n  if (abs(r) < 1e-10) {\n    return (x - xmin) / (xmax - xmin);\n  }\n  return (exp(r * x) - exp(r * xmin)) / (exp(r * xmax) - exp(r * xmin));\n}\nreal expgrowth_lcdf(real x, real xmin, real xmax, real r) {\n  if (x < xmin) {\n    return negative_infinity();\n  }\n  if (x > xmax) {\n    return 0;\n  }\n  return log(expgrowth_cdf(x | xmin, xmax, r));\n}\nvector primary_lcdf_vec(vector p, int primary_id,\n                        array[] real primary_params, data real pwindow) {\n  int N = num_elements(p);\n  vector[N] out;\n  for (i in 1:N) {\n    out[i] = primary_lcdf(p[i] | primary_id, primary_params, pwindow);\n  }\n  return out;\n}\nreal discretestep_lcdf(\n  data real d, vector boundaries, vector pmf,\n  int primary_id, array[] real primary_params, data real pwindow\n) {\n  int K = num_elements(pmf);\n  // Integration support in u = d - p for p in [0, pwindow]. It is not\n  // clipped at 0 so boundaries that start below zero (delays with negative\n  // support) are handled; for non-negative boundaries the per-bin clip to\n  // [boundaries[k], boundaries[k + 1]] below gives the same result.\n  real u_min = d - pwindow;\n  real u_max = d;\n\n  // Structural-zero short-circuit. Below `boundaries[2]` F_step is zero\n  // and the bin-1 contribution carries `cum_before = 0`, so the integral\n  // collapses to 0. Returning `negative_infinity()` directly keeps\n  // `log(0)` off the autodiff tape so downstream `log_diff_exp(a, -inf)`\n  // evaluates cleanly with a zero gradient w.r.t. `pmf`.\n  if (u_max <= boundaries[2]) return negative_infinity();\n\n  // Sub-interval endpoints in u-space, clipped to [u_min, u_max].\n  vector[K] lo = fmax(u_min, head(boundaries, K));\n  vector[K] hi = fmin(u_max, tail(boundaries, K));\n\n  // F_step is right-continuous and on [b_k, b_{k+1}) takes the value\n  // sum_{j < k} pmf[j] (mass before bin k). cumulative_sum(pmf) gives\n  // the mass through and including bin k, so we shift right by one.\n  vector[K] cum_before;\n  cum_before[1] = 0;\n  if (K > 1) cum_before[2:K] = head(cumulative_sum(pmf), K - 1);\n\n  // 0/1 mask drops bins with `hi <= lo` from the reduction without a\n  // branch in the inner expression. Built on `data`-level inputs.\n  vector[K] active;\n  for (k in 1:K) active[k] = hi[k] > lo[k] ? 1 : 0;\n\n  // F_primary at lo/hi via two vectorised calls; one masked subtraction\n  // gives the per-bin difference for the dot product.\n  vector[K] f_lo = primary_lcdf_vec(d - lo, primary_id, primary_params,\n                                    pwindow);\n  vector[K] f_hi = primary_lcdf_vec(d - hi, primary_id, primary_params,\n                                    pwindow);\n  vector[K] f_diff = (exp(f_lo) - exp(f_hi)) .* active;\n\n  real integral = dot_product(cum_before, f_diff);\n\n  // Tail region [boundaries[K+1], u_max]: F_step = 1, contributing\n  // F_primary(d - tail_start) - F_primary(d - u_max).\n  real tail_start = fmax(boundaries[K + 1], u_min);\n  if (tail_start < u_max) {\n    real fp_tail = exp(primary_lcdf(d - tail_start | primary_id,\n                                    primary_params, pwindow));\n    real fp_end = exp(primary_lcdf(d - u_max | primary_id,\n                                   primary_params, pwindow));\n    integral += fp_tail - fp_end;\n  }\n\n  return log(integral);\n}\nvector hazards_to_pmf(vector hazards) {\n  int K = num_elements(hazards);\n  vector[K] log_surv;\n  log_surv[1] = 0;\n  if (K > 1) {\n    log_surv[2:K] = cumulative_sum(log1m(hazards[1:(K - 1)]));\n  }\n  return hazards .* exp(log_surv);\n}\nreal discretehazard_lcdf(\n  data real d, vector boundaries, vector hazards,\n  int primary_id, array[] real primary_params, data real pwindow\n) {\n  return discretestep_lcdf(\n    d | boundaries, hazards_to_pmf(hazards), primary_id, primary_params,\n    pwindow\n  );\n}\nint check_uniform_terms_params(int dist_id, array[] real params) {\n  if (dist_id == 16) return params[1] > 1;\n  if (dist_id == 19 || dist_id == 22) return params[1] > 2;\n  return 1;\n}\nint check_for_uniform_terms(int dist_id, int primary_id) {\n  if (primary_id != 1) return 0;\n  return dist_id == 1 || dist_id == 2 || dist_id == 3 || dist_id == 4 ||\n    dist_id == 5 || dist_id == 9 || dist_id == 13 || dist_id == 16 ||\n    dist_id == 18 || dist_id == 19 || dist_id == 21 || dist_id == 22;\n}\nint check_for_analytical(int dist_id, int primary_id) {\n  // Delays with uniform primary terms and a Uniform primary, see\n  // check_for_uniform_terms()\n  if (check_for_uniform_terms(dist_id, primary_id)) return 1;\n  // Keep this primary list in sync with `primary_lcdf`; see the note above.\n  if (dist_id == 26 || dist_id == 27 || dist_id == 28) {\n    return primary_id == 1 || primary_id == 2;\n  }\n  return 0; // No analytical solution for other combinations\n}\nreal primarycensored_uniform_lcdf_from_terms(vector terms_d, vector terms_q,\n                                             data real pwindow) {\n  real log_A = log_sum_exp(terms_d[1], terms_q[2]);\n  real log_B = log_sum_exp(terms_q[1], terms_d[2]);\n  // Deep enough into the lower tail every term underflows together. Both\n  // are then `-inf` and `log_diff_exp` would give NaN, so return the limit\n  // directly.\n  if (log_A == negative_infinity() && log_B == negative_infinity()) {\n    return negative_infinity();\n  }\n  real log_cdf = log_diff_exp(log_A, log_B) - log(pwindow);\n  // Rounding can put the upper tail a few ulp above 1. Not fmin, which\n  // would turn a NaN into 0.\n  return log_cdf > 0 ? 0 : log_cdf;\n}\nvector primarycensored_gamma_uniform_terms(real t,\n                                           array[] real params) {\n  if (t <= 0) {\n    return rep_vector(negative_infinity(), 2);\n  }\n  real shape = params[1];\n  real rate = params[2];\n  real log_rate = gamma_log_rate(rate);\n  // log E where E = k * theta = shape / rate is the mean of the delay\n  real log_E = log(shape) - log_rate;\n  real log_t = log(t);\n  real log_y = log_t + log_rate;\n  real log_F_T_kp1 = gamma_lcdf_logx(log_y, shape + 1);\n  real log_lead = shape * log_y - exp(log_y) - lgamma(shape + 1);\n  return [\n    log_t + log_sum_exp(log_F_T_kp1, log_lead),\n    log_E + log_F_T_kp1\n  ]';\n}\nvector primarycensored_lognormal_uniform_terms(real t,\n                                               array[] real params) {\n  real mu = params[1];\n  real sigma = params[2];\n  real mu_sigma2 = mu + square(sigma);\n  // log E where E = exp(mu + sigma^2/2) is the mean of the delay\n  real log_E = mu + 0.5 * square(sigma);\n  real log_t_F_T = lognormal_lcdf_underflows(t, mu, sigma)\n                   ? negative_infinity()\n                   : log(t) + lognormal_lcdf(t | mu, sigma);\n  real log_E_tF_T = lognormal_lcdf_underflows(t, mu_sigma2, sigma)\n                    ? negative_infinity()\n                    : log_E + lognormal_lcdf(t | mu_sigma2, sigma);\n  return [log_t_F_T, log_E_tF_T]';\n}\nreal log_weibull_g(real t, real shape, real scale) {\n  real x = pow(t * inv(scale), shape);\n  real a = 1 + inv(shape);\n  // gamma_lcdf(x | a, 1) is log(gamma_p(a, x)), but reverse-mode gamma_p()\n  // returns zero gradients for x / a > 10\n  // (https://github.com/stan-dev/math/issues/2006). gamma_lcdf() has the\n  // same value and computes its own gradients without that cutoff.\n  return gamma_lcdf(x | a, 1) + lgamma(a);\n}\nvector primarycensored_weibull_uniform_terms(real t,\n                                             array[] real params) {\n  if (t <= 0) {\n    return rep_vector(negative_infinity(), 2);\n  }\n  real shape = params[1];\n  real scale = params[2];\n  return [\n    log(t) + weibull_lcdf(t | shape, scale),\n    log(scale) + log_weibull_g(t, shape, scale)\n  ]';\n}\nvector primarycensored_gengamma_uniform_terms(real t,\n                                              array[] real params) {\n  if (t <= 0) {\n    return rep_vector(negative_infinity(), 2);\n  }\n  real shape = params[1];\n  real scale = params[2];\n  real k = params[3];\n  real k_shift = k + inv(shape);\n  real log_E = log(scale) + lgamma(k_shift) - lgamma(k);\n  return [\n    log(t) + gengamma_lcdf(t | shape, scale, k),\n    log_E + gengamma_lcdf(t | shape, scale, k_shift)\n  ]';\n}\nvector primarycensored_exponential_uniform_terms(real t,\n                                                 array[] real params) {\n  if (t <= 0) {\n    return rep_vector(negative_infinity(), 2);\n  }\n  real rate = params[1];\n  real x = rate * t;\n  real log_x_shortfall;\n  if (x < 0.1) {\n    // x - 1 + exp(-x) = x^2 / 2 * sum_k 2 (-x)^k / (k + 2)!\n    real term = 1;\n    real series = 1;\n    for (k in 1:10) {\n      term *= -x / (k + 2);\n      series += term;\n    }\n    log_x_shortfall = 2 * log(x) - log2() + log(series);\n  } else {\n    log_x_shortfall = log_diff_exp(log(x), log1m_exp(-x));\n  }\n  return [log_x_shortfall - log(rate), negative_infinity()]';\n}\nvector primarycensored_normal_uniform_terms(real t, array[] real params) {\n  real mu = params[1];\n  real sigma = params[2];\n  real z = (t - mu) / sigma;\n  real log_g;\n  if (z < -10) {\n    // z * Phi(z) and phi(z) cancel, so use the asymptotic series\n    real y = inv(square(z));\n    real term = 1;\n    real series = 1;\n    for (k in 1:20) {\n      term *= -(2 * k + 1) * y;\n      series += term;\n    }\n    log_g = std_normal_lpdf(z) + log(y) + log(series);\n  } else {\n    // erfc has an exact autodiff derivative, unlike std_normal_lcdf\n    log_g = log(exp(std_normal_lpdf(z)) + z * 0.5 * erfc(-z / sqrt2()));\n  }\n  return [log(sigma) + log_g, negative_infinity()]';\n}\nvector primarycensored_beta_uniform_terms(real t, array[] real params) {\n  if (t <= 0) {\n    return rep_vector(negative_infinity(), 2);\n  }\n  real a = params[1];\n  real b = params[2];\n  real log_E = log(a) - log(a + b);\n  if (t >= 1) {\n    return [log(t), log_E]';\n  }\n  return [\n    log(t) + beta_lcdf(t | a, b),\n    log_E + beta_lcdf(t | a + 1, b)\n  ]';\n}\nreal primarycensored_inv_gamma_lcdf(real y, real alpha, real beta) {\n  real x = beta / y;\n  // Both need about 9 sqrt(alpha) terms at worst, near x = alpha + 1\n  real max_terms = 20 * sqrt(alpha) + 200;\n  if (x > alpha + 1) {\n    // Within 1e-9 of an integer the fraction terminates, or nearly so, at\n    // i = round(alpha). Its value has then converged but its shape gradient\n    // has not, so the loop runs 50 terms past alpha before it stops.\n    int near_int = abs(alpha - round(alpha)) < 1e-9;\n    if (near_int) max_terms += alpha + 50;\n    real b = x + 1 - alpha;\n    // Modified Lentz starts c at infinity and floors |c| and |d| at 1e-300\n    // to avoid dividing by zero\n    real c = 1e300;\n    real d = inv(b);\n    real h = d;\n    int converged = 0;\n    int i = 0;\n    while (!converged) {\n      i += 1;\n      if (i > max_terms) {\n        reject(\"primarycensored_inv_gamma_lcdf: continued fraction did not \",\n               \"converge for x = \", x, \" and alpha = \", alpha);\n      }\n      real an = -i * (i - alpha);\n      b += 2;\n      d = an * d + b;\n      d = inv(abs(d) < 1e-300 ? 1e-300 : d);\n      c = b + an / c;\n      c = abs(c) < 1e-300 ? 1e-300 : c;\n      h *= d * c;\n      // An update within 1e-15 of 1 leaves h unchanged in double precision\n      converged = abs(d * c - 1) < 1e-15 && (!near_int || i > alpha + 50);\n    }\n    return (alpha - 1) * log(x) - x - lgamma(alpha) + log(x) + log(h);\n  }\n  real term = 1;\n  real series = 1;\n  int n = 0;\n  // A term below 1e-17 of the sum leaves it unchanged in double precision\n  while (term >= 1e-17 * series) {\n    n += 1;\n    if (n > max_terms) {\n      reject(\"primarycensored_inv_gamma_lcdf: series did not converge for \",\n             \"x = \", x, \" and alpha = \", alpha);\n    }\n    term *= x / (alpha + n);\n    series += term;\n  }\n  return log1m_exp(alpha * log(x) - x - lgamma(alpha + 1) + log(series));\n}\nvector primarycensored_invgamma_uniform_terms(real t, array[] real params) {\n  real alpha = params[1];\n  real beta = params[2];\n  if (t <= 0) {\n    return rep_vector(negative_infinity(), 2);\n  }\n  real log_E = log(beta) - log(alpha - 1);\n  return [\n    log(t) + primarycensored_inv_gamma_lcdf(t | alpha, beta),\n    log_E + primarycensored_inv_gamma_lcdf(t | alpha - 1, beta)\n  ]';\n}\nvector primarycensored_pareto_uniform_terms(real t, array[] real params) {\n  real y_min = params[1];\n  real alpha = params[2];\n  if (t <= y_min) {\n    return rep_vector(negative_infinity(), 2);\n  }\n  real L = log(t) - log(y_min);\n  real e = alpha - 1;\n  real log_G_over_y_min;\n  if (L * fmax(1, abs(e)) < 0.1) {\n    // sum_{n >= 2} L^n / n! * (1 + (-1)^n e^(n - 1))\n    real power = L;\n    real e_power = 1;\n    real g = 0;\n    for (n in 2:14) {\n      power *= L / n;\n      e_power *= e;\n      g += power * (1 + (n % 2 == 0 ? e_power : -e_power));\n    }\n    log_G_over_y_min = log(g);\n  } else {\n    real u = -e * L;\n    // expm1(u) / u, with a series near 0 so the derivative is finite\n    real expm1_over_u = abs(u) < 1e-4\n      ? 1 + u * (0.5 + u * (inv(6) + u * inv(24))) : expm1(u) / u;\n    log_G_over_y_min = log(expm1(L) - L * expm1_over_u);\n  }\n  return [log(y_min) + log_G_over_y_min, negative_infinity()]';\n}\nvector primarycensored_uniform_terms(real t, data int dist_id,\n                                     array[] real params) {\n  if (dist_id == 2) {\n    return primarycensored_gamma_uniform_terms(t, params);\n  } else if (dist_id == 1) {\n    return primarycensored_lognormal_uniform_terms(t, params);\n  } else if (dist_id == 3) {\n    return primarycensored_weibull_uniform_terms(t, params);\n  } else if (dist_id == 5) {\n    return primarycensored_gengamma_uniform_terms(t, params);\n  } else if (dist_id == 4) {\n    return primarycensored_exponential_uniform_terms(t, params);\n  } else if (dist_id == 18) {\n    return primarycensored_normal_uniform_terms(t, params);\n  } else if (dist_id == 16) {\n    return primarycensored_invgamma_uniform_terms(t, params);\n  } else if (dist_id == 13) {\n    // Gamma with shape nu / 2 and rate 1 / 2\n    return primarycensored_gamma_uniform_terms(t, {params[1] / 2, 0.5});\n  } else if (dist_id == 19) {\n    // Inverse gamma with shape nu / 2 and scale 1 / 2\n    return primarycensored_invgamma_uniform_terms(\n      t, {params[1] / 2, 0.5}\n    );\n  } else if (dist_id == 22) {\n    // Inverse gamma with shape nu / 2 and scale nu * s^2 / 2\n    return primarycensored_invgamma_uniform_terms(\n      t, {params[1] / 2, params[1] * square(params[2]) / 2}\n    );\n  } else if (dist_id == 21) {\n    return primarycensored_pareto_uniform_terms(t, params);\n  } else if (dist_id == 9) {\n    return primarycensored_beta_uniform_terms(t, params);\n  }\n  reject(\"Invalid distribution identifier: \", dist_id);\n}\nreal primarycensored_analytical_lcdf_raw(data real d, int dist_id,\n                                         array[] real params,\n                                         data real pwindow,\n                                         int primary_id,\n                                         array[] real primary_params) {\n  // Delays on the reals are not clipped at 0\n  real q = dist_has_positive_support(dist_id) ? max({d - pwindow, 0.0})\n                                              : d - pwindow;\n\n  if (check_for_uniform_terms(dist_id, primary_id)) {\n    // The Beta CDF is exactly 1 once the window is above its support\n    if (dist_id == 9 && q >= 1) return 0;\n    return primarycensored_uniform_lcdf_from_terms(\n      primarycensored_uniform_terms(d, dist_id, params),\n      primarycensored_uniform_terms(q, dist_id, params), pwindow\n    );\n  } else if (dist_id == 26) {\n    // params = [boundaries (K+1), pmf (K)]; length 2*K + 1.\n    int K = (size(params) - 1) %/% 2;\n    return discretestep_lcdf(\n      d | to_vector(segment(params, 1, K + 1)),\n          to_vector(segment(params, K + 2, K)),\n          primary_id, primary_params, pwindow\n    );\n  } else if (dist_id == 27 || dist_id == 28) {\n    // params = [boundaries (K+1), hazards (K)]; length 2*K + 1. The last\n    // hazard must equal 1. RW (27) and RE (28) only differ in their\n    // prior so they share this likelihood dispatch.\n    int K = (size(params) - 1) %/% 2;\n    return discretehazard_lcdf(\n      d | to_vector(segment(params, 1, K + 1)),\n          to_vector(segment(params, K + 2, K)),\n          primary_id, primary_params, pwindow\n    );\n  }\n  return negative_infinity();\n}\nreal primarycensored_analytical_cdf(data real d, int dist_id,\n                                          array[] real params,\n                                          data real pwindow, data real L,\n                                          data real D, int primary_id,\n                                          array[] real primary_params) {\n  return exp(primarycensored_analytical_lcdf(d | dist_id, params, pwindow, L, D, primary_id, primary_params));\n}\nreal primarycensored_analytical_lcdf(data real d, int dist_id,\n                                           array[] real params,\n                                           data real pwindow, data real L,\n                                           data real D, int primary_id,\n                                           array[] real primary_params) {\n  if (d <= L) return negative_infinity();\n  if (d >= D) return 0;\n\n  real result = primarycensored_analytical_lcdf_raw(\n    d, dist_id, params, pwindow, primary_id, primary_params\n  );\n\n  // Apply truncation normalization. Skip when F(L) = 0 makes it a no-op\n  // (positive support, L <= 0), as in primarycensored_cdf().\n  if (!is_inf(D) || L > 0 ||\n      (!is_inf(L) && !dist_has_positive_support(dist_id))) {\n    vector[2] bounds = primarycensored_truncation_bounds(\n      L, D, dist_id, params, pwindow, primary_id, primary_params\n    );\n    real log_cdf_L = bounds[1];\n    real log_cdf_D = bounds[2];\n\n    real log_normalizer = primarycensored_log_normalizer(log_cdf_D, log_cdf_L, L);\n    result = primarycensored_apply_truncation(result, log_cdf_L, log_normalizer, L);\n  }\n\n  return result;\n}\nint dist_has_positive_support(data int dist_id) {\n  if (dist_id == 1) return 1;   // Lognormal\n  if (dist_id == 2) return 1;   // Gamma\n  if (dist_id == 3) return 1;   // Weibull\n  if (dist_id == 4) return 1;   // Exponential\n  if (dist_id == 5) return 1;   // Generalised gamma\n  if (dist_id == 9) return 1;   // Beta (support on [0, 1])\n  if (dist_id == 13) return 1;  // Chi-square\n  if (dist_id == 16) return 1;  // Inverse Gamma\n  if (dist_id == 19) return 1;  // Inverse Chi-square\n  if (dist_id == 21) return 1;  // Pareto\n  if (dist_id == 22) return 1;  // Scaled inverse Chi-square\n  return 0;\n}\nreal primary_lcdf(real p, int primary_id, array[] real primary_params,\n                  data real pwindow) {\n  if (primary_id == 1) {\n    // Uniform on [0, pwindow]: built-in uniform_lcdf matches the package\n    // primary semantics over [0, pwindow].\n    if (p <= 0) return negative_infinity();\n    if (p >= pwindow) return 0;\n    return uniform_lcdf(p | 0, pwindow);\n  } else if (primary_id == 2) {\n    return expgrowth_lcdf(p | 0, pwindow, primary_params[1]);\n  }\n  reject(\"primary_lcdf: unsupported primary_id \", primary_id);\n}\nreal gamma_lcdf_logx(real log_x, real a) {\n  if (!(a > 0) || is_inf(a)) {\n    reject(\"gamma_lcdf_logx: shape must be positive finite, found \", a);\n  }\n  real x = exp(log_x);\n  if (!(x < 0.5 * (a + 1))) {\n    return gamma_lcdf(x | a, 1);\n  }\n  // Each term is at most half the one before\n  real term = 1;\n  real total = 1;\n  real n = a;\n  while (term > 1e-16 * total) {\n    n += 1;\n    term *= x / n;\n    total += term;\n  }\n  return a * log_x - x - lgamma(a + 1) + log(total);\n}\nreal gamma_log_rate(real rate) {\n  if (!(rate > 0) || is_inf(rate)) {\n    reject(\"Gamma rate must be positive finite, found \", rate);\n  }\n  return log(rate);\n}\nint lognormal_lcdf_underflows(real y, real mu, real sigma) {\n  if (y <= 0) {\n    return 1;\n  }\n  return (log(y) - mu) / sigma < -38 ? 1 : 0;\n}\nreal gengamma_lcdf(real y, real shape, real scale, real k) {\n  if (y < 0) {\n    reject(\"gengamma_lcdf: y must be non-negative, found y = \", y);\n  }\n  if (y == 0) {\n    return negative_infinity();\n  }\n  return gamma_lcdf_logx(shape * (log(y) - log(scale)), k);\n}\nreal primarycensored_log_cdf_diff(real log_cdf_upper, real log_cdf_lower) {\n  if (log_cdf_upper <= log_cdf_lower) {\n    return negative_infinity();\n  }\n  return log_diff_exp(log_cdf_upper, log_cdf_lower);\n}\nreal primarycensored_log_normalizer(real log_cdf_D, real log_cdf_L, real L) {\n  if (!is_inf(L)) {\n    return log_diff_exp(log_cdf_D, log_cdf_L);\n  } else {\n    return log_cdf_D;\n  }\n}\nreal primarycensored_apply_truncation(real log_cdf, real log_cdf_L,\n                                      real log_normalizer, real L) {\n  if (!is_inf(L)) {\n    return log_diff_exp(log_cdf, log_cdf_L) - log_normalizer;\n  } else {\n    return log_cdf - log_normalizer;\n  }\n}\nvector primarycensored_truncation_bounds(\n  data real L, data real D,\n  data int dist_id, array[] real params, data real pwindow,\n  data int primary_id, array[] real primary_params\n) {\n  vector[2] result;\n  // Internal lower bound for the un-truncated distribution: 0 lets the\n  // `d <= L` early-exit in primarycensored_lcdf return -inf for delays below\n  // the natural support of positive-support distributions; -inf disables that\n  // short-circuit so distributions with support on the reals are integrated.\n  // Expression is inlined (rather than bound to a local) so Stan's data-flow\n  // checker recognises it as data-only.\n\n  // Get CDF at lower truncation point L\n  if (is_inf(L)) {\n    result[1] = negative_infinity();\n  } else {\n    result[1] = primarycensored_lcdf(\n      L | dist_id, params, pwindow,\n      dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),\n      positive_infinity(), primary_id, primary_params\n    );\n  }\n\n  // Get CDF at upper truncation point D\n  if (is_inf(D)) {\n    result[2] = 0;\n  } else {\n    result[2] = primarycensored_lcdf(\n      D | dist_id, params, pwindow,\n      dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),\n      positive_infinity(), primary_id, primary_params\n    );\n  }\n\n  return result;\n}\nreal primarycensored_cdf(data real d, data int dist_id, array[] real params,\n                               data real pwindow, data real L, data real D,\n                               data int primary_id,\n                               array[] real primary_params) {\n  real result;\n  if (d <= L) {\n    return 0;\n  }\n\n  if (d >= D) {\n    return 1;\n  }\n\n  // Check if an analytical solution exists\n  if (check_for_analytical(dist_id, primary_id) &&\n      check_uniform_terms_params(dist_id, params)) {\n    // Use analytical solution\n    result = primarycensored_analytical_cdf(\n      d | dist_id, params, pwindow, L, D, primary_id, primary_params\n    );\n  } else {\n    // Use numerical integration for other cases. The integration variable\n    // ranges over the primary-event time, so the natural lower bound is\n    // d - pwindow. For positive-support delays the integrand `F_delay(t)` is\n    // 0 for t <= 0, so an unclipped lower bound just adds a flat zero region\n    // for negative t. Distributions with support on the reals also accept the\n    // unclipped lower bound directly.\n    real lower_bound = d - pwindow;\n    int n_params = num_elements(params);\n    int n_primary_params = num_elements(primary_params);\n    array[n_params + n_primary_params] real theta = append_array(params, primary_params);\n    array[4] int ids = {dist_id, primary_id, n_params, n_primary_params};\n\n    vector[1] y0 = rep_vector(0.0, 1);\n    result = ode_rk45(primarycensored_ode, y0, lower_bound, {d}, theta, {d, pwindow}, ids)[1, 1];\n\n    // Apply truncation normalization on log scale for numerical stability.\n    // Skip when F(L) = 0 makes it a no-op (positive support, L <= 0).\n    if (!is_inf(D) || L > 0 ||\n        (!is_inf(L) && !dist_has_positive_support(dist_id))) {\n      real log_result = log(result);\n      vector[2] bounds = primarycensored_truncation_bounds(\n        L, D, dist_id, params, pwindow, primary_id, primary_params\n      );\n      real log_cdf_L = bounds[1];\n      real log_cdf_D = bounds[2];\n\n      real log_normalizer = primarycensored_log_normalizer(log_cdf_D, log_cdf_L, L);\n      log_result = primarycensored_apply_truncation(\n        log_result, log_cdf_L, log_normalizer, L\n      );\n      result = exp(log_result);\n    }\n  }\n\n  return result;\n}\nreal primarycensored_lcdf(data real d, data int dist_id, array[] real params,\n                                data real pwindow, data real L, data real D,\n                                data int primary_id,\n                                array[] real primary_params) {\n  real result;\n\n  if (d <= L) {\n    return negative_infinity();\n  }\n\n  if (d >= D) {\n    return 0;\n  }\n\n  // Check if an analytical solution exists. The internal lower bound is 0 for\n  // positive-support delays (lets the d <= L early-exit return -inf for d <= 0)\n  // and -inf for distributions with support on the reals.\n  if (check_for_analytical(dist_id, primary_id) &&\n      check_uniform_terms_params(dist_id, params)) {\n    result = primarycensored_analytical_lcdf(\n      d | dist_id, params, pwindow,\n      dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),\n      positive_infinity(), primary_id, primary_params\n    );\n  } else {\n    // Use numerical integration\n    result = log(primarycensored_cdf(\n      d | dist_id, params, pwindow,\n      dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),\n      positive_infinity(), primary_id, primary_params\n    ));\n  }\n\n  // Handle truncation normalization. Skip when F(L) = 0 makes it a no-op\n  // (positive support, L <= 0) to avoid the cancelling log_diff_exp.\n  if (!is_inf(D) || L > 0 ||\n      (!is_inf(L) && !dist_has_positive_support(dist_id))) {\n    vector[2] bounds = primarycensored_truncation_bounds(\n      L, D, dist_id, params, pwindow, primary_id, primary_params\n    );\n    real log_cdf_L = bounds[1];\n    real log_cdf_D = bounds[2];\n\n    real log_normalizer = primarycensored_log_normalizer(log_cdf_D, log_cdf_L, L);\n    result = primarycensored_apply_truncation(result, log_cdf_L, log_normalizer, L);\n  }\n\n  return result;\n}\nreal primarycensored_lpmf(data int d, data int dist_id, array[] real params,\n                                data real pwindow, data real d_upper,\n                                data real L, data real D, data int primary_id,\n                                array[] real primary_params) {\n  if (d_upper > D) {\n    reject(\"Upper truncation point is greater than D. It is \", d_upper,\n           \" and D is \", D, \". Resolve this by increasing D to be greater or equal to d + swindow or decreasing swindow.\");\n  }\n  if (d_upper <= d) {\n    reject(\"Upper truncation point is less than or equal to d. It is \", d_upper,\n           \" and d is \", d, \". Resolve this by increasing d to be less than d_upper.\");\n  }\n  if (d < L) {\n    return negative_infinity();\n  }\n  real log_cdf_upper = primarycensored_lcdf(\n    d_upper | dist_id, params, pwindow,\n    dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),\n    positive_infinity(), primary_id, primary_params\n  );\n  real log_cdf_lower = primarycensored_lcdf(\n    d | dist_id, params, pwindow,\n    dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),\n    positive_infinity(), primary_id, primary_params\n  );\n\n  // Apply truncation normalization: log((F(d_upper) - F(d)) / (F(D) - F(L))).\n  // Skip when F(L) = 0 makes it a no-op (positive support, L <= 0).\n  if (!is_inf(D) || L > 0 ||\n      (!is_inf(L) && !dist_has_positive_support(dist_id))) {\n    real log_cdf_D;\n    real log_cdf_L;\n\n    // Get CDF at lower truncation point L\n    if (is_inf(L)) {\n      // No left truncation (L = -inf sentinel)\n      log_cdf_L = negative_infinity();\n    } else if (d == L) {\n      // Reuse already computed CDF at d\n      log_cdf_L = log_cdf_lower;\n    } else {\n      // Compute CDF at L directly\n      log_cdf_L = primarycensored_lcdf(\n        L | dist_id, params, pwindow,\n        dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),\n        positive_infinity(), primary_id, primary_params\n      );\n    }\n\n    // Get CDF at upper truncation point D\n    if (d_upper == D) {\n      log_cdf_D = log_cdf_upper;\n    } else if (is_inf(D)) {\n      log_cdf_D = 0;\n    } else {\n      log_cdf_D = primarycensored_lcdf(\n        D | dist_id, params, pwindow,\n        dist_has_positive_support(dist_id) ? 0.0 : negative_infinity(),\n        positive_infinity(), primary_id, primary_params\n      );\n    }\n\n    real log_normalizer = primarycensored_log_normalizer(log_cdf_D, log_cdf_L, L);\n    return primarycensored_log_cdf_diff(log_cdf_upper, log_cdf_lower)\n           - log_normalizer;\n  } else {\n    return primarycensored_log_cdf_diff(log_cdf_upper, log_cdf_lower);\n  }\n}"

This is useful when you want to extract a self-contained set of Stan
functions for use in another project.

#### 2.2.2 Exploring function dependencies

To understand which functions a particular Stan function depends on, use
[`pcd_stan_function_deps()`](https://primarycensored.epinowcast.org/reference/pcd_stan_function_deps.md):

``` r

pcd_stan_function_deps("primarycensored_lpmf")
```

    ##  [1] "dist_has_positive_support"                
    ##  [2] "primarycensored_log_cdf_diff"             
    ##  [3] "primarycensored_log_normalizer"           
    ##  [4] "check_uniform_terms_params"               
    ##  [5] "check_for_uniform_terms"                  
    ##  [6] "check_for_analytical"                     
    ##  [7] "expgrowth_cdf"                            
    ##  [8] "expgrowth_lcdf"                           
    ##  [9] "primary_lcdf"                             
    ## [10] "primary_lcdf_vec"                         
    ## [11] "discretestep_lcdf"                        
    ## [12] "hazards_to_pmf"                           
    ## [13] "discretehazard_lcdf"                      
    ## [14] "primarycensored_uniform_lcdf_from_terms"  
    ## [15] "gamma_lcdf_logx"                          
    ## [16] "gamma_log_rate"                           
    ## [17] "primarycensored_gamma_uniform_terms"      
    ## [18] "lognormal_lcdf_underflows"                
    ## [19] "primarycensored_lognormal_uniform_terms"  
    ## [20] "log_weibull_g"                            
    ## [21] "primarycensored_weibull_uniform_terms"    
    ## [22] "gengamma_lcdf"                            
    ## [23] "primarycensored_gengamma_uniform_terms"   
    ## [24] "primarycensored_exponential_uniform_terms"
    ## [25] "primarycensored_normal_uniform_terms"     
    ## [26] "primarycensored_beta_uniform_terms"       
    ## [27] "primarycensored_inv_gamma_lcdf"           
    ## [28] "primarycensored_invgamma_uniform_terms"   
    ## [29] "primarycensored_pareto_uniform_terms"     
    ## [30] "primarycensored_uniform_terms"            
    ## [31] "primarycensored_analytical_lcdf_raw"      
    ## [32] "primarycensored_apply_truncation"         
    ## [33] "primarycensored_truncation_bounds"        
    ## [34] "primarycensored_analytical_cdf"           
    ## [35] "primarycensored_cdf"                      
    ## [36] "primarycensored_analytical_lcdf"          
    ## [37] "primarycensored_lcdf"                     
    ## [38] "primarycensored_lpmf"

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

    ## Stan functions written to: /tmp/RtmpR7KkGu/expgrowth_rng.stan

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
