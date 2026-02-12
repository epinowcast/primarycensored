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
    ##  [6] "check_for_analytical"                  
    ##  [7] "primarycensored_gamma_uniform_lcdf"    
    ##  [8] "primarycensored_lognormal_uniform_lcdf"
    ##  [9] "log_weibull_g"                         
    ## [10] "primarycensored_weibull_uniform_lcdf"  
    ## [11] "primarycensored_analytical_lcdf_raw"   
    ## [12] "primarycensored_analytical_lcdf"       
    ## [13] "primarycensored_analytical_cdf"        
    ## [14] "dist_lcdf"                             
    ## [15] "primary_lpdf"                          
    ## [16] "primarycensored_ode"                   
    ## [17] "primarycensored_log_normalizer"        
    ## [18] "primarycensored_apply_truncation"      
    ## [19] "primarycensored_truncation_bounds"     
    ## [20] "primarycensored_cdf"                   
    ## [21] "primarycensored_lcdf"                  
    ## [22] "primarycensored_lpmf"                  
    ## [23] "primarycensored_pmf"                   
    ## [24] "primarycensored_sone_lpmf_vectorized"  
    ## [25] "primarycensored_sone_pmf_vectorized"

### 2.2 Accessing Stan functions

Stan functions are accessed using the
[`pcd_load_stan_functions()`](https://primarycensored.epinowcast.org/dev/reference/pcd_load_stan_functions.md)
function. This function takes the name of the function as an argument
and returns the function as a string. It can additionally write the
functions to a file and wrap them in a `functions{}` block.

``` r
pcd_load_stan_functions("primarycensored_lpmf")
```

    ## [1] "// Stan functions from primarycensored version 1.3.0.9000\nreal primarycensored_lpmf(data int d, int dist_id, array[] real params,\n                                data real pwindow, data real d_upper,\n                                data real L, data real D, int primary_id,\n                                array[] real primary_params) {\n  if (d_upper > D) {\n    reject(\"Upper truncation point is greater than D. It is \", d_upper,\n           \" and D is \", D, \". Resolve this by increasing D to be greater or equal to d + swindow or decreasing swindow.\");\n  }\n  if (d_upper <= d) {\n    reject(\"Upper truncation point is less than or equal to d. It is \", d_upper,\n           \" and d is \", d, \". Resolve this by increasing d to be less than d_upper.\");\n  }\n  if (d < L) {\n    return negative_infinity();\n  }\n  real log_cdf_upper = primarycensored_lcdf(\n    d_upper | dist_id, params, pwindow, 0, positive_infinity(), primary_id, primary_params\n  );\n  real log_cdf_lower = primarycensored_lcdf(\n    d | dist_id, params, pwindow, 0, positive_infinity(), primary_id, primary_params\n  );\n\n  // Apply truncation normalization: log((F(d_upper) - F(d)) / (F(D) - F(L)))\n  if (!is_inf(D) || L > 0) {\n    real log_cdf_D;\n    real log_cdf_L;\n\n    // Get CDF at lower truncation point L\n    if (L <= 0) {\n      // No left truncation\n      log_cdf_L = negative_infinity();\n    } else if (d == L) {\n      // Reuse already computed CDF at d\n      log_cdf_L = log_cdf_lower;\n    } else {\n      // Compute CDF at L directly\n      log_cdf_L = primarycensored_lcdf(\n        L | dist_id, params, pwindow, 0, positive_infinity(),\n        primary_id, primary_params\n      );\n    }\n\n    // Get CDF at upper truncation point D\n    if (d_upper == D) {\n      log_cdf_D = log_cdf_upper;\n    } else if (is_inf(D)) {\n      log_cdf_D = 0;\n    } else {\n      log_cdf_D = primarycensored_lcdf(\n        D | dist_id, params, pwindow, 0, positive_infinity(),\n        primary_id, primary_params\n      );\n    }\n\n    real log_normalizer = primarycensored_log_normalizer(log_cdf_D, log_cdf_L, L);\n    return log_diff_exp(log_cdf_upper, log_cdf_lower) - log_normalizer;\n  } else {\n    return log_diff_exp(log_cdf_upper, log_cdf_lower);\n  }\n}"

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

    ## [1] "// Stan functions from primarycensored version 1.3.0.9000\nint check_for_analytical(int dist_id, int primary_id) {\n  if (dist_id == 2 && primary_id == 1) return 1; // Gamma delay with Uniform primary\n  if (dist_id == 1 && primary_id == 1) return 1; // Lognormal delay with Uniform primary\n  if (dist_id == 3 && primary_id == 1) return 1; // Weibull delay with Uniform primary\n  return 0; // No analytical solution for other combinations\n}\nreal primarycensored_gamma_uniform_lcdf(data real d, real q, array[] real params, data real pwindow) {\n  real shape = params[1];\n  real rate = params[2];\n  real shape_1 = shape + 1;\n  real log_window = log(pwindow);\n\n  real log_F_T = gamma_lcdf(d | shape, rate);\n  real log_F_T_kp1 = gamma_lcdf(d | shape_1, rate);\n\n  real log_delta_F_T_kp1;\n  real log_delta_F_T_k;\n  real log_F_Splus;\n\n  if (q != 0) {\n    real log_F_T_q = gamma_lcdf(q | shape, rate);\n    real log_F_T_q_kp1 = gamma_lcdf(q | shape_1, rate);\n\n    // Ensure that the first argument is greater than the second\n    log_delta_F_T_kp1 = log_diff_exp(log_F_T_kp1, log_F_T_q_kp1);\n    log_delta_F_T_k = log_diff_exp(log_F_T, log_F_T_q);\n\n    log_F_Splus = log_diff_exp(\n      log_F_T,\n      log_diff_exp(\n        log(shape * inv(rate)) + log_delta_F_T_kp1,\n        log(d - pwindow) + log_delta_F_T_k\n      ) - log_window\n    );\n  } else {\n    log_delta_F_T_kp1 = log_F_T_kp1;\n    log_delta_F_T_k = log_F_T;\n\n    log_F_Splus = log_diff_exp(\n      log_F_T,\n      log_sum_exp(\n        log(shape * inv(rate)) + log_delta_F_T_kp1,\n        log(pwindow - d) + log_delta_F_T_k\n      ) - log_window\n    );\n  }\n\n  return log_F_Splus;\n}\nreal primarycensored_lognormal_uniform_lcdf(data real d, real q, array[] real params, data real pwindow) {\n  real mu = params[1];\n  real sigma = params[2];\n  real mu_sigma2 = mu + square(sigma);\n  real log_window = log(pwindow);\n\n  real log_F_T = lognormal_lcdf(d | mu, sigma);\n  real log_F_T_mu_sigma2 = lognormal_lcdf(d | mu_sigma2, sigma);\n\n  real log_delta_F_T_mu_sigma;\n  real log_delta_F_T;\n  real log_F_Splus;\n\n  if (q != 0) {\n    real log_F_T_q = lognormal_lcdf(q | mu, sigma);\n    real log_F_T_q_mu_sigma2 = lognormal_lcdf(q | mu_sigma2, sigma);\n\n    // Ensure that the first argument is greater than the second\n    log_delta_F_T_mu_sigma = log_diff_exp(\n      log_F_T_mu_sigma2, log_F_T_q_mu_sigma2\n    );\n    log_delta_F_T = log_diff_exp(log_F_T, log_F_T_q);\n\n    log_F_Splus = log_diff_exp(\n      log_F_T,\n      log_diff_exp(\n        (mu + 0.5 * square(sigma)) + log_delta_F_T_mu_sigma,\n        log(d - pwindow) + log_delta_F_T\n      ) - log_window\n    );\n  } else {\n    log_delta_F_T_mu_sigma = log_F_T_mu_sigma2;\n    log_delta_F_T = log_F_T;\n\n    log_F_Splus = log_diff_exp(\n      log_F_T,\n      log_sum_exp(\n        (mu + 0.5 * square(sigma)) + log_delta_F_T_mu_sigma,\n        log(pwindow - d) + log_delta_F_T\n      ) - log_window\n    );\n  }\n\n  return log_F_Splus;\n}\nreal log_weibull_g(real t, real shape, real scale) {\n  real x = pow(t * inv(scale), shape);\n  real a = 1 + inv(shape);\n  return log(gamma_p(a, x)) + lgamma(a);\n}\nreal primarycensored_weibull_uniform_lcdf(data real d, real q, array[] real params, data real pwindow) {\n  real shape = params[1];\n  real scale = params[2];\n  real log_window = log(pwindow);\n\n  real log_F_T = weibull_lcdf(d | shape, scale);\n\n  real log_delta_g;\n  real log_delta_F_T;\n  real log_F_Splus;\n\n  if (q != 0) {\n    real log_F_T_q = weibull_lcdf(q | shape, scale);\n\n    log_delta_g = log_diff_exp(\n      log_weibull_g(d, shape, scale),\n      log_weibull_g(q, shape, scale)\n    );\n    log_delta_F_T = log_diff_exp(log_F_T, log_F_T_q);\n\n    log_F_Splus = log_diff_exp(\n      log_F_T,\n      log_diff_exp(\n        log(scale) + log_delta_g,\n        log(d - pwindow) + log_delta_F_T\n      ) - log_window\n    );\n  } else {\n    log_delta_g = log_weibull_g(d, shape, scale);\n    log_delta_F_T = log_F_T;\n\n    log_F_Splus = log_diff_exp(\n      log_F_T,\n      log_sum_exp(\n        log(scale) + log_delta_g,\n        log(pwindow - d) + log_delta_F_T\n      ) - log_window\n    );\n  }\n\n  return log_F_Splus;\n}\nreal primarycensored_analytical_lcdf_raw(data real d, int dist_id,\n                                         array[] real params,\n                                         data real pwindow,\n                                         int primary_id) {\n  real q = max({d - pwindow, 0});\n\n  if (dist_id == 2 && primary_id == 1) {\n    return primarycensored_gamma_uniform_lcdf(d | q, params, pwindow);\n  } else if (dist_id == 1 && primary_id == 1) {\n    return primarycensored_lognormal_uniform_lcdf(d | q, params, pwindow);\n  } else if (dist_id == 3 && primary_id == 1) {\n    return primarycensored_weibull_uniform_lcdf(d | q, params, pwindow);\n  }\n  return negative_infinity();\n}\nreal primarycensored_analytical_lcdf(data real d, int dist_id,\n                                           array[] real params,\n                                           data real pwindow, data real L,\n                                           data real D, int primary_id,\n                                           array[] real primary_params) {\n  if (d <= L) return negative_infinity();\n  if (d >= D) return 0;\n\n  real result = primarycensored_analytical_lcdf_raw(\n    d, dist_id, params, pwindow, primary_id\n  );\n\n  // Apply truncation normalization\n  if (!is_inf(D) || L > 0) {\n    vector[2] bounds = primarycensored_truncation_bounds(\n      L, D, dist_id, params, pwindow, primary_id, primary_params\n    );\n    real log_cdf_L = bounds[1];\n    real log_cdf_D = bounds[2];\n\n    real log_normalizer = primarycensored_log_normalizer(log_cdf_D, log_cdf_L, L);\n    result = primarycensored_apply_truncation(result, log_cdf_L, log_normalizer, L);\n  }\n\n  return result;\n}\nreal primarycensored_analytical_cdf(data real d, int dist_id,\n                                          array[] real params,\n                                          data real pwindow, data real L,\n                                          data real D, int primary_id,\n                                          array[] real primary_params) {\n  return exp(primarycensored_analytical_lcdf(d | dist_id, params, pwindow, L, D, primary_id, primary_params));\n}\nreal primarycensored_log_normalizer(real log_cdf_D, real log_cdf_L, real L) {\n  if (L > 0) {\n    return log_diff_exp(log_cdf_D, log_cdf_L);\n  } else {\n    return log_cdf_D;\n  }\n}\nreal primarycensored_apply_truncation(real log_cdf, real log_cdf_L,\n                                      real log_normalizer, real L) {\n  if (L > 0) {\n    return log_diff_exp(log_cdf, log_cdf_L) - log_normalizer;\n  } else {\n    return log_cdf - log_normalizer;\n  }\n}\nvector primarycensored_truncation_bounds(\n  data real L, data real D,\n  int dist_id, array[] real params, data real pwindow,\n  int primary_id, array[] real primary_params\n) {\n  vector[2] result;\n\n  // Get CDF at lower truncation point L\n  if (L <= 0) {\n    result[1] = negative_infinity();\n  } else {\n    result[1] = primarycensored_lcdf(\n      L | dist_id, params, pwindow, 0, positive_infinity(),\n      primary_id, primary_params\n    );\n  }\n\n  // Get CDF at upper truncation point D\n  if (is_inf(D)) {\n    result[2] = 0;\n  } else {\n    result[2] = primarycensored_lcdf(\n      D | dist_id, params, pwindow, 0, positive_infinity(),\n      primary_id, primary_params\n    );\n  }\n\n  return result;\n}\nreal primarycensored_cdf(data real d, int dist_id, array[] real params,\n                               data real pwindow, data real L, data real D,\n                               int primary_id,\n                               array[] real primary_params) {\n  real result;\n  if (d <= L) {\n    return 0;\n  }\n\n  if (d >= D) {\n    return 1;\n  }\n\n  // Check if an analytical solution exists\n  if (check_for_analytical(dist_id, primary_id)) {\n    // Use analytical solution\n    result = primarycensored_analytical_cdf(\n      d | dist_id, params, pwindow, L, D, primary_id, primary_params\n    );\n  } else {\n    // Use numerical integration for other cases\n    real lower_bound = max({d - pwindow, 1e-6});\n    int n_params = num_elements(params);\n    int n_primary_params = num_elements(primary_params);\n    array[n_params + n_primary_params] real theta = append_array(params, primary_params);\n    array[4] int ids = {dist_id, primary_id, n_params, n_primary_params};\n\n    vector[1] y0 = rep_vector(0.0, 1);\n    result = ode_rk45(primarycensored_ode, y0, lower_bound, {d}, theta, {d, pwindow}, ids)[1, 1];\n\n    // Apply truncation normalization on log scale for numerical stability\n    if (!is_inf(D) || L > 0) {\n      real log_result = log(result);\n      vector[2] bounds = primarycensored_truncation_bounds(\n        L, D, dist_id, params, pwindow, primary_id, primary_params\n      );\n      real log_cdf_L = bounds[1];\n      real log_cdf_D = bounds[2];\n\n      real log_normalizer = primarycensored_log_normalizer(log_cdf_D, log_cdf_L, L);\n      log_result = primarycensored_apply_truncation(\n        log_result, log_cdf_L, log_normalizer, L\n      );\n      result = exp(log_result);\n    }\n  }\n\n  return result;\n}\nreal primarycensored_lcdf(data real d, int dist_id, array[] real params,\n                                data real pwindow, data real L, data real D,\n                                int primary_id,\n                                array[] real primary_params) {\n  real result;\n\n  if (d <= L) {\n    return negative_infinity();\n  }\n\n  if (d >= D) {\n    return 0;\n  }\n\n  // Check if an analytical solution exists\n  if (check_for_analytical(dist_id, primary_id)) {\n    result = primarycensored_analytical_lcdf(\n      d | dist_id, params, pwindow, 0, positive_infinity(), primary_id, primary_params\n    );\n  } else {\n    // Use numerical integration\n    result = log(primarycensored_cdf(\n      d | dist_id, params, pwindow, 0, positive_infinity(), primary_id, primary_params\n    ));\n  }\n\n  // Handle truncation normalization\n  if (!is_inf(D) || L > 0) {\n    vector[2] bounds = primarycensored_truncation_bounds(\n      L, D, dist_id, params, pwindow, primary_id, primary_params\n    );\n    real log_cdf_L = bounds[1];\n    real log_cdf_D = bounds[2];\n\n    real log_normalizer = primarycensored_log_normalizer(log_cdf_D, log_cdf_L, L);\n    result = primarycensored_apply_truncation(result, log_cdf_L, log_normalizer, L);\n  }\n\n  return result;\n}\nreal primarycensored_lpmf(data int d, int dist_id, array[] real params,\n                                data real pwindow, data real d_upper,\n                                data real L, data real D, int primary_id,\n                                array[] real primary_params) {\n  if (d_upper > D) {\n    reject(\"Upper truncation point is greater than D. It is \", d_upper,\n           \" and D is \", D, \". Resolve this by increasing D to be greater or equal to d + swindow or decreasing swindow.\");\n  }\n  if (d_upper <= d) {\n    reject(\"Upper truncation point is less than or equal to d. It is \", d_upper,\n           \" and d is \", d, \". Resolve this by increasing d to be less than d_upper.\");\n  }\n  if (d < L) {\n    return negative_infinity();\n  }\n  real log_cdf_upper = primarycensored_lcdf(\n    d_upper | dist_id, params, pwindow, 0, positive_infinity(), primary_id, primary_params\n  );\n  real log_cdf_lower = primarycensored_lcdf(\n    d | dist_id, params, pwindow, 0, positive_infinity(), primary_id, primary_params\n  );\n\n  // Apply truncation normalization: log((F(d_upper) - F(d)) / (F(D) - F(L)))\n  if (!is_inf(D) || L > 0) {\n    real log_cdf_D;\n    real log_cdf_L;\n\n    // Get CDF at lower truncation point L\n    if (L <= 0) {\n      // No left truncation\n      log_cdf_L = negative_infinity();\n    } else if (d == L) {\n      // Reuse already computed CDF at d\n      log_cdf_L = log_cdf_lower;\n    } else {\n      // Compute CDF at L directly\n      log_cdf_L = primarycensored_lcdf(\n        L | dist_id, params, pwindow, 0, positive_infinity(),\n        primary_id, primary_params\n      );\n    }\n\n    // Get CDF at upper truncation point D\n    if (d_upper == D) {\n      log_cdf_D = log_cdf_upper;\n    } else if (is_inf(D)) {\n      log_cdf_D = 0;\n    } else {\n      log_cdf_D = primarycensored_lcdf(\n        D | dist_id, params, pwindow, 0, positive_infinity(),\n        primary_id, primary_params\n      );\n    }\n\n    real log_normalizer = primarycensored_log_normalizer(log_cdf_D, log_cdf_L, L);\n    return log_diff_exp(log_cdf_upper, log_cdf_lower) - log_normalizer;\n  } else {\n    return log_diff_exp(log_cdf_upper, log_cdf_lower);\n  }\n}"

This is useful when you want to extract a self-contained set of Stan
functions for use in another project.

#### 2.2.2 Exploring function dependencies

To understand which functions a particular Stan function depends on, use
[`pcd_stan_function_deps()`](https://primarycensored.epinowcast.org/dev/reference/pcd_stan_function_deps.md):

``` r
pcd_stan_function_deps("primarycensored_lpmf")
```

    ##  [1] "primarycensored_log_normalizer"        
    ##  [2] "check_for_analytical"                  
    ##  [3] "primarycensored_gamma_uniform_lcdf"    
    ##  [4] "primarycensored_lognormal_uniform_lcdf"
    ##  [5] "log_weibull_g"                         
    ##  [6] "primarycensored_weibull_uniform_lcdf"  
    ##  [7] "primarycensored_analytical_lcdf_raw"   
    ##  [8] "primarycensored_apply_truncation"      
    ##  [9] "primarycensored_truncation_bounds"     
    ## [10] "primarycensored_analytical_lcdf"       
    ## [11] "primarycensored_analytical_cdf"        
    ## [12] "primarycensored_cdf"                   
    ## [13] "primarycensored_lcdf"                  
    ## [14] "primarycensored_lpmf"

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

    ## Stan functions written to: /tmp/RtmpdLMOZZ/expgrowth_rng.stan

This can now be compiled and used in the same way as any other
`cmdstanr` model.

``` r
model <- cmdstan_model(expgrowth_rng_file)
model
```

    ## functions {
    ## // Stan functions from primarycensored version 1.3.0.9000
    ## real expgrowth_rng(real xmin, real xmax, real r) {
    ##   real u = uniform_rng(0, 1);
    ##   if (abs(r) < 1e-10) {
    ##     return xmin + u * (xmax - xmin);
    ##   }
    ##   return xmin + log(u * (exp(r * xmax) - exp(r * xmin)) + exp(r * xmin)) / r;
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
