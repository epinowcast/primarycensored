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
    ##  [6] "check_for_analytical"                  
    ##  [7] "primarycensored_gamma_uniform_lcdf"    
    ##  [8] "primarycensored_lognormal_uniform_lcdf"
    ##  [9] "log_weibull_g"                         
    ## [10] "primarycensored_weibull_uniform_lcdf"  
    ## [11] "primarycensored_analytical_lcdf"       
    ## [12] "primarycensored_analytical_cdf"        
    ## [13] "dist_lcdf"                             
    ## [14] "primary_lpdf"                          
    ## [15] "primarycensored_ode"                   
    ## [16] "primarycensored_cdf"                   
    ## [17] "primarycensored_lcdf"                  
    ## [18] "primarycensored_lpmf"                  
    ## [19] "primarycensored_pmf"                   
    ## [20] "primarycensored_sone_lpmf_vectorized"  
    ## [21] "primarycensored_sone_pmf_vectorized"

### 2.2 Accessing Stan functions

Stan functions are accessed using the
[`pcd_load_stan_functions()`](https://primarycensored.epinowcast.org/reference/pcd_load_stan_functions.md)
function. This function takes the name of the function as an argument
and returns the function as a string. It can additionally write the
functions to a file and wrap them in a `functions{}` block.

``` r
pcd_load_stan_functions("primarycensored_lpmf")
```

    ## [1] "// Stan functions from primarycensored version 1.3.0\nreal primarycensored_lpmf(data int d, int dist_id, array[] real params,\n                                data real pwindow, data real d_upper,\n                                data real D, int primary_id,\n                                array[] real primary_params) {\n  if (d_upper > D) {\n    reject(\"Upper truncation point is greater than D. It is \", d_upper,\n           \" and D is \", D, \". Resolve this by increasing D to be greater or equal to d + swindow or decreasing swindow.\");\n  }\n  if (d_upper <= d) {\n    reject(\"Upper truncation point is less than or equal to d. It is \", d_upper,\n           \" and d is \", d, \". Resolve this by increasing d to be less than d_upper.\");\n  }\n  real log_cdf_upper = primarycensored_lcdf(\n    d_upper | dist_id, params, pwindow, positive_infinity(), primary_id, primary_params\n  );\n  real log_cdf_lower = primarycensored_lcdf(\n    d | dist_id, params, pwindow, positive_infinity(), primary_id, primary_params\n  );\n  if (!is_inf(D)) {\n    real log_cdf_D;\n\n    if (d_upper == D) {\n      log_cdf_D = log_cdf_upper;\n    } else {\n      log_cdf_D = primarycensored_lcdf(\n        D | dist_id, params, pwindow, positive_infinity(), primary_id, primary_params\n      );\n    }\n    return log_diff_exp(log_cdf_upper, log_cdf_lower) - log_cdf_D;\n  } else {\n    return log_diff_exp(log_cdf_upper, log_cdf_lower);\n  }\n}"

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

    ## Stan functions written to: /tmp/RtmpOnP8lF/expgrowth_rng.stan

This can now be compiled and used in the same way as any other
`cmdstanr` model.

``` r
model <- cmdstan_model(expgrowth_rng_file)
model
```

    ## functions {
    ## // Stan functions from primarycensored version 1.3.0
    ## real expgrowth_rng(real min, real max, real r) {
    ##   real u = uniform_rng(0, 1);
    ##   if (abs(r) < 1e-10) {
    ##     return min + u * (max - min);
    ##   }
    ##   return min + log(u * (exp(r * max) - exp(r * min)) + exp(r * min)) / r;
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
