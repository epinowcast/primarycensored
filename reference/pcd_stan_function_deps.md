# Get dependencies for a Stan function

Returns all Stan functions that the specified function depends on, in
topological order (dependencies before the functions that use them).

## Usage

``` r
pcd_stan_function_deps(
  function_name,
  stan_path = primarycensored::pcd_stan_path()
)
```

## Arguments

- function_name:

  Character string, the name of the Stan function.

- stan_path:

  Character string specifying the path to the directory containing Stan
  files. Defaults to the Stan path of the primarycensored package.

## Value

A character vector of function names that the specified function depends
on, ordered so that dependencies come before functions that use them.
The requested function itself is included as the last element.

## See also

Tools for working with package Stan functions
[`pcd_load_stan_functions()`](https://primarycensored.epinowcast.org/reference/pcd_load_stan_functions.md),
[`pcd_stan_dist_id()`](https://primarycensored.epinowcast.org/reference/pcd_stan_dist_id.md),
[`pcd_stan_files()`](https://primarycensored.epinowcast.org/reference/pcd_stan_files.md),
[`pcd_stan_functions()`](https://primarycensored.epinowcast.org/reference/pcd_stan_functions.md),
[`pcd_stan_path()`](https://primarycensored.epinowcast.org/reference/pcd_stan_path.md)

## Examples

``` r
# See what primarycensored_lpmf depends on
pcd_stan_function_deps("primarycensored_lpmf")
#>  [1] "primarycensored_log_normalizer"        
#>  [2] "primarycensored_apply_truncation"      
#>  [3] "dist_has_positive_support"             
#>  [4] "primarycensored_truncation_bounds"     
#>  [5] "check_for_analytical"                  
#>  [6] "primarycensored_gamma_uniform_lcdf"    
#>  [7] "primarycensored_lognormal_uniform_lcdf"
#>  [8] "log_weibull_g"                         
#>  [9] "primarycensored_weibull_uniform_lcdf"  
#> [10] "primarycensored_analytical_lcdf_raw"   
#> [11] "primarycensored_analytical_lcdf"       
#> [12] "primarycensored_analytical_cdf"        
#> [13] "primarycensored_cdf"                   
#> [14] "primarycensored_lcdf"                  
#> [15] "primarycensored_lpmf"                  

# A function with no dependencies
pcd_stan_function_deps("expgrowth_pdf")
#> [1] "expgrowth_pdf"
```
