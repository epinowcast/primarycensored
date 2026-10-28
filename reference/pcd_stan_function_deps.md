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
#>  [1] "primarycensored_log_cdf_diff"             
#>  [2] "primarycensored_log_normalizer"           
#>  [3] "primarycensored_apply_truncation"         
#>  [4] "dist_has_positive_support"                
#>  [5] "primarycensored_truncation_bounds"        
#>  [6] "check_uniform_terms_params"               
#>  [7] "check_for_uniform_terms"                  
#>  [8] "check_for_analytical"                     
#>  [9] "expgrowth_cdf"                            
#> [10] "expgrowth_lcdf"                           
#> [11] "primary_lcdf"                             
#> [12] "primary_lcdf_vec"                         
#> [13] "discretestep_lcdf"                        
#> [14] "hazards_to_pmf"                           
#> [15] "discretehazard_lcdf"                      
#> [16] "primarycensored_uniform_lcdf_from_terms"  
#> [17] "gamma_lcdf_logx"                          
#> [18] "gamma_log_rate"                           
#> [19] "primarycensored_gamma_uniform_terms"      
#> [20] "lognormal_lcdf_underflows"                
#> [21] "primarycensored_lognormal_uniform_terms"  
#> [22] "log_weibull_g"                            
#> [23] "primarycensored_weibull_uniform_terms"    
#> [24] "gengamma_lcdf"                            
#> [25] "primarycensored_gengamma_uniform_terms"   
#> [26] "primarycensored_exponential_uniform_terms"
#> [27] "primarycensored_normal_uniform_terms"     
#> [28] "primarycensored_beta_uniform_terms"       
#> [29] "primarycensored_inv_gamma_lcdf"           
#> [30] "primarycensored_invgamma_uniform_terms"   
#> [31] "primarycensored_pareto_uniform_terms"     
#> [32] "primarycensored_uniform_terms"            
#> [33] "primarycensored_analytical_lcdf_raw"      
#> [34] "primarycensored_analytical_lcdf"          
#> [35] "primarycensored_analytical_cdf"           
#> [36] "primarycensored_cdf"                      
#> [37] "primarycensored_lcdf"                     
#> [38] "primarycensored_lpmf"                     

# A function with no dependencies
pcd_stan_function_deps("expgrowth_pdf")
#> [1] "expgrowth_pdf"
```
