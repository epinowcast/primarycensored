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
#>  [5] "check_for_uniform_terms"                
#>  [6] "check_for_analytical"                   
#>  [7] "expgrowth_cdf"                          
#>  [8] "expgrowth_lcdf"                         
#>  [9] "primary_lcdf"                           
#> [10] "primary_lcdf_vec"                       
#> [11] "discretestep_lcdf"                      
#> [12] "hazards_to_pmf"                         
#> [13] "discretehazard_lcdf"                    
#> [14] "primarycensored_uniform_lcdf_from_terms"
#> [15] "primarycensored_gamma_uniform_terms"    
#> [16] "primarycensored_gamma_uniform_lcdf"     
#> [17] "lognormal_lcdf_underflows"              
#> [18] "primarycensored_lognormal_uniform_terms"
#> [19] "primarycensored_lognormal_uniform_lcdf" 
#> [20] "log_weibull_g"                          
#> [21] "primarycensored_weibull_uniform_terms"  
#> [22] "primarycensored_weibull_uniform_lcdf"   
#> [23] "gengamma_lcdf"                          
#> [24] "primarycensored_gengamma_uniform_terms" 
#> [25] "primarycensored_gengamma_uniform_lcdf"  
#> [26] "primarycensored_analytical_lcdf_raw"    
#> [27] "primarycensored_analytical_lcdf"        
#> [28] "primarycensored_analytical_cdf"         
#> [29] "primarycensored_cdf"                    
#> [30] "primarycensored_lcdf"                   
#> [31] "primarycensored_lpmf"                   

# A function with no dependencies
pcd_stan_function_deps("expgrowth_pdf")
#> [1] "expgrowth_pdf"
```
