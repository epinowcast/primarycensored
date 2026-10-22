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
#>  [5] "primarycensored_numeric_cdf"                         
#>  [6] "check_for_uniform_terms"                             
#>  [7] "check_for_analytical"                                
#>  [8] "check_for_tilt_transform"                            
#>  [9] "check_for_exptilt"                                   
#> [10] "check_for_analytical_params"                         
#> [11] "expgrowth_cdf"                                       
#> [12] "expgrowth_lcdf"                                      
#> [13] "primary_lcdf"                                        
#> [14] "primary_lcdf_vec"                                    
#> [15] "discretestep_lcdf"                                   
#> [16] "hazards_to_pmf"                                      
#> [17] "discretehazard_lcdf"                                 
#> [18] "primarycensored_uniform_lcdf_from_terms"             
#> [19] "gamma_log_lead_logx"                                 
#> [20] "gamma_log_pq_logx"                                   
#> [21] "gamma_lcdf_logx"                                     
#> [22] "gamma_lcdf_logx_pair"                                
#> [23] "gamma_log_rate"                                      
#> [24] "primarycensored_gamma_uniform_terms"                 
#> [25] "primarycensored_gamma_uniform_lcdf"                  
#> [26] "lognormal_lcdf_underflows"                           
#> [27] "primarycensored_lognormal_uniform_terms"             
#> [28] "primarycensored_lognormal_uniform_lcdf"              
#> [29] "log_weibull_g"                                       
#> [30] "primarycensored_weibull_uniform_terms"               
#> [31] "primarycensored_weibull_uniform_lcdf"                
#> [32] "primarycensored_gengamma_uniform_terms"              
#> [33] "primarycensored_gengamma_uniform_lcdf"               
#> [34] "primarycensored_log_std_normal_cdf"                  
#> [35] "log_tilt_transform_pair"                             
#> [36] "primarycensored_log_diff_exp"                        
#> [37] "primarycensored_moments_from_partial"                
#> [38] "primarycensored_gamma_tilt_moments"                  
#> [39] "primarycensored_tilt_moments"                        
#> [40] "primarycensored_exptilt_form"                        
#> [41] "primarycensored_tail_diff"                           
#> [42] "primarycensored_exptilt_lcdf_from_terms"             
#> [43] "primarycensored_exptilt_small_window_lcdf_from_terms"
#> [44] "primarycensored_exptilt_small_delay_lcdf_from_terms" 
#> [45] "pstep_lcdf"                                          
#> [46] "phazard_lcdf"                                        
#> [47] "gengamma_lcdf"                                       
#> [48] "dist_lcdf"                                           
#> [49] "primarycensored_exptilt_lcdf"                        
#> [50] "primarycensored_analytical_lcdf_raw"                 
#> [51] "primarycensored_analytical_lcdf"                     
#> [52] "primarycensored_analytical_cdf"                      
#> [53] "primarycensored_cdf"                                 
#> [54] "primarycensored_lcdf"                                
#> [55] "primarycensored_lpmf"                                

# A function with no dependencies
pcd_stan_function_deps("expgrowth_pdf")
#> [1] "expgrowth_pdf"
```
