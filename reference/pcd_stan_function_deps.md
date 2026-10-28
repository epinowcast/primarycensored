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
#>   [1] "primarycensored_log_normalizer"                      
#>   [2] "primarycensored_apply_truncation"                    
#>   [3] "dist_has_positive_support"                           
#>   [4] "primarycensored_truncation_bounds"                   
#>   [5] "gumbel_log_sum_exp"                                  
#>   [6] "primarycensored_log_std_normal_cdf"                  
#>   [7] "gamma_lcdf_logx"                                     
#>   [8] "gamma_lcdf_lccdf_logx"                               
#>   [9] "gamma_log_rate"                                      
#>  [10] "lognormal_lcdf_underflows"                           
#>  [11] "pstep_lcdf"                                          
#>  [12] "hazards_to_pmf"                                      
#>  [13] "phazard_lcdf"                                        
#>  [14] "primarycensored_loglogistic_lcdf"                    
#>  [15] "gengamma_lcdf"                                       
#>  [16] "dist_lcdf"                                           
#>  [17] "gumbel_delay_lcdf"                                   
#>  [18] "gumbel_numeric_delay"                                
#>  [19] "tgumbel_log_delta_window"                            
#>  [20] "tgumbel_log1m_exp_neg_exp"                           
#>  [21] "tgumbel_log_norm"                                    
#>  [22] "tgumbel_lpdf"                                        
#>  [23] "gumbel_numeric_point"                                
#>  [24] "primarycensored_gumbel_ode"                          
#>  [25] "gumbel_numeric_piece_lcdf"                           
#>  [26] "gumbel_numeric_spike_lcdf"                           
#>  [27] "gumbel_numeric_z_limits"                             
#>  [28] "gumbel_z_log_integrand"                              
#>  [29] "gumbel_numeric_z_half"                               
#>  [30] "gumbel_numeric_z_lcdf"                               
#>  [31] "primarycensored_gumbel_numeric_lcdf"                 
#>  [32] "primarycensored_gumbel_numeric_cdf"                  
#>  [33] "expgrowth_lpdf"                                      
#>  [34] "primarycensored_log_diff_exp"                        
#>  [35] "tlogis_log_mass"                                     
#>  [36] "tlogis_lpdf"                                         
#>  [37] "primary_lpdf"                                        
#>  [38] "primarycensored_ode"                                 
#>  [39] "tlogis_spike_times"                                  
#>  [40] "primarycensored_numeric_cdf"                         
#>  [41] "check_for_uniform_terms"                             
#>  [42] "check_for_analytical"                                
#>  [43] "check_for_tilt_transform"                            
#>  [44] "check_for_exptilt"                                   
#>  [45] "check_for_gumbel"                                    
#>  [46] "gumbel_n_terms"                                      
#>  [47] "check_for_gumbel_params"                             
#>  [48] "check_for_tlogis"                                    
#>  [49] "tlogis_log_bound"                                    
#>  [50] "tlogis_series_terms"                                 
#>  [51] "tlogis_plan"                                         
#>  [52] "tlogis_is_valid"                                     
#>  [53] "tlogis_n_pos"                                        
#>  [54] "tlogis_n_neg"                                        
#>  [55] "tlogis_plan_applies"                                 
#>  [56] "check_for_tlogis_params"                             
#>  [57] "check_for_analytical_params"                         
#>  [58] "check_for_tilt_transform_at"                         
#>  [59] "check_for_analytical_delay"                          
#>  [60] "expgrowth_cdf"                                       
#>  [61] "expgrowth_lcdf"                                      
#>  [62] "tgumbel_lcdf"                                        
#>  [63] "tlogis_lcdf"                                         
#>  [64] "primary_lcdf"                                        
#>  [65] "primary_lcdf_vec"                                    
#>  [66] "discretestep_lcdf"                                   
#>  [67] "discretehazard_lcdf"                                 
#>  [68] "primarycensored_uniform_lcdf_from_terms"             
#>  [69] "primarycensored_gamma_uniform_terms"                 
#>  [70] "primarycensored_gamma_uniform_lcdf"                  
#>  [71] "primarycensored_lognormal_uniform_terms"             
#>  [72] "primarycensored_lognormal_uniform_lcdf"              
#>  [73] "log_weibull_g"                                       
#>  [74] "primarycensored_weibull_uniform_terms"               
#>  [75] "primarycensored_weibull_uniform_lcdf"                
#>  [76] "primarycensored_gengamma_uniform_terms"              
#>  [77] "primarycensored_gengamma_uniform_lcdf"               
#>  [78] "loglogistic_log_error_scale"                         
#>  [79] "loglogistic_window_ill_conditioned"                  
#>  [80] "primarycensored_loglogistic_lccdf"                   
#>  [81] "loglogistic_numeric_ode"                             
#>  [82] "loglogistic_numeric_lcdf"                            
#>  [83] "loglogistic_ratio_small"                             
#>  [84] "loglogistic_partial_at_one"                          
#>  [85] "loglogistic_tau_integral"                            
#>  [86] "loglogistic_tail_term"                               
#>  [87] "loglogistic_ratio_large"                             
#>  [88] "loglogistic_moment_ratio"                            
#>  [89] "primarycensored_loglogistic_uniform_terms"           
#>  [90] "primarycensored_loglogistic_uniform_lcdf"            
#>  [91] "primarycensored_gauss_legendre_nodes"                
#>  [92] "primarycensored_gauss_legendre_log_weights"          
#>  [93] "primarycensored_lognormal_tilt_one_panel"            
#>  [94] "primarycensored_lognormal_tilt_panel"                
#>  [95] "primarycensored_lognormal_tilt_bump"                 
#>  [96] "log_tilt_transform_context"                          
#>  [97] "primarycensored_exptilt_form"                        
#>  [98] "primarycensored_moments_from_partial"                
#>  [99] "primarycensored_gamma_tilt_moments"                  
#> [100] "primarycensored_stacy_tilt_moments"                  
#> [101] "primarycensored_tilt_moments"                        
#> [102] "primarycensored_stacy_series"                        
#> [103] "primarycensored_lognormal_tilt_quadrature"           
#> [104] "primarycensored_lognormal_tilt_series"               
#> [105] "primarycensored_lognormal_tilt_pair"                 
#> [106] "log_tilt_transform_pair_shared"                      
#> [107] "log_tilt_transform_pair"                             
#> [108] "primarycensored_exptilt_terms"                       
#> [109] "primarycensored_tail_diff"                           
#> [110] "primarycensored_exptilt_lcdf_from_terms"             
#> [111] "primarycensored_exptilt_small_window_lcdf_from_terms"
#> [112] "primarycensored_exptilt_small_delay_lcdf_from_terms" 
#> [113] "primarycensored_exptilt_lcdf_shared"                 
#> [114] "primarycensored_exptilt_lcdf"                        
#> [115] "gumbel_series_accepted"                              
#> [116] "primarycensored_gumbel_terms"                        
#> [117] "gumbel_finite_abs"                                   
#> [118] "primarycensored_gumbel_lcdf_from_terms"              
#> [119] "primarycensored_gumbel_lcdf"                         
#> [120] "tlogis_weights"                                      
#> [121] "primarycensored_tlogis_terms"                        
#> [122] "tlogis_terms_diff"                                   
#> [123] "tlogis_weighted_sum"                                 
#> [124] "tlogis_part_before_window"                           
#> [125] "tlogis_part_above_location"                          
#> [126] "tlogis_part_below_location"                          
#> [127] "primarycensored_tlogis_lcdf_from_terms"              
#> [128] "tlogis_is_small_delay"                               
#> [129] "primarycensored_tlogis_small_delay_lcdf"             
#> [130] "primarycensored_tlogis_lcdf"                         
#> [131] "primarycensored_analytical_lcdf_raw"                 
#> [132] "primarycensored_analytical_lcdf"                     
#> [133] "primarycensored_analytical_cdf"                      
#> [134] "primarycensored_cdf"                                 
#> [135] "primarycensored_lcdf"                                
#> [136] "primarycensored_exptilt_lcdf_pair"                   
#> [137] "primarycensored_lpmf"                                

# A function with no dependencies
pcd_stan_function_deps("expgrowth_pdf")
#> [1] "expgrowth_pdf"
```
