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
#>   [7] "gamma_log_lead_logx"                                 
#>   [8] "gamma_log_pq_logx"                                   
#>   [9] "gamma_lcdf_logx"                                     
#>  [10] "gamma_log_rate"                                      
#>  [11] "lognormal_lcdf_underflows"                           
#>  [12] "pstep_lcdf"                                          
#>  [13] "hazards_to_pmf"                                      
#>  [14] "phazard_lcdf"                                        
#>  [15] "primarycensored_loglogistic_lcdf"                    
#>  [16] "gengamma_lcdf"                                       
#>  [17] "dist_lcdf"                                           
#>  [18] "gumbel_delay_lcdf"                                   
#>  [19] "gumbel_numeric_delay"                                
#>  [20] "tgumbel_log_delta_window"                            
#>  [21] "tgumbel_log1m_exp_neg_exp"                           
#>  [22] "tgumbel_log_norm"                                    
#>  [23] "tgumbel_lpdf"                                        
#>  [24] "gumbel_numeric_point"                                
#>  [25] "primarycensored_gumbel_ode"                          
#>  [26] "gumbel_numeric_piece_lcdf"                           
#>  [27] "gumbel_numeric_spike_lcdf"                           
#>  [28] "gumbel_numeric_z_limits"                             
#>  [29] "gumbel_z_log_integrand"                              
#>  [30] "gumbel_numeric_z_half"                               
#>  [31] "gumbel_numeric_z_lcdf"                               
#>  [32] "primarycensored_gumbel_numeric_lcdf"                 
#>  [33] "primarycensored_gumbel_numeric_cdf"                  
#>  [34] "expgrowth_lpdf"                                      
#>  [35] "primarycensored_log_diff_exp"                        
#>  [36] "tlogis_log_mass"                                     
#>  [37] "tlogis_lpdf"                                         
#>  [38] "primary_lpdf"                                        
#>  [39] "primarycensored_ode"                                 
#>  [40] "tlogis_spike_times"                                  
#>  [41] "primarycensored_numeric_cdf"                         
#>  [42] "check_for_uniform_terms"                             
#>  [43] "check_for_analytical"                                
#>  [44] "check_for_tilt_transform"                            
#>  [45] "check_for_exptilt"                                   
#>  [46] "check_for_gumbel"                                    
#>  [47] "gumbel_n_terms"                                      
#>  [48] "check_for_gumbel_params"                             
#>  [49] "check_for_tlogis"                                    
#>  [50] "tlogis_log_bound"                                    
#>  [51] "tlogis_series_terms"                                 
#>  [52] "tlogis_plan"                                         
#>  [53] "tlogis_is_valid"                                     
#>  [54] "tlogis_n_pos"                                        
#>  [55] "tlogis_n_neg"                                        
#>  [56] "tlogis_plan_applies"                                 
#>  [57] "check_for_tlogis_params"                             
#>  [58] "check_for_analytical_params"                         
#>  [59] "check_for_tilt_transform_at"                         
#>  [60] "check_for_analytical_delay"                          
#>  [61] "expgrowth_cdf"                                       
#>  [62] "expgrowth_lcdf"                                      
#>  [63] "tgumbel_lcdf"                                        
#>  [64] "tlogis_lcdf"                                         
#>  [65] "primary_lcdf"                                        
#>  [66] "primary_lcdf_vec"                                    
#>  [67] "discretestep_lcdf"                                   
#>  [68] "discretehazard_lcdf"                                 
#>  [69] "primarycensored_uniform_lcdf_from_terms"             
#>  [70] "gamma_lcdf_logx_pair"                                
#>  [71] "primarycensored_gamma_uniform_terms"                 
#>  [72] "primarycensored_gamma_uniform_lcdf"                  
#>  [73] "primarycensored_lognormal_uniform_terms"             
#>  [74] "primarycensored_lognormal_uniform_lcdf"              
#>  [75] "log_weibull_g"                                       
#>  [76] "primarycensored_weibull_uniform_terms"               
#>  [77] "primarycensored_weibull_uniform_lcdf"                
#>  [78] "primarycensored_gengamma_uniform_terms"              
#>  [79] "primarycensored_gengamma_uniform_lcdf"               
#>  [80] "loglogistic_log_error_scale"                         
#>  [81] "loglogistic_window_ill_conditioned"                  
#>  [82] "primarycensored_loglogistic_lccdf"                   
#>  [83] "loglogistic_numeric_ode"                             
#>  [84] "loglogistic_numeric_lcdf"                            
#>  [85] "loglogistic_ratio_small"                             
#>  [86] "loglogistic_partial_at_one"                          
#>  [87] "loglogistic_tau_integral"                            
#>  [88] "loglogistic_tail_term"                               
#>  [89] "loglogistic_ratio_large"                             
#>  [90] "loglogistic_moment_ratio"                            
#>  [91] "primarycensored_loglogistic_uniform_terms"           
#>  [92] "primarycensored_loglogistic_uniform_lcdf"            
#>  [93] "primarycensored_gauss_legendre_nodes"                
#>  [94] "primarycensored_gauss_legendre_log_weights"          
#>  [95] "primarycensored_lognormal_tilt_one_panel"            
#>  [96] "primarycensored_lognormal_tilt_panel"                
#>  [97] "primarycensored_lognormal_tilt_bump"                 
#>  [98] "log_tilt_transform_context"                          
#>  [99] "primarycensored_exptilt_form"                        
#> [100] "primarycensored_moments_from_partial"                
#> [101] "primarycensored_gamma_tilt_moments"                  
#> [102] "primarycensored_stacy_tilt_moments"                  
#> [103] "primarycensored_tilt_moments"                        
#> [104] "primarycensored_stacy_series"                        
#> [105] "primarycensored_lognormal_tilt_quadrature"           
#> [106] "primarycensored_lognormal_tilt_series"               
#> [107] "primarycensored_lognormal_tilt_pair"                 
#> [108] "log_tilt_transform_pair_shared"                      
#> [109] "log_tilt_transform_pair"                             
#> [110] "primarycensored_exptilt_terms"                       
#> [111] "primarycensored_tail_diff"                           
#> [112] "primarycensored_exptilt_lcdf_from_terms"             
#> [113] "primarycensored_exptilt_small_window_lcdf_from_terms"
#> [114] "primarycensored_exptilt_small_delay_lcdf_from_terms" 
#> [115] "primarycensored_exptilt_lcdf_shared"                 
#> [116] "primarycensored_exptilt_lcdf"                        
#> [117] "gumbel_series_accepted"                              
#> [118] "primarycensored_gumbel_terms"                        
#> [119] "gumbel_finite_abs"                                   
#> [120] "primarycensored_gumbel_lcdf_from_terms"              
#> [121] "primarycensored_gumbel_lcdf"                         
#> [122] "tlogis_weights"                                      
#> [123] "primarycensored_tlogis_terms"                        
#> [124] "tlogis_terms_diff"                                   
#> [125] "tlogis_weighted_sum"                                 
#> [126] "tlogis_part_before_window"                           
#> [127] "tlogis_part_above_location"                          
#> [128] "tlogis_part_below_location"                          
#> [129] "primarycensored_tlogis_lcdf_from_terms"              
#> [130] "tlogis_is_small_delay"                               
#> [131] "primarycensored_tlogis_small_delay_lcdf"             
#> [132] "primarycensored_tlogis_lcdf"                         
#> [133] "primarycensored_analytical_lcdf_raw"                 
#> [134] "primarycensored_analytical_lcdf"                     
#> [135] "primarycensored_analytical_cdf"                      
#> [136] "primarycensored_cdf"                                 
#> [137] "primarycensored_lcdf"                                
#> [138] "primarycensored_exptilt_lcdf_pair"                   
#> [139] "primarycensored_lpmf"                                

# A function with no dependencies
pcd_stan_function_deps("expgrowth_pdf")
#> [1] "expgrowth_pdf"
```
