# Package index

## Primary event censored distribution functions

Functions for generating, evaluating density, and computing cumulative
probabilities of primary event censored distributions

- [`dprimarycensored()`](https://primarycensored.epinowcast.org/dev/reference/dprimarycensored.md)
  [`dpcens()`](https://primarycensored.epinowcast.org/dev/reference/dprimarycensored.md)
  : Compute the primary event censored PMF for delays
- [`pprimarycensored()`](https://primarycensored.epinowcast.org/dev/reference/pprimarycensored.md)
  [`ppcens()`](https://primarycensored.epinowcast.org/dev/reference/pprimarycensored.md)
  : Compute the primary event censored CDF for delays
- [`qprimarycensored()`](https://primarycensored.epinowcast.org/dev/reference/qprimarycensored.md)
  [`qpcens()`](https://primarycensored.epinowcast.org/dev/reference/qprimarycensored.md)
  : Compute quantiles corresponding to target probabilities for primary
  event censored delays
- [`rprimarycensored()`](https://primarycensored.epinowcast.org/dev/reference/rprimarycensored.md)
  [`rpcens()`](https://primarycensored.epinowcast.org/dev/reference/rprimarycensored.md)
  : Generate random samples from a primary event censored distribution

## Primary event distributions

Probability density and random generation functions for primary event
distributions

- [`dexpgrowth()`](https://primarycensored.epinowcast.org/dev/reference/expgrowth.md)
  [`pexpgrowth()`](https://primarycensored.epinowcast.org/dev/reference/expgrowth.md)
  [`rexpgrowth()`](https://primarycensored.epinowcast.org/dev/reference/expgrowth.md)
  : Exponential growth distribution functions

## Primary censored distribution class and methods

S3 class and methods for computing primary event censored distributions,
focusing on the internal machinery used by the package. Unlike the
primary event distributions section which deals with specific
distribution functions, this section covers the general framework for
handling censored distributions.

- [`new_pcens()`](https://primarycensored.epinowcast.org/dev/reference/new_pcens.md)
  : S3 class for primary event censored distribution computation
- [`pcens_cdf()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.md)
  : Compute primary event censored CDF
- [`pcens_cdf(`*`<default>`*`)`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.default.md)
  : Default method for computing primary event censored CDF
- [`pcens_cdf(`*`<pcens_pgamma_dunif>`*`)`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_pgamma_dunif.md)
  : Method for Gamma delay with uniform primary
- [`pcens_cdf(`*`<pcens_plnorm_dunif>`*`)`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_plnorm_dunif.md)
  : Method for Log-Normal delay with uniform primary
- [`pcens_cdf(`*`<pcens_pweibull_dunif>`*`)`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_pweibull_dunif.md)
  : Method for Weibull delay with uniform primary
- [`pcens_quantile()`](https://primarycensored.epinowcast.org/dev/reference/pcens_quantile.md)
  : Compute primary event censored quantiles
- [`pcens_quantile(`*`<default>`*`)`](https://primarycensored.epinowcast.org/dev/reference/pcens_quantile.default.md)
  : Default method for computing primary event censored quantiles

## Utility functions

Utility functions for working withe package

- [`add_name_attribute()`](https://primarycensored.epinowcast.org/dev/reference/add_name_attribute.md)
  : Helper method for custom distributions
- [`pcd_dist_name()`](https://primarycensored.epinowcast.org/dev/reference/pcd_dist_name.md)
  : Get distribution function cdf or pdf name
- [`pcd_distributions`](https://primarycensored.epinowcast.org/dev/reference/pcd_distributions.md)
  : Supported delay distributions
- [`pcd_primary_distributions`](https://primarycensored.epinowcast.org/dev/reference/pcd_primary_distributions.md)
  : Supported primary event distributions

## Distribution checking functions

Functions to validate cumulative distribution functions (CDFs) and
probability density functions (PDFs)

- [`check_dprimary()`](https://primarycensored.epinowcast.org/dev/reference/check_dprimary.md)
  : Check if a function is a valid bounded probability density function
  (PDF)
- [`check_pdist()`](https://primarycensored.epinowcast.org/dev/reference/check_pdist.md)
  : Check if a function is a valid cumulative distribution function
  (CDF)
- [`check_truncation()`](https://primarycensored.epinowcast.org/dev/reference/check_truncation.md)
  : Check if truncation time is appropriate relative to the maximum
  delay

## Tools for working with package Stan functions

Utility functions for interfacing with Stan models and extracting
results

- [`pcd_load_stan_functions()`](https://primarycensored.epinowcast.org/dev/reference/pcd_load_stan_functions.md)
  : Load Stan functions as a string
- [`pcd_stan_dist_id()`](https://primarycensored.epinowcast.org/dev/reference/pcd_stan_dist_id.md)
  : Get distribution stan ID by name
- [`pcd_stan_files()`](https://primarycensored.epinowcast.org/dev/reference/pcd_stan_files.md)
  : Get Stan files containing specified functions
- [`pcd_stan_functions()`](https://primarycensored.epinowcast.org/dev/reference/pcd_stan_functions.md)
  : Get Stan function names from Stan files
- [`pcd_stan_path()`](https://primarycensored.epinowcast.org/dev/reference/pcd_stan_path.md)
  : Get the path to the Stan code

## Wrappers facilitating the use of other modelling packages

Functions that wrap around external packages like fitdistrplus to fit
distributions to doubly censored data

- [`fitdistdoublecens()`](https://primarycensored.epinowcast.org/dev/reference/fitdistdoublecens.md)
  : Fit a distribution to doubly censored data
- [`pcd_as_stan_data()`](https://primarycensored.epinowcast.org/dev/reference/pcd_as_stan_data.md)
  : Prepare data for primarycensored Stan model
- [`pcd_cmdstan_model()`](https://primarycensored.epinowcast.org/dev/reference/pcd_cmdstan_model.md)
  : Create a CmdStanModel with primarycensored Stan functions
