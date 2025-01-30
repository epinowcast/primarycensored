
# primarycensored 1.0.0.1000

Development release.

## Package

- Updated the CI so that stan code is also tested on Windows and Mac. This is to ensure that the stan code is compatible with these platforms and in response to a CI bug in `epidist`.
- Revised approach to automatic discovery of distribution functions. This soft deprecates the `pdist_name` and `dprimary_name` arguments throughout. Users wishing to pass distribution names (i.e. to potentially leverage analytical solutions) are advised to use the newly introduced `add_name_attribute()` function. Adds transient dependency on `lifecycle` and `rlang` packages. See #188 by @pearsonca.

## Documentation

- Added a CRAN downloads badge to the README.
- Corrected how to specify an empty array in the docs of `primarycensored_lcdf()`.

## Bug fixes

- Added a missing `@family` tag to the `pcens` functions. This omission resulted in the Weibull analytical solution not being visible in the package documentation.
- Added precalculation of vector sizes to the `primarycensored_cdf()` stan function, avoiding errors on some platforms due to narrowing conversions in aggregate initialisation. 
- Changed `D` to be of type real in `pcens_model.stan` in order to support infinite `relative_obs_time`.

# primarycensored 1.0.0

This is the first major release of `primarycensored` and has been submitted to CRAN.

## Package

- Fix internal package misspelling of `primary_lpdf`.
- Move to "stable" lifecycle status.
- Added `rhub` checks to the `Github Actions` workflow.
- Added `dependencies: "hard"` to the `R-CMD-check` workflow to ensure checks pass without optional dependencies.
- Improved handling of examples that use optional dependencies.
- Check all URLs for redirects.
- Ensure that all functions have documented return values.

# primarycensored 0.6.0

This release renames the package to `primarycensored` from `primarycensoredist` and also renames many of the functions to remove the `dist` in their name. This was done to make the package name and the functions more consistent and to remove the need to use the `dist` suffix. It also aligns it with the new `PrimaryCensored.jl` package in our [Julia ecosystem](https://github.com/Epiaware).

Aside from name changes, this release also adds an analytical solution for the weibull distribution with uniform primary censoring, removes the need to assign functions to the global environment for `fitdistdoublecens()` by using `withr`, and adds a `check_truncation()` function to check if the truncation time is larger than the maximum observed delay. This is used in `fitdistdoublecens()` and `pcd_as_stan_data()` to ensure that the truncation time is appropriate to maximise computational efficiency.

## Package

* Removed the need to assign functions to the global environment for `fitdistdoublecens()` by using `withr`.
* Added a `check_truncation()` function to check if the truncation time is larger than the maximum observed delay. This is used in `fitdistdoublecens()` and `pcd_as_stan_data()` to ensure that the truncation time is appropriate to maximise computational efficiency.
* `pcd_as_cmdstan_data()` has been renamed to `pcd_as_stan_data()` to better reflect that it is used for `Stan` models in general rather than just the `CmdStan` models.
* The stan code has been refactored into a folder of functions within the current `stan` folder and the `stan` model has been moved into the `stan` folder. All paths to the stan code have been updated to reflect this.
* Added R and stan implementations of the primary censored cdf for the weibull distribution with uniform primary censoring.
* The package has been renamed to `primarycensored` as have all functions that use "dist" in their name.

## Documentation

* Simplified the "Analytic solutions" vignette by removing verbose derivation details.
* Added links between vignettes to make it easier to navigate the documentation.
* Added explicit usage of `pdist`, `dprimary`, `rdist`, and `rprimary` arguments in the getting started vignette to make it easier to link to mathematical details.
* Fixed error in "Analytic solutions" vignette where the Weibull density was not being treated as zero for negative delays.
* Split "Why it works" vignette into two separate vignettes, "Why it works" and "Analytic solutions for censored delay distributions".

# primarycensored 0.5.0

This release adds a new `{touchstone}` based benchmark suite to the package. It also adds a new "How it works" vignette which aims to give the reader more details into how the primary censored distributions work.

As part of the "How it works" we (@SamuelBrand1) found analytical solutions for the gamma, lognormal, and weibull distributions with uniform primary censoring. These are now implemented for the lognormal and gamma distributions in the `R` and `stan` code providing significant speedups to the fitting process (~10-20 times faster). The Weibull will be added in the next release.

## Package

* Add `{touchstone}` based benchmarks for benchmarking R utility functions, and fitting the `stan` and `fitdistplus` models.
* Added a "How it works" vignette.
* Added R infrastructure for analytical solutions via the `primarycensored` S3 class.
* Added Weibull analytical solution to "How it works" vignette.
* Added analytical solutions for the gamma and lognormal distributions with uniform primary censoring to both the `R` and `stan` code.
* Added numerical protection to ensure that CDFs for delays greater than the maximum truncation are exactly 1.

# primarycensored 0.4.0

In this release, we have added a new package `stan` model for fitting distributions using the `cmdstanr` package. We have also added a new function `fitdistdoublecens()` to allow for fitting of double censored and truncated data using the `fitdistrplus` package. As well as these functionality improvements this release focuses on improving the stability of the `stan` model and improving the speed of the `primarycensored_ode` function.

## Package

* Added a new function `fitdistdoublecens()` to allow for fitting of double censored and truncated data using the `fitdistrplus` package.
* Added low level tests for the Stan `primarycensored_ode` function.
* Rephrased the stan code to use a ODE solver rather than a numerical integration method. This allows for much faster and more stable computation of the likelihood
* Added a `CmdStan` model for fitting distributions using the `cmdstanr` package.
* Added helpers functions for working with the new `CmdStan` model and added an example to the vignette.
* Added parameter recovery tests for the new `CmdStan` model which tests the `primarycensored_lpmf` function when used with NUTS based fitting.

# primarycensored 0.3.0

This release fixes and improves truncation handling across the code base. It also adds a new vignette showcasing how to use the `primarycensored` and `fitdistrplus` packages together to fit distributions.

## Package

* Updated the approach to truncation to be outside the primary censored distribution integral.
* Improved tests that compare random sampling and probability mass/density functions between R and Stan.
* Improved cross-testing between R and Stan implementations of the primary censored distributions.
* Worked on improving the stability of the `primarycensored_lpmf` when used for NUTS based fitting (i.e. in Stan).

## Documentation

* @athowes improved the getting started vignette by catching a few grammar errors and simplifying language.
* Added a new vignette showcasing how to use the `primarycensored` and `fitdistrplus` packages together to fit distributions.

# primarycensored 0.2.0

This release puts in place initial documentation and vignettes. It also includes a new primary censored distribution interface to allow for non-secondary event censored distributions. Development of this release as identified some numerical issues in the
gradient evaluations for the primary censored distributions which may lead to breaking
interface changes in `0.3.0` for the Stan code.

## Package

* Added support for `swindow = 0` to `rprimarycensored` to allow for non-secondary event censored distributions.
* Adapted `rprimarycensored` so that truncation is based on the primary censored distribution before secondary events are censored. This better matches the generative process.
* Added a new Stan interface tool to enable finding which files functions are implemented in the Stan code.

## Documentation

* Added a getting started vignette.
* Added a vignette showcasing how to use the package Stan code with `cmdstanr`.
* Added a vignette showcasing how to fit distributions using the `cmdstanr` package.

# primarycensored 0.1.0

This is the initial `primarycensored` release and includes R and stan tools for dealing with potentially truncated primary event censored delay distributions. We expect all current features to work but the UI may change as the package matures over the next few versions.

## Package

* Added package skeleton.
* Added checking input functions.
* Added stan functions for primary censored and truncated distributions.
* Added R functions for primary censored and truncated distributions.
* Add R function to facilitate working with the Stan code.
* Added tests for primary censored and truncated distributions.
* Added tests to compare R and Stan implementations.
* Added tests for the R functions that facilitate working with the Stan code.
* Resolved R CMD check errors, warnings and notes.
* Added a hexsticker.
* Added vignette skeletons in preparation for `0.2.0` release.
