# primarycensoreddist 0.4.0.1000

This is the current development version.

## Package

* Add `{touchstone}` based benchmarks for benchmarking R utility functions, and fitting the `stan` and `fitdistplus` models.

# primarycensoreddist 0.4.0

In this release, we have added a new package `stan` model for fitting distributions using the `cmdstanr` package. We have also added a new function `fitdistdoublecens()` to allow for fitting of double censored and truncated data using the `fitdistrplus` package. As well as these functionality improvements this release focuses on improving the stability of the `stan` model and improving the speed of the `primary_censored_ode` function.

## Package

* Added a new function `fitdistdoublecens()` to allow for fitting of double censored and truncated data using the `fitdistrplus` package.
* Added low level tests for the Stan `primary_censored_ode` function.
* Rephrased the stan code to use a ODE solver rather than a numerical integration method. This allows for much faster and more stable computation of the likelihood
* Added a `CmdStan` model for fitting distributions using the `cmdstanr` package.
* Added helpers functions for working with the new `CmdStan` model and added an example to the vignette.
* Added parameter recovery tests for the new `CmdStan` model which tests the `primary_censored_dist_lpmf` function when used with NUTS based fitting.

# primarycensoreddist 0.3.0

This release fixes and improves truncation handling across the code base. It also adds a new vignette showcasing how to use the `primarycensoreddist` and `fitdistrplus` packages together to fit distributions.

## Package

* Updated the approach to truncation to be outside the primary censored distribution integral.
* Improved tests that compare random sampling and probability mass/density functions between R and Stan.
* Improved cross-testing between R and Stan implementations of the primary censored distributions.
* Worked on improving the stability of the `primary_censored_dist_lpmf` when used for NUTS based fitting (i.e. in Stan).

## Documentation

* @athowes improved the getting started vignette by catching a few grammar errors and simplifying language.
* Added a new vignette showcasing how to use the `primarycensoreddist` and `fitdistrplus` packages together to fit distributions.

# primarycensoreddist 0.2.0

This release puts in place initial documentation and vignettes. It also includes a new primary censored distribution interface to allow for non-secondary event censored distributions. Development of this release as identified some numerical issues in the
gradient evaluations for the primary censored distributions which may lead to breaking
interface changes in `0.3.0` for the Stan code.

## Package

* Added support for `swindow = 0` to `rprimarycensoreddist` to allow for non-secondary event censored distributions.
* Adapted `rprimarycensoreddist` so that truncation is based on the primary censored distribution before secondary events are censored. This better matches the generative process.
* Added a new Stan interface tool to enable finding which files functions are implemented in the Stan code.

## Documentation

* Added a getting started vignette.
* Added a vignette showcasing how to use the package Stan code with `cmdstanr`.
* Added a vignette showcasing how to fit distributions using the `cmdstanr` package.

# primarycensoreddist 0.1.0

This is the initial `primarycensoreddist` release and includes R and stan tools for dealing with potentially truncated primary event censored delay distributions. We expect all current features to work but the UI may change as the package matures over the next few versions.

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
