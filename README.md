
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Primary Event Censored Distributions in R and Stan <a href="https://primarycensoreddist.epinowcast.org/"><img src="man/figures/logo.png" align="right" height="139" alt="primarycensoreddist website" /></a>

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![R-CMD-check](https://github.com/epinowcast/primarycensoreddist/workflows/R-CMD-check/badge.svg)](https://github.com/epinowcast/primarycensoreddist/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/epinowcast/primarycensoreddist/branch/main/graph/badge.svg)](https://app.codecov.io/gh/epinowcast/primarycensoreddist)
[![Universe](https://epinowcast.r-universe.dev/badges/primarycensoreddist)](https://epinowcast.r-universe.dev/primarycensoreddist)
[![MIT
license](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/epinowcast/primarycensoreddist/blob/master/LICENSE.md/)
[![GitHub
contributors](https://img.shields.io/github/contributors/epinowcast/primarycensoreddist)](https://github.com/epinowcast/primarycensoreddist/graphs/contributors)
<!-- badges: end -->

## Summary

This package provides both R functions for working with primary event
censored distributions and Stan implementations for use in Bayesian
modeling. Primary event censored distributions are useful for modeling
delayed reporting scenarios in epidemiology and other fields. It
provides support for arbitrary delay distributions, a range of common
primary distributions, and allows for truncation and secondary event
censoring to be accounted for.

## Installation

<details>
<summary>
Installing the package
</summary>

You can install the latest released version using the normal `R`
function, though you need to point to `r-universe` instead of CRAN:

``` r
install.packages(
  "primarycensoreddist",
  repos = "https://epinowcast.r-universe.dev"
)
```

Alternatively, you can use the [`remotes`
package](https://remotes.r-lib.org/) to install the development version
from Github (warning! this version may contain breaking changes and/or
bugs):

``` r
remotes::install_github(
  "epinowcast/primarycensoreddist",
  dependencies = TRUE
)
```

Similarly, you can install historical versions by specifying the release
tag (e.g. this installs
[`0.2.0`](https://github.com/epinowcast/primarycensoreddist/releases/tag/v0.2.0)):

``` r
remotes::install_github(
  "epinowcast/primarycensoreddist",
  dependencies = TRUE, ref = "v0.2.0"
)
```

*Note: You can also use that last approach to install a specific commit
if needed, e.g. if you want to try out a specific unreleased feature,
but not the absolute latest developmental version.*

</details>
<details>
<summary>
Installing CmdStan (optional for Stan functionality)
</summary>

If you wish to use the Stan functions, you will need to install
[CmdStan](https://mc-stan.org/users/interfaces/cmdstan), which also
entails having a suitable C++ toolchain setup. We recommend using the
[`cmdstanr` package](https://mc-stan.org/cmdstanr/). The Stan team
provides instructions in the [*Getting started with
`cmdstanr`*](https://mc-stan.org/cmdstanr/articles/cmdstanr.html)
vignette, with other details and support at the [package
site](https://mc-stan.org/cmdstanr/) along with some key instructions
available in the [Stan resources package
vignette](https://package.epinowcast.org/articles/stan-help.html#toolchain),
but the brief version is:

``` r
# if you not yet installed `primarycensoreddist`, or you installed it without
# `Suggests` dependencies
install.packages(
  "cmdstanr",
  repos = c("https://stan-dev.r-universe.dev", getOption("repos"))
)
# once `cmdstanr` is installed:
cmdstanr::install_cmdstan()
```

*Note: You can speed up CmdStan installation using the `cores` argument.
If you are installing a particular version of `epinowcast`, you may also
need to install a past version of CmdStan, which you can do with the
`version` argument.*

</details>

## Resources

We provide a range of other documentation, case studies, and community
spaces to ask (and answer!) questions:

<details>
<summary>
Package Website
</summary>

The [`primarycensoreddist`
website](https://primarycensoreddist.epinowcast.org/) includes a
function reference, model outline, and case studies using the package.
The site mainly concerns the release version, but you can also find
documentation for [the latest development
version](https://primarycensoreddist.epinowcast.org/dev/).

</details>
<details>
<summary>
Vignettes
</summary>

We have created [package
vignettes](https://primarycensoreddist.epinowcast.org/articles) to help
you get started with primarycensoreddist and to highlight other features
with case studies.

</details>
<details>
<summary>
Organisation Website
</summary>

Our [organisation website](https://www.epinowcast.org/) includes links
to other resources, [guest posts](https://www.epinowcast.org/blog.html),
and [seminar schedule](https://www.epinowcast.org/seminars.html) for
both upcoming and past recordings.

</details>
<details>
<summary>
Community Forum
</summary>

Our [community forum](https://community.epinowcast.org/) has areas for
[question and answer](https://community.epinowcast.org/c/interface/15)
and [considering new methods and
tools](https://community.epinowcast.org/c/projects/11), among others. If
you are generally interested in real-time analysis of infectious
disease, you may find this useful even if you do not use
`primarycensoreddist`.

</details>

## Contributing

We welcome contributions and new contributors! We particularly
appreciate help on [identifying and identified
issues](https://github.com/epinowcast/primarycensoreddist/issues).
Please check and add to the issues, and/or add a [pull
request](https://github.com/epinowcast/primarycensoreddist/pulls) and
see our [contributing
guide](https://github.com/epinowcast/.github/blob/main/CONTRIBUTING.md)
for more information.

If you need a different underlying model for your work:
`primarycensoreddist` provides a flexible framework for censored
distributions in both R and Stan. If you implement new distributions or
censoring mechanisms that expand the overall flexibility or improve the
defaults, please let us know either here or on the [community
forum](https://community.epinowcast.org/). We always like to hear about
new use-cases and extensions to the package.

### How to make a bug report or feature request

Please briefly describe your problem and what output you expect in an
[issue](https://github.com/epinowcast/primarycensoreddist/issues). If
you have a question, please don’t open an issue. Instead, ask on our [Q
and A
page](https://github.com/epinowcast/primarycensoreddist/discussions/categories/q-a).
See our [contributing
guide](https://github.com/epinowcast/.github/blob/main/CONTRIBUTING.md)
for more information.

### Code of Conduct

Please note that the `primarycensoreddist` project is released with a
[Contributor Code of
Conduct](https://github.com/epinowcast/.github/blob/main/CODE_OF_CONDUCT.md).
By contributing to this project, you agree to abide by its terms.

## Citation

If making use of our methodology or the methodology on which ours is
based, please cite the relevant papers from our [methods
outline](https://primarycensoreddist.epinowcast.org/articles/methods.html).
If you use `primarycensoreddist` in your work, please consider citing it
with `citation("primarycensoreddist")`.

## Contributors

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->

All contributions to this project are gratefully acknowledged using the
[`allcontributors`
package](https://github.com/ropenscilabs/allcontributors) following the
[all-contributors](https://allcontributors.org) specification.
Contributions of any kind are welcome!

### Code

<a href="https://github.com/epinowcast/primarycensoreddist/commits?author=seabbs">seabbs</a>

### Issues

<a href="https://github.com/epinowcast/primarycensoreddist/issues?q=is%3Aissue+author%3Azsusswein">zsusswein</a>,
<a href="https://github.com/epinowcast/primarycensoreddist/issues?q=is%3Aissue+author%3ASamuelBrand1">SamuelBrand1</a>

<!-- markdownlint-enable -->
<!-- prettier-ignore-end -->
<!-- ALL-CONTRIBUTORS-LIST:END -->
