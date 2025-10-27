# Supported delay distributions

A dataset containing information about the supported delay distributions
in primarycensored. Includes both distributions with base R
implementations and those only available in Stan. Distributions beyond
these are not supported in the stan code but any user functions can be
used in the R code.

## Usage

``` r
pcd_distributions
```

## Format

A data.frame with 17 rows and 4 columns:

- name:

  Distribution name

- pdist:

  R distribution function name (e.g. plnorm), NA if there is no base R
  implementation

- aliases:

  Alternative names/identifiers

- stan_id:

  Stan distribution ID used in the stan code

## See also

Utility functions for working with distributions
[`add_name_attribute()`](https://primarycensored.epinowcast.org/reference/add_name_attribute.md),
[`pcd_dist_name()`](https://primarycensored.epinowcast.org/reference/pcd_dist_name.md),
[`pcd_primary_distributions`](https://primarycensored.epinowcast.org/reference/pcd_primary_distributions.md)
