# Supported primary event distributions

A dataset containing information about the supported primary event
distributions in primarycensored. Distributions beyond these are not
supported in the stan code but any user functions can be used in the R
code.

## Usage

``` r
pcd_primary_distributions
```

## Format

A data.frame with 2 rows and 4 columns:

- name:

  Distribution name

- dprimary:

  R density function name

- aliases:

  Alternative names/identifiers

- stan_id:

  Stan distribution ID used in the stan code

## See also

Utility functions for working with distributions
[`add_name_attribute()`](https://primarycensored.epinowcast.org/reference/add_name_attribute.md),
[`pcd_dist_name()`](https://primarycensored.epinowcast.org/reference/pcd_dist_name.md),
[`pcd_distributions`](https://primarycensored.epinowcast.org/reference/pcd_distributions.md)
