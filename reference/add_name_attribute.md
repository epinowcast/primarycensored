# Helper method for custom distributions

[`pprimarycensored()`](https://primarycensored.epinowcast.org/reference/pprimarycensored.md)
and related functions can identify which distributions are provided via
the `pdist` and `dprimary` arguments when those are base R functions
(e.g. `punif`, `dexp`) via the `name` attribute.

## Usage

``` r
add_name_attribute(func, name)
```

## Arguments

- func:

  Function, for example the `p`- or `d`- form of a distribution
  function.

- name:

  Character string, starting with "p" or "d" indicating the underlying
  distribution.

## Value

Function, with a "name" attribute added

## Details

If you need to use a non-base R implementation, but know the
distribution name, you can use this helper function to set it in a way
that will be detected by
[`pprimarycensored()`](https://primarycensored.epinowcast.org/reference/pprimarycensored.md)
and related functions.

This is useful as it enables the automatic use of analytical solutions
for distributions where they exist. You can check which analytical
solutions are available using `methods(pcens_cdf)` and check
distribution names using
[`pcd_dist_name()`](https://primarycensored.epinowcast.org/reference/pcd_dist_name.md).

## See also

Utility functions for working with distributions
[`pcd_dist_name()`](https://primarycensored.epinowcast.org/reference/pcd_dist_name.md),
[`pcd_distributions`](https://primarycensored.epinowcast.org/reference/pcd_distributions.md),
[`pcd_primary_distributions`](https://primarycensored.epinowcast.org/reference/pcd_primary_distributions.md)

## Examples

``` r
dist <- add_name_attribute(pnorm, "hello")
attr(dist, "name")
#> [1] "hello"
```
