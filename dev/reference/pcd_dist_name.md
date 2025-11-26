# Get distribution function cdf or pdf name

Get distribution function cdf or pdf name

## Usage

``` r
pcd_dist_name(name, type = c("delay", "primary"))
```

## Arguments

- name:

  String. Distribution name or alias

- type:

  String. "delay" or "primary" corresponding to the type of distribution
  to use as the look up. If delay then
  [`pcd_distributions()`](https://primarycensored.epinowcast.org/dev/reference/pcd_distributions.md)
  is used, if primary then
  [`pcd_primary_distributions()`](https://primarycensored.epinowcast.org/dev/reference/pcd_primary_distributions.md)
  is used.

## Value

String distribution function name or NA if no base R implementation

## See also

Utility functions for working with distributions
[`add_name_attribute()`](https://primarycensored.epinowcast.org/dev/reference/add_name_attribute.md),
[`pcd_distributions`](https://primarycensored.epinowcast.org/dev/reference/pcd_distributions.md),
[`pcd_primary_distributions`](https://primarycensored.epinowcast.org/dev/reference/pcd_primary_distributions.md)

## Examples

``` r
pcd_dist_name("lnorm")
#> [1] "plnorm"
pcd_dist_name("lognormal")
#> [1] "plnorm"
pcd_dist_name("gamma")
#> [1] "pgamma"
pcd_dist_name("weibull")
#> [1] "pweibull"
pcd_dist_name("exp")
#> [1] "pexp"
pcd_dist_name("unif", type = "primary")
#> [1] "dunif"
pcd_dist_name("expgrowth", type = "primary")
#> [1] "dexpgrowth"
```
