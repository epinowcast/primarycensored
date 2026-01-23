# Get distribution stan ID by name

Get distribution stan ID by name

## Usage

``` r
pcd_stan_dist_id(name, type = c("delay", "primary"))
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

Numeric distribution ID

## See also

Tools for working with package Stan functions
[`pcd_load_stan_functions()`](https://primarycensored.epinowcast.org/dev/reference/pcd_load_stan_functions.md),
[`pcd_stan_files()`](https://primarycensored.epinowcast.org/dev/reference/pcd_stan_files.md),
[`pcd_stan_function_deps()`](https://primarycensored.epinowcast.org/dev/reference/pcd_stan_function_deps.md),
[`pcd_stan_functions()`](https://primarycensored.epinowcast.org/dev/reference/pcd_stan_functions.md),
[`pcd_stan_path()`](https://primarycensored.epinowcast.org/dev/reference/pcd_stan_path.md)

## Examples

``` r
pcd_stan_dist_id("lnorm")
#> [1] 1
pcd_stan_dist_id("lognormal")
#> [1] 1
pcd_stan_dist_id("gamma")
#> [1] 2
pcd_stan_dist_id("weibull")
#> [1] 3
pcd_stan_dist_id("exp")
#> [1] 4
pcd_stan_dist_id("unif", type = "primary")
#> [1] 1
```
