# Check whether two names refer to the same registry primary distribution

Check whether two names refer to the same registry primary distribution

## Usage

``` r
.same_primary(d_name, p_name)
```

## Arguments

- d_name, p_name:

  Names of a primary density and CDF. Each may be a name, alias, density
  or CDF name from
  [pcd_primary_distributions](https://primarycensored.epinowcast.org/dev/reference/pcd_primary_distributions.md).

## Value

`TRUE` if both names match the same row of
[pcd_primary_distributions](https://primarycensored.epinowcast.org/dev/reference/pcd_primary_distributions.md),
otherwise `FALSE`.
