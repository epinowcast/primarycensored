# Name a stats function found in the distribution registries

Takes the C routine called at the end of the body of a `stats` function
(for example `C_pgamma` for
[`stats::pgamma()`](https://rdrr.io/r/stats/GammaDist.html)) as the
candidate name. The name is returned if it is in
[pcd_distributions](https://primarycensored.epinowcast.org/reference/pcd_distributions.md)
or
[pcd_primary_distributions](https://primarycensored.epinowcast.org/reference/pcd_primary_distributions.md)
and `func` is identical to the `stats` function of that name.

## Usage

``` r
.registry_name(func)
```

## Arguments

- func:

  Function, for example the `p`- or `d`- form of a distribution
  function.

## Value

The registry name of `func`, or `NULL` if it is not found.
