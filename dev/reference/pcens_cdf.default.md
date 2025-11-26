# Default method for computing primary event censored CDF

This method serves as a fallback for combinations of delay and primary
event distributions that don't have specific implementations. It uses a
numeric integration method.

## Usage

``` r
# Default S3 method
pcens_cdf(object, q, pwindow, use_numeric = FALSE)
```

## Arguments

- object:

  A `primarycensored` object as created by
  [`new_pcens()`](https://primarycensored.epinowcast.org/dev/reference/new_pcens.md).

- q:

  Vector of quantiles

- pwindow:

  Primary event window

- use_numeric:

  Logical, if TRUE forces use of numeric integration even for
  distributions with analytical solutions. This is primarily useful for
  testing purposes or for settings where the analytical solution breaks
  down.

## Value

Vector of computed primary event censored CDFs

## Details

This method implements the numerical integration approach for computing
the primary event censored CDF. It uses the same mathematical
formulation as described in the details section of
[`pprimarycensored()`](https://primarycensored.epinowcast.org/dev/reference/pprimarycensored.md),
but applies numerical integration instead of analytical solutions.

## See also

[`pprimarycensored()`](https://primarycensored.epinowcast.org/dev/reference/pprimarycensored.md)
for the mathematical details of the primary event censored CDF
computation.

Low level primary event censored distribution objects and methods
[`new_pcens()`](https://primarycensored.epinowcast.org/dev/reference/new_pcens.md),
[`pcens_cdf()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.md),
[`pcens_cdf.pcens_pgamma_dunif()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_pgamma_dunif.md),
[`pcens_cdf.pcens_plnorm_dunif()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_plnorm_dunif.md),
[`pcens_cdf.pcens_pweibull_dunif()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_pweibull_dunif.md),
[`pcens_quantile()`](https://primarycensored.epinowcast.org/dev/reference/pcens_quantile.md),
[`pcens_quantile.default()`](https://primarycensored.epinowcast.org/dev/reference/pcens_quantile.default.md)

## Examples

``` r
# Create a primarycensored object with gamma delay and uniform primary
pcens_obj <- new_pcens(
  pdist = pgamma,
  dprimary = dunif,
  dprimary_args = list(min = 0, max = 1),
  shape = 3,
  scale = 2
)

# Compute CDF for a single value
pcens_cdf(pcens_obj, q = 9, pwindow = 1)
#> [1] 0.7955788

# Compute CDF for multiple values
pcens_cdf(pcens_obj, q = c(4, 6, 8), pwindow = 1)
#> [1] 0.2564303 0.5178596 0.7221287
```
