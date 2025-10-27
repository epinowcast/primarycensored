# Method for Log-Normal delay with uniform primary

Method for Log-Normal delay with uniform primary

## Usage

``` r
# S3 method for class 'pcens_plnorm_dunif'
pcens_cdf(object, q, pwindow, use_numeric = FALSE)
```

## Arguments

- object:

  A `primarycensored` object as created by
  [`new_pcens()`](https://primarycensored.epinowcast.org/reference/new_pcens.md).

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

## See also

Low level primary event censored distribution objects and methods
[`new_pcens()`](https://primarycensored.epinowcast.org/reference/new_pcens.md),
[`pcens_cdf()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.md),
[`pcens_cdf.default()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.default.md),
[`pcens_cdf.pcens_pgamma_dunif()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_pgamma_dunif.md),
[`pcens_cdf.pcens_pweibull_dunif()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_pweibull_dunif.md),
[`pcens_quantile()`](https://primarycensored.epinowcast.org/reference/pcens_quantile.md),
[`pcens_quantile.default()`](https://primarycensored.epinowcast.org/reference/pcens_quantile.default.md)
