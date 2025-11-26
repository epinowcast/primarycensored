# Compute primary event censored quantiles

This function inverts the primary event censored CDF to compute
quantiles. It uses numerical optimisation via optim to find the value q
such that
[`pcens_cdf()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.md)
is close to the specified probability. Currently, only the default
numerical inversion method is implemented. Future analytical solutions
may be added.

## Usage

``` r
pcens_quantile(object, p, pwindow, D = Inf, use_numeric = FALSE, ...)
```

## Arguments

- object:

  A `primarycensored` object as created by
  [`new_pcens()`](https://primarycensored.epinowcast.org/dev/reference/new_pcens.md).

- p:

  A vector of probabilities at which to compute the quantiles.

- pwindow:

  Primary event window

- D:

  Maximum delay (truncation point). If finite, the distribution is
  truncated at D. If set to Inf, no truncation is applied. Defaults to
  Inf.

- use_numeric:

  Logical; if TRUE forces the use of numeric inversion even if an
  analytical solution is available (not yet implemented).

- ...:

  Additional arguments to be passed to pdist

## Value

Vector of primary event censored quantiles.

## See also

Low level primary event censored distribution objects and methods
[`new_pcens()`](https://primarycensored.epinowcast.org/dev/reference/new_pcens.md),
[`pcens_cdf()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.md),
[`pcens_cdf.default()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.default.md),
[`pcens_cdf.pcens_pgamma_dunif()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_pgamma_dunif.md),
[`pcens_cdf.pcens_plnorm_dunif()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_plnorm_dunif.md),
[`pcens_cdf.pcens_pweibull_dunif()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_pweibull_dunif.md),
[`pcens_quantile.default()`](https://primarycensored.epinowcast.org/dev/reference/pcens_quantile.default.md)
