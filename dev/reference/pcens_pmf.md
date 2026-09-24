# Compute primary event censored PMF

Computes the primary event censored PMF for a `pcens` object as created
by
[`new_pcens()`](https://primarycensored.epinowcast.org/dev/reference/new_pcens.md).
Secondary event windows and truncation are handled as in
[`dprimarycensored()`](https://primarycensored.epinowcast.org/dev/reference/dprimarycensored.md).

## Usage

``` r
pcens_pmf(object, x, pwindow, swindow = 1, L = -Inf, D = Inf, log = FALSE, ...)
```

## Arguments

- object:

  A `pcens` object as created by
  [`new_pcens()`](https://primarycensored.epinowcast.org/dev/reference/new_pcens.md).

- x:

  Vector of quantiles

- pwindow:

  Primary event window

- swindow:

  Secondary event window (default: 1)

- L:

  Minimum delay (lower truncation point). Defaults to `-Inf`, meaning no
  left truncation. For any finite value of L the distribution is
  left-truncated at L.

- D:

  Maximum delay (upper truncation point). If finite, the distribution is
  truncated at D. If set to Inf, no upper truncation is applied.
  Defaults to Inf.

- log:

  Logical; if TRUE, probabilities p are given as log(p)

- ...:

  Additional arguments passed to methods.

## Value

Vector of primary event censored PMFs, normalized over \[L, D\] if
truncation is applied

## See also

Low level primary event censored distribution objects and methods
[`new_pcens()`](https://primarycensored.epinowcast.org/dev/reference/new_pcens.md),
[`pcens_cdf()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.md),
[`pcens_cdf.default()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.default.md),
[`pcens_cdf.pcens_pdiscretehazard()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_pdiscretehazard.md),
[`pcens_cdf.pcens_pdiscretestep()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_pdiscretestep.md),
[`pcens_cdf.pcens_pgamma_dunif()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_pgamma_dunif.md),
[`pcens_cdf.pcens_pgengamma.orig_dunif()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_pgengamma.orig_dunif.md),
[`pcens_cdf.pcens_pgengamma_dunif()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_pgengamma_dunif.md),
[`pcens_cdf.pcens_plnorm_dunif()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_plnorm_dunif.md),
[`pcens_cdf.pcens_pweibull_dunif()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_pweibull_dunif.md),
[`pcens_pmf.default()`](https://primarycensored.epinowcast.org/dev/reference/pcens_pmf.default.md),
[`pcens_quantile()`](https://primarycensored.epinowcast.org/dev/reference/pcens_quantile.md),
[`pcens_quantile.default()`](https://primarycensored.epinowcast.org/dev/reference/pcens_quantile.default.md),
[`update.pcens()`](https://primarycensored.epinowcast.org/dev/reference/update.pcens.md)
