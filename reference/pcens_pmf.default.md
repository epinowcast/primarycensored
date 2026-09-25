# Default method for computing primary event censored PMF

Computes the PMF by differencing
[`pcens_cdf()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.md)
at `x` and `min(x + swindow, D)`, normalised over \[L, D\]. Where
`swindow = 0` the primary event censored density at `x` is returned
instead. See
[`dprimarycensored()`](https://primarycensored.epinowcast.org/reference/dprimarycensored.md)
for the details.

## Usage

``` r
# Default S3 method
pcens_pmf(object, x, pwindow, swindow = 1, L = -Inf, D = Inf, log = FALSE, ...)
```

## Arguments

- object:

  A `pcens` object as created by
  [`new_pcens()`](https://primarycensored.epinowcast.org/reference/new_pcens.md).

- x:

  Vector of quantiles

- pwindow:

  Primary event window. Use `pwindow = 0` for an exactly observed
  primary event, in which case the delay CDF is used directly.

- swindow:

  Secondary event window (default: 1). Use `swindow = 0` for an exactly
  observed secondary event, in which case a density is returned rather
  than a probability (see Details).

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
[`new_pcens()`](https://primarycensored.epinowcast.org/reference/new_pcens.md),
[`pcens_cdf()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.md),
[`pcens_cdf.default()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.default.md),
[`pcens_cdf.pcens_pdiscretehazard()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_pdiscretehazard.md),
[`pcens_cdf.pcens_pdiscretestep()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_pdiscretestep.md),
[`pcens_cdf.pcens_pgamma_dunif()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_pgamma_dunif.md),
[`pcens_cdf.pcens_pgengamma.orig_dunif()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_pgengamma.orig_dunif.md),
[`pcens_cdf.pcens_pgengamma_dunif()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_pgengamma_dunif.md),
[`pcens_cdf.pcens_plnorm_dunif()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_plnorm_dunif.md),
[`pcens_cdf.pcens_pweibull_dunif()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_pweibull_dunif.md),
[`pcens_pmf()`](https://primarycensored.epinowcast.org/reference/pcens_pmf.md),
[`pcens_quantile()`](https://primarycensored.epinowcast.org/reference/pcens_quantile.md),
[`pcens_quantile.default()`](https://primarycensored.epinowcast.org/reference/pcens_quantile.default.md),
[`update.pcens()`](https://primarycensored.epinowcast.org/reference/update.pcens.md)

## Examples

``` r
obj <- new_pcens(
  pdist = pgamma, dprimary = dunif,
  primary_args = list(min = 0, max = 1),
  shape = 3, scale = 2
)
pcens_pmf(obj, x = 0:9, pwindow = 1, D = 10)
#>  [1] 0.004551244 0.045675031 0.105784553 0.144941751 0.157178437 0.149641077
#>  [7] 0.131158433 0.108576455 0.086202811 0.066290209
```
