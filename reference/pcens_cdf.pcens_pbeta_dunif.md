# Method for Beta delay with uniform primary

Analytical solution for the beta distribution, which has support on
\\\[0, 1\]\\. A non-central beta (`ncp` not zero) uses the numerical
[`pcens_cdf.default()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.default.md)
method.

## Usage

``` r
# S3 method for class 'pcens_pbeta_dunif'
pcens_cdf(object, q, pwindow, use_numeric = FALSE)
```

## Arguments

- object:

  A `primarycensored` object as created by
  [`new_pcens()`](https://primarycensored.epinowcast.org/reference/new_pcens.md).

- q:

  Vector of quantiles

- pwindow:

  Primary event window. Use `pwindow = 0` for an exactly observed
  primary event, in which case the delay CDF is used directly.

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
[`pcens_cdf.pcens_pchisq_dunif()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_pchisq_dunif.md),
[`pcens_cdf.pcens_pdiscretehazard()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_pdiscretehazard.md),
[`pcens_cdf.pcens_pdiscretestep()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_pdiscretestep.md),
[`pcens_cdf.pcens_pexp_dunif()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_pexp_dunif.md),
[`pcens_cdf.pcens_pgamma_dunif()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_pgamma_dunif.md),
[`pcens_cdf.pcens_pgengamma.orig_dunif()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_pgengamma.orig_dunif.md),
[`pcens_cdf.pcens_pgengamma_dunif()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_pgengamma_dunif.md),
[`pcens_cdf.pcens_pinvchisq_dunif()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_pinvchisq_dunif.md),
[`pcens_cdf.pcens_pinvgamma_dunif()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_pinvgamma_dunif.md),
[`pcens_cdf.pcens_plnorm_dunif()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_plnorm_dunif.md),
[`pcens_cdf.pcens_pnorm_dunif()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_pnorm_dunif.md),
[`pcens_cdf.pcens_ppareto_dunif()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_ppareto_dunif.md),
[`pcens_cdf.pcens_pweibull_dunif()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_pweibull_dunif.md),
[`pcens_pmf()`](https://primarycensored.epinowcast.org/reference/pcens_pmf.md),
[`pcens_pmf.default()`](https://primarycensored.epinowcast.org/reference/pcens_pmf.default.md),
[`pcens_quantile()`](https://primarycensored.epinowcast.org/reference/pcens_quantile.md),
[`pcens_quantile.default()`](https://primarycensored.epinowcast.org/reference/pcens_quantile.default.md),
[`update.pcens()`](https://primarycensored.epinowcast.org/reference/update.pcens.md)

## Examples

``` r
pcens_obj <- new_pcens(
  pdist = pbeta,
  dprimary = dunif,
  primary_args = list(min = 0, max = 1),
  shape1 = 2,
  shape2 = 3
)
pcens_cdf(pcens_obj, q = c(0.2, 0.6, 1.5), pwindow = 0.5)
#> [1] 0.025984 0.435300 1.000000
```
