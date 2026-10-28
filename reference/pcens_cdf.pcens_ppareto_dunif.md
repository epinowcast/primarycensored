# Method for Pareto delay with uniform primary

Analytical solution for the type I Pareto distribution of
[`extraDistr::ppareto()`](https://rdrr.io/pkg/extraDistr/man/Pareto.html),
with shape `a` and scale (minimum) `b`. The `extraDistr` package is
needed to create the `pcens` object. The solution holds for every shape,
including shapes at or below 1 where the mean is infinite. See
[`vignette("analytic-solutions")`](https://primarycensored.epinowcast.org/articles/analytic-solutions.md)
for the derivation.

## Usage

``` r
# S3 method for class 'pcens_ppareto_dunif'
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
[`pcens_cdf.pcens_pbeta_dunif()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_pbeta_dunif.md),
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
[`pcens_cdf.pcens_pweibull_dunif()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_pweibull_dunif.md),
[`pcens_pmf()`](https://primarycensored.epinowcast.org/reference/pcens_pmf.md),
[`pcens_pmf.default()`](https://primarycensored.epinowcast.org/reference/pcens_pmf.default.md),
[`pcens_quantile()`](https://primarycensored.epinowcast.org/reference/pcens_quantile.md),
[`pcens_quantile.default()`](https://primarycensored.epinowcast.org/reference/pcens_quantile.default.md),
[`update.pcens()`](https://primarycensored.epinowcast.org/reference/update.pcens.md)

## Examples

``` r
pcens_obj <- new_pcens(
  pdist = extraDistr::ppareto,
  dprimary = dunif,
  primary_args = list(min = 0, max = 1),
  a = 2.5,
  b = 0.5
)
pcens_cdf(pcens_obj, q = c(1, 4, 8), pwindow = 1)
#> [1] 0.2845178 0.9920509 0.9988450
```
