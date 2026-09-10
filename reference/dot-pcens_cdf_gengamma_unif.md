# Analytical primary event censored CDF for the generalised gamma

Shared implementation for the Stacy parameterisation used by both
generalised gamma
[`pcens_cdf()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.md)
methods.

## Usage

``` r
.pcens_cdf_gengamma_unif(q, pwindow, shape, scale, k)
```

## Arguments

- q:

  Vector of quantiles

- pwindow:

  Primary event window

- shape, scale, k:

  Generalised gamma parameters in the Stacy parameterisation of
  [`flexsurv::pgengamma.orig()`](http://chjackson.github.io/flexsurv-dev/reference/GenGamma.orig.md).

## Value

Vector of computed primary event censored CDFs
