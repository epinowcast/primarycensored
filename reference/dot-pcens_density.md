# Primary event censored density

Computes the density of the primary event censored delay, the derivative
of
[`pcens_cdf()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.md)
in `x`. This is the contribution of an observation with a zero-width
secondary window.

## Usage

``` r
.pcens_density(object, x, pwindow)
```

## Arguments

- object:

  A `pcens` object as created by
  [`new_pcens()`](https://primarycensored.epinowcast.org/reference/new_pcens.md).

- x:

  Vector of points at which to evaluate the density.

- pwindow:

  Primary event window.

## Value

Vector of densities, not normalised for truncation.

## Details

With `pwindow = 0` this is the delay density. With a uniform primary
event distribution it is \\(F(x) - F(x - pwindow)) / pwindow\\, which
only needs the delay CDF. Otherwise the delay density is integrated
against the primary event density over \[0, pwindow\].
