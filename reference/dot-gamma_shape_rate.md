# Gamma delay parameters of a pcens object

The rate is 1 if neither `rate` nor `scale` is given, as in
[`stats::pgamma()`](https://rdrr.io/r/stats/GammaDist.html).

## Usage

``` r
.gamma_shape_rate(object)
```

## Arguments

- object:

  A `pcens` object.

## Value

A list with `shape` and `rate`.
