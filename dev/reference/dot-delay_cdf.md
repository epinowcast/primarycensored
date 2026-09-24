# Delay CDF of a pcens object

Evaluates the delay distribution CDF of a `pcens` object with its stored
parameters. This is the primary event censored CDF when `pwindow = 0`.

## Usage

``` r
.delay_cdf(object, q)
```

## Arguments

- object:

  A `pcens` object as created by
  [`new_pcens()`](https://primarycensored.epinowcast.org/dev/reference/new_pcens.md).

- q:

  Vector of quantiles.

## Value

Vector of delay CDF values.
