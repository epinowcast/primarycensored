# Truncated Gumbel primary censored CDF for a finite positive window

Uses the series where `series` is `TRUE` and its estimated error is
below `.gumbel_tol`, and the numerical method otherwise.

## Usage

``` r
.gumbel_cdf(object, q, pwindow, mu, beta, series = FALSE)
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

- mu:

  Location of the Gumbel distribution before truncation, on the same
  scale as `x`.

- beta:

  Scale of the Gumbel distribution before truncation, positive.

- series:

  Logical; whether to try the series.

## Value

Vector of CDFs.
