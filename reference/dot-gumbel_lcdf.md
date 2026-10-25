# Truncated Gumbel log CDF from the series

Truncated Gumbel log CDF from the series

## Usage

``` r
.gumbel_lcdf(object, q, pwindow, mu, beta)
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

## Value

List with the log CDFs `log_cdf` and their estimated relative errors
`error`.
