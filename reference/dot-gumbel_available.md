# Whether the truncated Gumbel series can be used

The choice depends on the parameters only. The rounding error is about
the machine precision times the size of the log transform at the largest
tilt times the terms \\e^{s_0}\\, so the series is skipped where that is
above `.gumbel_screen_tol`.

## Usage

``` r
.gumbel_available(object, mu, beta)
```

## Arguments

- object:

  A `primarycensored` object as created by
  [`new_pcens()`](https://primarycensored.epinowcast.org/reference/new_pcens.md).

- mu:

  Location of the Gumbel distribution before truncation, on the same
  scale as `x`.

- beta:

  Scale of the Gumbel distribution before truncation, positive.

## Value

`TRUE` if the series can be used.
