# Sum the truncated Gumbel series

Sum the truncated Gumbel series

## Usage

``` r
.gumbel_combine(log_delta, log_f_y, scale_t, q, pwindow, mu, beta)
```

## Arguments

- log_delta:

  Matrix of the log of \\\Delta_n(q)\\, with a row per `q` and a column
  per tilt, \\n = 0, \dots, N\\.

- log_f_y:

  Log delay CDF at `q - pwindow`.

- scale_t:

  Matrix, shaped as `log_delta`, of the size of the log transforms, for
  the rounding estimate.

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
