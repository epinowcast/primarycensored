# Log differences of the truncated Gumbel window

With \\s(x) = \exp\\-(x - \mu) / \beta\\\\,
`.tgumbel_log_delta_window()` gives the log of \\s(min) - s(max)\\ and
`.tgumbel_log_deltas()` the logs of \\s(min) - s(x)\\ and \\s(x) -
s(max)\\.

## Usage

``` r
.tgumbel_log_delta_window(min, max, mu, beta)

.tgumbel_log_deltas(x, min, max, mu, beta)
```

## Arguments

- min:

  Minimum value of the distribution range. Default is 0.

- max:

  Maximum value of the distribution range. Default is 1.

- mu:

  Location of the Gumbel distribution before truncation, on the same
  scale as `x`.

- beta:

  Scale of the Gumbel distribution before truncation, positive.

## Value

A numeric vector for `.tgumbel_log_delta_window()`, and a list with
`lower` and `upper` for `.tgumbel_log_deltas()`.
