# Log-scale arithmetic helpers

Vectorised. `.log_diff_exp()` gives `-Inf` rather than `NaN` for a zero
or negative difference.

## Usage

``` r
.log1m_exp(x)

.log_diff_exp(a, b)

.log_sum_exp(a, b)
```

## Arguments

- a, b, x:

  Numeric vectors on the log scale.

## Value

`log(1 - exp(x))`, `log(exp(a) - exp(b))` and `log(exp(a) + exp(b))` for
`.log1m_exp()`, `.log_diff_exp()` and `.log_sum_exp()`.
