# Number of terms of the truncated Gumbel series

Number of terms of the truncated Gumbel series

## Usage

``` r
.gumbel_n_terms(log_s0)
```

## Arguments

- log_s0:

  Log of \\s_0 = e^{\mu / \beta}\\, that is `mu / beta`.

## Value

The smallest integer `N`, at least 2, for which \\s_0^n / n!\\ is below
`1e-18 * max(s0, 1)` after `N` terms.
