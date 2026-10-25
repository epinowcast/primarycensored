# Log of `1 - exp(-exp(x))`

Returns `x` below -37, where `exp(x)` would be lost to rounding.

## Usage

``` r
.log1m_exp_neg_exp(x)
```

## Arguments

- x:

  Numeric vector.

## Value

Numeric vector of `log(1 - exp(-exp(x)))`.
