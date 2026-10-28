# Evaluate `z * pnorm(z) + dnorm(z)`, stably for very negative `z`

Evaluate `z * pnorm(z) + dnorm(z)`, stably for very negative `z`

## Usage

``` r
.norm_shortfall(z)
```

## Arguments

- z:

  Numeric vector of standardised values.

## Value

Vector of `z * pnorm(z) + dnorm(z)`.
