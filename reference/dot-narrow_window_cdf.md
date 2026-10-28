# Mean of the delay CDF over narrow uniform primary windows

Uses a 5 point Gauss-Legendre rule over \[d - w, d\].

## Usage

``` r
.narrow_window_cdf(d, w, delay_cdf, params)
```

## Arguments

- d:

  Finite delays above `lower`.

- w:

  Primary event windows greater than zero, one per delay.

- delay_cdf:

  Function of delays and `params` returning the delay CDF.

- params:

  Named list of delay parameters, including the delay mean \\E\\ as
  `mean`.

## Value

Numeric vector of CDF values, one per delay.
