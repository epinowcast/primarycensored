# Log moments of a gamma delay about a point

The log of \\G_k(t) = \int_0^t (t - u)^k f(u) du\\ for `k = 1, 2, 3`,
from gamma CDFs with the shape raised by `k`.

## Usage

``` r
.gamma_moments(t, shape, rate)
```

## Arguments

- t:

  Numeric vector of points.

- shape, rate:

  Gamma delay parameters.

## Value

As for
[`.pcens_tilt_moments()`](https://primarycensored.epinowcast.org/reference/pcens_tilt.md).
