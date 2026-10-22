# Log moments about a point from log partial moments

The log of \\G_k(t) = \int_0^t (t - u)^k f(u) du\\ for `k = 1, 2, 3`
from the partial moments \\m_k(t) = \int_0^t u^k f(u) du\\, with \\G_1 =
t m_0 - m_1\\, \\G_2 = t G_1 - (t m_1 - m_2)\\ and \\G_3 = t G_2 - (t (t
m_1 - m_2) - (t m_2 - m_3))\\, each a difference of positive integrals.

## Usage

``` r
.moments_from_partial(log_t, log_m0, log_m1, log_m2, log_m3)
```

## Arguments

- log_t:

  Log of the points.

- log_m0, log_m1, log_m2, log_m3:

  Log partial moments at each point.

## Value

As for
[`.pcens_tilt_moments()`](https://primarycensored.epinowcast.org/reference/pcens_tilt.md).
