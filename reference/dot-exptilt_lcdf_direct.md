# Direct form of the exponentially tilted log CDF

The expression of
[pcens_cdf_exptilt](https://primarycensored.epinowcast.org/reference/pcens_cdf_exptilt.md)
on the log scale. The result is `NaN` where the transform at either
endpoint is, and where the transform has no upper tail and the terms of
the numerator exceed the smaller of the CDF and its survival function by
more than `.exptilt_max_loss`.

## Usage

``` r
.exptilt_lcdf_direct(object, q, pwindow, rho, lower)
```

## Arguments

- object:

  A `pcens` object.

- q:

  Vector of quantiles

- pwindow:

  Primary event window. Use `pwindow = 0` for an exactly observed
  primary event, in which case the delay CDF is used directly.

- lower:

  Lower end of the support of the delay, 0 or `-Inf`.

## Value

Vector of log CDFs.
