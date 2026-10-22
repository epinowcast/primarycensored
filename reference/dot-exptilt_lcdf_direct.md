# Direct form of the exponentially tilted log CDF

The expression of
[pcens_cdf_exptilt](https://primarycensored.epinowcast.org/reference/pcens_cdf_exptilt.md)
on the log scale.

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
