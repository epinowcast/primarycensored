# Small delay form of the exponentially tilted log CDF

See
[pcens_cdf_exptilt](https://primarycensored.epinowcast.org/reference/pcens_cdf_exptilt.md).

## Usage

``` r
.exptilt_lcdf_small_delay(object, q, pwindow, rho)
```

## Arguments

- object:

  A `pcens` object.

- q:

  Vector of quantiles

- pwindow:

  Primary event window. Use `pwindow = 0` for an exactly observed
  primary event, in which case the delay CDF is used directly.

## Value

Vector of log CDFs.
