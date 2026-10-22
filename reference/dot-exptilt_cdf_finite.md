# Exponentially tilted CDF at finite points

Chooses the form at each point as in the Stan function
`primarycensored_exptilt_form()`, see
[pcens_cdf_exptilt](https://primarycensored.epinowcast.org/reference/pcens_cdf_exptilt.md).

## Usage

``` r
.exptilt_cdf_finite(object, q, pwindow, rho)
```

## Arguments

- object:

  A `primarycensored` object as created by
  [`new_pcens()`](https://primarycensored.epinowcast.org/reference/new_pcens.md).

- q:

  Vector of quantiles

- pwindow:

  Primary event window. Use `pwindow = 0` for an exactly observed
  primary event, in which case the delay CDF is used directly.

- rho:

  The tilt, the `r` of the exponential growth primary.

## Value

Vector of CDFs, clamped to \[0, 1\].
