# Log of `abs(exp(rho * pwindow) - 1)`

Log of `abs(exp(rho * pwindow) - 1)`

## Usage

``` r
.exptilt_log_den(rho, pwindow)
```

## Arguments

- rho:

  The tilt, the `r` of the exponential growth primary.

- pwindow:

  Primary event window. Use `pwindow = 0` for an exactly observed
  primary event, in which case the delay CDF is used directly.

## Value

The log of the normalising constant of the tilted window.
