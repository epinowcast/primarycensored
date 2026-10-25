# Primary event censored CDF for an exponentially tilted primary

Shared implementation of the
[pcens_cdf_exptilt](https://primarycensored.epinowcast.org/reference/pcens_cdf_exptilt.md)
methods.

## Usage

``` r
.pcens_cdf_exptilt(
  object,
  q,
  pwindow,
  use_numeric = FALSE,
  min_q = 0L,
  min_xw = Inf
)
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

- use_numeric:

  Logical, if TRUE forces use of numeric integration even for
  distributions with analytical solutions. This is primarily useful for
  testing purposes or for settings where the analytical solution breaks
  down.

- min_q:

  Number of quantiles below which the numerical method is used, for
  delays whose transform has a fixed cost. The default 0 always uses the
  closed forms.

- min_xw:

  Largest \\\|\rho\| w\\ at which `min_q` applies. The numerical method
  loses accuracy for a larger tilt times the window, so the closed forms
  are used there for any number of quantiles.

## Value

Vector of CDFs.
