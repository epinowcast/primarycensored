# Tilt transforms of delay distributions

Internal generics for the delay part of the analytical solutions with an
exponentially tilted primary event window, see
[pcens_cdf_exptilt](https://primarycensored.epinowcast.org/reference/pcens_cdf_exptilt.md).
They are built from the tilt transform \\T_f(\xi; \tau) = \int e^{\xi u}
f(u) du\\ up to \\\tau\\, from 0 for delays on the non-negative reals. A
delay is added with a method for each generic, dispatched on the delay
class of a `pcens` object. `.pcens_tilt_fits()` and `.pcens_tilt_pair()`
have defaults. The Stan equivalents are `check_for_tilt_transform()`,
`log_tilt_transform_pair()` and `primarycensored_tilt_moments()`.

## Usage

``` r
.pcens_tilt_transform(object, t, xi, upper = FALSE)

.pcens_tilt_available(object, xi)

.pcens_tilt_lower(object)

.pcens_tilt_moments(object, t)

.pcens_tilt_fits(object, xi, t, pwindow = 0)

.pcens_tilt_pair(object, t, xi)
```

## Arguments

- object:

  A `pcens` object.

- t:

  Numeric vector of points.

- xi:

  The tilt, a single number or, for `.pcens_tilt_transform()`, a vector
  of the length of `t`. A window with tilt \\\rho\\ needs `xi = -rho`.

- upper:

  Logical. If `TRUE`, the transform over `(t, Inf)` rather than up to
  `t`.

- pwindow:

  The primary event window.

## Value

- `.pcens_tilt_transform()`: the log of the transform at each `t`.

- `.pcens_tilt_available()`: whether the transform is closed form and
  the tilted delay exists for `xi`. If not, the numerical method is
  used.

- `.pcens_tilt_fits()`: `TRUE` at each `t` where the transform is used
  for a window `pwindow`, otherwise the numerical method. Only called
  where `.pcens_tilt_available()` is `TRUE`. The default is `TRUE`.

- `.pcens_tilt_pair()`: a matrix with columns `lower` and `upper`, the
  log transform up to and beyond each `t`. The default calls
  `.pcens_tilt_transform()` twice and keeps the `log_loss` attribute of
  the lower transform.

- `.pcens_tilt_lower()`: the lower end of the support, 0 or `-Inf`.

- `.pcens_tilt_moments()`: a matrix with columns `G1`, `G2` and `G3`,
  the log of \\G_k(t) = \int (t - u)^k f(u) du\\ up to `t`, `-Inf` below
  the support.
