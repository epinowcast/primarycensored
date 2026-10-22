# Endpoints and the positions of `q` and `q - pwindow` among them

Each term of the exponentially tilted CDF depends on one endpoint, so it
is computed once per unique endpoint. Endpoints at or below a finite
`lower` share one entry.

## Usage

``` r
.exptilt_layout(q, pwindow, lower)
```

## Arguments

- q:

  Vector of quantiles

- pwindow:

  Primary event window. Use `pwindow = 0` for an exactly observed
  primary event, in which case the delay CDF is used directly.

- lower:

  Lower end of the support of the delay, 0 or `-Inf`.

## Value

List with the unique `endpoints` and the positions `q` and `y` of `q`
and `q - pwindow` among them.
