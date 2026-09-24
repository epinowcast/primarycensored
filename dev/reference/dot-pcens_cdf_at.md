# Primary event censored CDF at a truncation point

Returns the primary event censored CDF, before truncation, at a single
truncation point, reusing an already computed value where possible.

## Usage

``` r
.pcens_cdf_at(object, bound, pwindow, points, cdfs, inf_value)
```

## Arguments

- object:

  A `pcens` object.

- bound:

  Numeric truncation point (`L` or `D`).

- pwindow:

  Primary event window.

- points:

  Numeric vector of points at which `cdfs` was computed.

- cdfs:

  Numeric vector of CDF values at `points`.

- inf_value:

  CDF value to return when `bound` is infinite.

## Value

A single numeric CDF value.
