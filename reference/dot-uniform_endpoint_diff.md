# Difference of terms between the ends of uniform primary windows

Returns `fn(d, params) - fn(lo, params)`. With `share = TRUE` and at
least 64 delays, lower ends that equal a delay reuse the terms at that
delay, so `fn` is evaluated once per distinct window end. This is tried
only when at least 3 of 6 sampled lower ends match a delay, which skips
the full match when few ends coincide, and used only when at most half
as many unmatched distinct lower ends as delays remain.

## Usage

``` r
.uniform_endpoint_diff(fn, d, lo, params, share)
```

## Arguments

- fn:

  Function of delays and `params`, such as `terms_fn` or `upper_fn`.

- d:

  Positive finite delays.

- lo:

  Lower window ends, one per delay.

- params:

  Named list of delay parameters, including the delay mean \\E\\ as
  `mean`.

- share:

  Whether terms may be shared between window ends, as in
  `.uniform_endpoint_diff()`. Only valid for a scalar window and scalar
  `params`.

## Value

Numeric vector of differences, one per delay.
