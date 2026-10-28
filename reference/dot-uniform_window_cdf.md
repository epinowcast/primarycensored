# Uniform primary CDF for finite delays and positive windows

Uniform primary CDF for finite delays and positive windows

## Usage

``` r
.uniform_window_cdf(d, w, terms_fn, upper_fn, delay_cdf, params, share, lower)
```

## Arguments

- d:

  Finite delays above `lower`.

- w:

  Primary event windows greater than zero, one per delay.

- terms_fn:

  Function of delays greater than or equal to `lower` and `params`
  returning \\G\\.

- upper_fn:

  Function of delays greater than or equal to `lower` and `params`
  returning \\H(t) = t S(t) - E \tilde S(t)\\, where \\S\\ and \\\tilde
  S\\ are the upper tails of \\F\\ and \\\tilde F\\.

- delay_cdf:

  Function of delays and `params` returning the delay CDF.

- params:

  Named list of delay parameters, including the delay mean \\E\\ as
  `mean`.

- share:

  Whether terms may be shared between window ends, as in
  [`.uniform_endpoint_diff()`](https://primarycensored.epinowcast.org/reference/dot-uniform_endpoint_diff.md).
  Only valid for a scalar window and scalar `params`.

- lower:

  Lower end of the delay support, \\L\\: 0, or `-Inf` for a delay on the
  reals.

## Value

Numeric vector of CDF values, not clamped to \[0, 1\].
