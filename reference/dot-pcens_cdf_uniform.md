# Uniform primary event censored CDF from terms at the window ends

For a delay with mean \\E\\ and a uniform primary event window of width
\\w\\, the CDF at \\d\\ is \\(G(d) - G(q)) / w\\ with \\q = \max(d - w,
0)\\ and \\G(t) = t F(t) - E \tilde F(t)\\. \\F\\ is the delay CDF and
\\\tilde F\\ is the CDF of the partial expectation distribution.

## Usage

``` r
.pcens_cdf_uniform(q, pwindow, terms_fn, upper_fn, delay_cdf, params)
```

## Arguments

- q:

  Vector of quantiles

- pwindow:

  Primary event window. Use `pwindow = 0` for an exactly observed
  primary event, in which case the delay CDF is used directly.

- terms_fn:

  Function of delays greater than or equal to zero and `params`
  returning \\G\\.

- upper_fn:

  Function of delays greater than or equal to zero and `params`
  returning \\H(t) = t S(t) - E \tilde S(t)\\, where \\S\\ and \\\tilde
  S\\ are the upper tails of \\F\\ and \\\tilde F\\.

- delay_cdf:

  Function of delays and `params` returning the delay CDF.

- params:

  Named list of delay parameters, including the delay mean \\E\\ as
  `mean`.

## Value

Numeric vector of CDF values in \[0, 1\].

## Details

Delays at or below zero give 0, delays of `Inf` give 1 and missing
delays are an error. Zero-width windows give the delay CDF. `q`,
`pwindow` and the elements of `params` are recycled against each other.
