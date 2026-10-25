# Numerical truncated Gumbel primary censored CDF

The integral is taken in \\u = s(z) - s(w)\\, with \\s(z) = \exp\\(\mu -
z) / \beta\\\\, for \\\mu\\ at or above the window width \\w\\, and in
\\z\\ otherwise, so a narrow spike of the density is resolved.

## Usage

``` r
.gumbel_numeric(object, q, pwindow, mu, beta)
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

- mu:

  Location of the Gumbel distribution before truncation, on the same
  scale as `x`.

- beta:

  Scale of the Gumbel distribution before truncation, positive.

## Value

Vector of CDFs at finite `q`.
