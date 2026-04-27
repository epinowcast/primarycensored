# Default method for computing primary event censored quantiles

This method inverts the primary event censored CDF by root-finding via
[`stats::uniroot()`](https://rdrr.io/r/stats/uniroot.html) with
`extendInt = "upX"`. The censored CDF is monotone in `q`, so
[`stats::uniroot()`](https://rdrr.io/r/stats/uniroot.html) extends its
starting bracket outward as needed and handles infinite `L` or `D`
without special casing.

## Usage

``` r
# Default S3 method
pcens_quantile(
  object,
  p,
  pwindow,
  L = -Inf,
  D = Inf,
  use_numeric = FALSE,
  init = 5,
  tol = 1e-08,
  max_iter = 10000,
  ...
)
```

## Arguments

- object:

  A `primarycensored` object as created by
  [`new_pcens()`](https://primarycensored.epinowcast.org/dev/reference/new_pcens.md).

- p:

  A vector of probabilities at which to compute the quantiles.

- pwindow:

  Primary event window

- L:

  Minimum delay (lower truncation point). Defaults to `-Inf`, meaning no
  left truncation. For any finite value of L the distribution is
  left-truncated at L.

- D:

  Maximum delay (upper truncation point). If finite, the distribution is
  truncated at D. If set to Inf, no upper truncation is applied.
  Defaults to Inf.

- use_numeric:

  Logical; if TRUE forces the use of numeric inversion even if an
  analytical solution is available (not yet implemented).

- init:

  Half-width of the initial search interval used when one or both
  truncation bounds are infinite. The starting interval is taken as
  `[-init, init]` (or `[L, L + 2 * init]` / `[D - 2 * init, D]` when
  only one bound is finite) and then extended outward by
  [`stats::uniroot()`](https://rdrr.io/r/stats/uniroot.html) as needed.
  Defaults to 5, which brackets the bulk of most commonly used delay
  distributions.

- tol:

  Numeric tolerance passed to
  [`stats::uniroot()`](https://rdrr.io/r/stats/uniroot.html).

- max_iter:

  Maximum number of
  [`stats::uniroot()`](https://rdrr.io/r/stats/uniroot.html) iterations.

- ...:

  Additional arguments passed to underlying functions.

## Value

A numeric vector containing the computed primary event censored
quantiles.

## See also

Low level primary event censored distribution objects and methods
[`new_pcens()`](https://primarycensored.epinowcast.org/dev/reference/new_pcens.md),
[`pcens_cdf()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.md),
[`pcens_cdf.default()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.default.md),
[`pcens_cdf.pcens_pgamma_dunif()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_pgamma_dunif.md),
[`pcens_cdf.pcens_plnorm_dunif()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_plnorm_dunif.md),
[`pcens_cdf.pcens_pweibull_dunif()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_pweibull_dunif.md),
[`pcens_quantile()`](https://primarycensored.epinowcast.org/dev/reference/pcens_quantile.md)

## Examples

``` r
# Create a primarycensored object with gamma delay and uniform primary
pcens_obj <- new_pcens(
  pdist = pgamma,
  dprimary = dunif,
  dprimary_args = list(min = 0, max = 1),
  shape = 3,
  scale = 2
)

# Compute quantile for a single probability
pcens_quantile(pcens_obj, p = 0.8, pwindow = 1)
#> [1] 9.069147

# Compute quantiles for multiple probabilities
pcens_quantile(pcens_obj, p = c(0.25, 0.5, 0.75), pwindow = 1)
#> [1] 3.951257 5.853359 8.351001

# Compute quantiles for multiple probabilities with truncation
pcens_quantile(pcens_obj, p = c(0.25, 0.5, 0.75), pwindow = 1, D = 10)
#> [1] 3.666380 5.270768 7.092986

# Compute quantiles with left truncation
pcens_quantile(pcens_obj, p = c(0.25, 0.5, 0.75), pwindow = 1, L = 1, D = 10)
#> [1] 3.689124 5.285634 7.102640
```
