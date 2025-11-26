# Default method for computing primary event censored quantiles

This method inverts the primary event censored CDF using numerical
optimisation via optim. For each probability value, it searches for the
delay such that the CDF computed by
[`pcens_cdf()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.md)
approximates the target probability.

## Usage

``` r
# Default S3 method
pcens_quantile(
  object,
  p,
  pwindow,
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

- D:

  Maximum delay (truncation point). If finite, the distribution is
  truncated at D. If set to Inf, no truncation is applied. Defaults to
  Inf.

- use_numeric:

  Logical; if TRUE forces the use of numeric inversion even if an
  analytical solution is available (not yet implemented).

- init:

  Initial guess for the delay. By default, 5.

- tol:

  Numeric tolerance for the convergence criterion in the optimisation
  routine.

- max_iter:

  Integer specifying the maximum number of iterations allowed during
  optimisation.

- ...:

  Additional arguments passed to underlying functions.

## Value

A numeric vector containing the computed primary event censored
quantiles.

## Details

The quantile is computed by minimising the squared difference between
the computed CDF and the target probability.

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
#> [1] 3.666380 5.270768 7.092987
```
