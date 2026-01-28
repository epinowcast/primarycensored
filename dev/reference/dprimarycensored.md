# Compute the primary event censored PMF for delays

This function computes the primary event censored probability mass
function (PMF) for a given set of quantiles. It adjusts the PMF of the
primary event distribution by accounting for the delay distribution and
potential truncation at a maximum delay (D) and minimum delay (L). The
function allows for custom primary event distributions and delay
distributions.

## Usage

``` r
dprimarycensored(
  x,
  pdist,
  pwindow = 1,
  swindow = 1,
  L = 0,
  D = Inf,
  dprimary = stats::dunif,
  dprimary_args = list(),
  log = FALSE,
  ...
)

dpcens(
  x,
  pdist,
  pwindow = 1,
  swindow = 1,
  L = 0,
  D = Inf,
  dprimary = stats::dunif,
  dprimary_args = list(),
  log = FALSE,
  ...
)
```

## Arguments

- x:

  Vector of quantiles

- pdist:

  Distribution function (CDF). The package can identify base R
  distributions for potential analytical solutions. For non-base R
  functions, users can apply
  [`add_name_attribute()`](https://primarycensored.epinowcast.org/dev/reference/add_name_attribute.md)
  to yield properly tagged functions if they wish to leverage the
  analytical solutions.

- pwindow:

  Primary event window

- swindow:

  Secondary event window (default: 1)

- L:

  Minimum delay (lower truncation point). If greater than 0, the
  distribution is left-truncated at L. This is useful for modelling
  generation intervals where day 0 is excluded, particularly when used
  in renewal models. Defaults to 0 (no left truncation).

- D:

  Maximum delay (upper truncation point). If finite, the distribution is
  truncated at D. If set to Inf, no upper truncation is applied.
  Defaults to Inf.

- dprimary:

  Function to generate the probability density function (PDF) of primary
  event times. This function should take a value `x` and a `pwindow`
  parameter, and return a probability density. It should be normalized
  to integrate to 1 over \[0, pwindow\]. Defaults to a uniform
  distribution over \[0, pwindow\]. Users can provide custom functions
  or use helper functions like `dexpgrowth` for an exponential growth
  distribution. See
  [`pcd_primary_distributions()`](https://primarycensored.epinowcast.org/dev/reference/pcd_primary_distributions.md)
  for examples. The package can identify base R distributions for
  potential analytical solutions. For non-base R functions, users can
  apply
  [`add_name_attribute()`](https://primarycensored.epinowcast.org/dev/reference/add_name_attribute.md)
  to yield properly tagged functions if they wish to leverage analytical
  solutions.

- dprimary_args:

  List of additional arguments to be passed to dprimary. For example,
  when using `dexpgrowth`, you would pass
  `list(min = 0, max = pwindow, r = 0.2)` to set the minimum, maximum,
  and rate parameters

- log:

  Logical; if TRUE, probabilities p are given as log(p)

- ...:

  Additional arguments to be passed to the distribution function

## Value

Vector of primary event censored PMFs, normalized over \[L, D\] if
truncation is applied

## Details

The primary event censored PMF is computed by taking the difference of
the primary event censored cumulative distribution function (CDF) at two
points, \\d + \text{swindow}\\ and \\d\\. The primary event censored
PMF, \\f\_{\text{cens}}(d)\\, is given by: \$\$ f\_{\text{cens}}(d) =
F\_{\text{cens}}(d + \text{swindow}) - F\_{\text{cens}}(d) \$\$ where
\\F\_{\text{cens}}\\ is the primary event censored CDF.

The function first computes the CDFs for all unique points (including
both \\d\\ and \\d + \text{swindow}\\) using
[`pprimarycensored()`](https://primarycensored.epinowcast.org/dev/reference/pprimarycensored.md).
It then creates a lookup table for these CDFs to efficiently calculate
the PMF for each input value. For delays less than L, the function
returns 0.

If truncation is applied (finite D or L \> 0), the PMF is normalized to
ensure it sums to 1 over the range \[L, D\\. This normalization uses:
\$\$ f\_{\text{cens,norm}}(d) = \frac{f\_{\text{cens}}(d)}{
F\_{\text{cens}}(D) - F\_{\text{cens}}(L)} \$\$ where
\\f\_{\text{cens,norm}}(d)\\ is the normalized PMF. For the explanation
and mathematical details of the CDF, refer to the documentation of
[`pprimarycensored()`](https://primarycensored.epinowcast.org/dev/reference/pprimarycensored.md).

## See also

Primary event censored distribution functions
[`pprimarycensored()`](https://primarycensored.epinowcast.org/dev/reference/pprimarycensored.md),
[`qprimarycensored()`](https://primarycensored.epinowcast.org/dev/reference/qprimarycensored.md),
[`rprimarycensored()`](https://primarycensored.epinowcast.org/dev/reference/rprimarycensored.md)

## Examples

``` r
# Example: Weibull distribution with uniform primary events
dprimarycensored(c(0.1, 0.5, 1), pweibull, shape = 1.5, scale = 2.0)
#> [1] 0.1577965 0.2735269 0.3463199

# Example: Weibull distribution with exponential growth primary events
dprimarycensored(
  c(0.1, 0.5, 1), pweibull,
  dprimary = dexpgrowth,
  dprimary_args = list(r = 0.2), shape = 1.5, scale = 2.0
)
#> [1] 0.1522796 0.2691280 0.3459055

# Example: Left-truncated distribution (e.g., for generation intervals)
dprimarycensored(1:9, pweibull, L = 1, D = 10, shape = 1.5, scale = 2.0)
#> [1] 0.3967387124 0.3138303103 0.1723520068 0.0760439783 0.0283706839
#> [6] 0.0091967620 0.0026354003 0.0006757134 0.0001564326
```
