# Compute the primary event censored PMF for delays

This function computes the primary event censored probability mass
function (PMF) for a given set of quantiles. It adjusts the PMF of the
primary event distribution by accounting for the delay distribution and
potential truncation at a maximum delay (D). The function allows for
custom primary event distributions and delay distributions.

## Usage

``` r
dprimarycensored(
  x,
  pdist,
  pwindow = 1,
  swindow = 1,
  D = Inf,
  dprimary = stats::dunif,
  dprimary_args = list(),
  log = FALSE,
  pdist_name = lifecycle::deprecated(),
  dprimary_name = lifecycle::deprecated(),
  ...
)

dpcens(
  x,
  pdist,
  pwindow = 1,
  swindow = 1,
  D = Inf,
  dprimary = stats::dunif,
  dprimary_args = list(),
  log = FALSE,
  pdist_name = lifecycle::deprecated(),
  dprimary_name = lifecycle::deprecated(),
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

- D:

  Maximum delay (truncation point). If finite, the distribution is
  truncated at D. If set to Inf, no truncation is applied. Defaults to
  Inf.

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

- pdist_name:

  **\[deprecated\]** this argument will be ignored in future versions;
  use
  [`add_name_attribute()`](https://primarycensored.epinowcast.org/dev/reference/add_name_attribute.md)
  on `pdist` instead

- dprimary_name:

  **\[deprecated\]** this argument will be ignored in future versions;
  use
  [`add_name_attribute()`](https://primarycensored.epinowcast.org/dev/reference/add_name_attribute.md)
  on `dprimary` instead

- ...:

  Additional arguments to be passed to the distribution function

## Value

Vector of primary event censored PMFs, normalized by D if finite
(truncation adjustment)

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
the PMF for each input value. For non-positive delays, the function
returns 0.

If a finite maximum delay \\D\\ is specified, the PMF is normalized to
ensure it sums to 1 over the range \[0, D\]. This normalization can be
expressed as: \$\$ f\_{\text{cens,norm}}(d) =
\frac{f\_{\text{cens}}(d)}{\sum\_{i=0}^{D-1} f\_{\text{cens}}(i)} \$\$
where \\f\_{\text{cens,norm}}(d)\\ is the normalized PMF and
\\f\_{\text{cens}}(d)\\ is the unnormalized PMF. For the explanation and
mathematical details of the CDF, refer to the documentation of
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
```
