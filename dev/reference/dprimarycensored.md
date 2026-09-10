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
  L = -Inf,
  D = Inf,
  dprimary = dunif,
  primary_args = NULL,
  pprimary = NULL,
  dprimary_args = NULL,
  log = FALSE,
  ...
)

dpcens(
  x,
  pdist,
  pwindow = 1,
  swindow = 1,
  L = -Inf,
  D = Inf,
  dprimary = dunif,
  primary_args = NULL,
  pprimary = NULL,
  dprimary_args = NULL,
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

  Minimum delay (lower truncation point). Defaults to `-Inf`, meaning no
  left truncation. For any finite value of L the distribution is
  left-truncated at L.

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

- primary_args:

  List of additional arguments to be passed to dprimary (and the
  matching primary CDF). For example, when using `dexpgrowth`, you would
  pass `list(min = 0, max = pwindow, r = 0.2)` to set the minimum,
  maximum, and rate parameters. Replaces the deprecated `dprimary_args`;
  defaults to `NULL`.

- pprimary:

  Optional CDF for the primary event distribution. May be a function or
  a character string naming a primary distribution in
  `pcd_primary_distributions`. Defaults to `NULL`, in which case the
  primary CDF is looked up automatically from the registry using the
  `"name"` attribute of `dprimary`. When both `dprimary` and `pprimary`
  carry a `"name"` attribute (or are base R functions whose name can be
  inferred), the two names must agree on everything other than the
  leading `d`/`p` prefix; mismatches such as `dunif` + `pexpgrowth`
  raise an error. Supplying `pprimary` explicitly is mainly useful when
  using a custom primary distribution whose CDF is not in the registry.

- dprimary_args:

  \[Deprecated\] Use `primary_args` instead.

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

When the secondary censoring interval extends past the upper truncation
point (\\d + \text{swindow} \> D\\) but the lower endpoint satisfies \\d
\< D\\, the upper endpoint is internally clipped to \\D\\ before
evaluating the CDF. The likelihood for such an observation is \\P(X \in
\[d, \min(d + \text{swindow}, D)\] \mid L \le X \le D)\\, which equals
the usual interval probability when \\d + \text{swindow} \le D\\. This
avoids erroring when an observation's secondary window straddles the
truncation point (relevant for non-parametric delays such as
[`pdiscretestep()`](https://primarycensored.epinowcast.org/dev/reference/pdiscretestep.md)).

Observations with \\d \ge D\\ are rejected with an error: under the
truncation \\X \le D\\, no event with latent value \\d \ge D\\ is
observable, and accepting such inputs would otherwise yield a 0/0
likelihood.

The PMF is normalised to ensure it sums to 1 over the range \[L, D\\.
This normalization uses: \$\$ f\_{\text{cens,norm}}(d) =
\frac{f\_{\text{cens}}(d)}{ F\_{\text{cens}}(D) - F\_{\text{cens}}(L)}
\$\$ where \\f\_{\text{cens,norm}}(d)\\ is the normalized PMF. For the
explanation and mathematical details of the CDF, refer to the
documentation of
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
  primary_args = list(r = 0.2), shape = 1.5, scale = 2.0
)
#> [1] 0.1522796 0.2691280 0.3459055

# Example: Left-truncated distribution (e.g., for generation intervals)
dprimarycensored(1:9, pweibull, L = 1, D = 10, shape = 1.5, scale = 2.0)
#> [1] 0.3967387124 0.3138303103 0.1723520068 0.0760439783 0.0283706839
#> [6] 0.0091967620 0.0026354003 0.0006757134 0.0001564326
```
