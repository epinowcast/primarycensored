# Compute quantiles corresponding to target probabilities for primary event censored delays

This function computes the quantiles (delay values) that correspond to
specified probabilities in the primary event censored distribution. For
a given probability p, it computes the delay value q such that the
cumulative probability up to q equals p in the primary event censored
distribution. The distribution accounts for both the delay distribution
and the primary event timing distribution.

## Usage

``` r
qprimarycensored(
  p,
  pdist,
  pwindow = 1,
  L = -Inf,
  D = Inf,
  dprimary = dunif,
  primary_args = NULL,
  pprimary = NULL,
  dprimary_args = NULL,
  ...
)

qpcens(
  p,
  pdist,
  pwindow = 1,
  L = -Inf,
  D = Inf,
  dprimary = dunif,
  primary_args = NULL,
  pprimary = NULL,
  dprimary_args = NULL,
  ...
)
```

## Arguments

- p:

  Vector of probabilities between 0 and 1 for which to compute
  corresponding quantiles

- pdist:

  Distribution function (CDF). The package can identify base R
  distributions for potential analytical solutions. For non-base R
  functions, users can apply
  [`add_name_attribute()`](https://primarycensored.epinowcast.org/dev/reference/add_name_attribute.md)
  to yield properly tagged functions if they wish to leverage the
  analytical solutions.

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

- ...:

  Additional arguments to be passed to pdist

## Value

Vector of delay values (quantiles) corresponding to the input
probabilities

## Details

For each probability, the function finds the delay value where that
proportion of events have occurred by that time in the primary event
censored distribution. This is done by inverting the cumulative
distribution function.

The function creates a `primarycensored` object using
[`new_pcens()`](https://primarycensored.epinowcast.org/dev/reference/new_pcens.md)
and then computes the quantiles using
[`pcens_quantile()`](https://primarycensored.epinowcast.org/dev/reference/pcens_quantile.md).
This approach allows for analytical solutions when available, falling
back to numerical methods when necessary.

For example, if p = 0.5, the function returns the median delay
(truncated over \[L, D\] if specified) where 50% of censored events
occur by this time and 50% occur after.

See `methods(pcens_quantile)` for which combinations have analytical
solutions implemented.

## See also

[`new_pcens()`](https://primarycensored.epinowcast.org/dev/reference/new_pcens.md)
and
[`pcens_quantile()`](https://primarycensored.epinowcast.org/dev/reference/pcens_quantile.md)

Primary event censored distribution functions
[`dprimarycensored()`](https://primarycensored.epinowcast.org/dev/reference/dprimarycensored.md),
[`pprimarycensored()`](https://primarycensored.epinowcast.org/dev/reference/pprimarycensored.md),
[`rprimarycensored()`](https://primarycensored.epinowcast.org/dev/reference/rprimarycensored.md)

## Examples

``` r
# Compute delays where 25%, 50%, and 75% of events occur by (quartiles)
# Using lognormal delays with uniform primary events
qprimarycensored(c(0.25, 0.5, 0.75), plnorm, meanlog = 0, sdlog = 1)
#> [1] 1.022949 1.540771 2.498358

# Same quartiles but with exponential growth in primary events
qprimarycensored(
  c(0.25, 0.5, 0.75), plnorm,
  dprimary = dexpgrowth,
  primary_args = list(r = 0.2), meanlog = 0, sdlog = 1
)
#> [1] 1.041285 1.557111 2.514701

# Same quartiles but with truncation at 10
qprimarycensored(
  c(0.25, 0.5, 0.75), plnorm,
  dprimary = dexpgrowth,
  primary_args = list(r = 0.2), meanlog = 0, sdlog = 1, D = 10
)
#> [1] 1.035312 1.541788 2.459511

# Left-truncated distribution (e.g., for generation intervals)
qprimarycensored(
  c(0.25, 0.5, 0.75), plnorm,
  L = 1, D = 10, meanlog = 0, sdlog = 1
)
#> [1] 1.368467 1.872762 2.856596
```
