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
  D = Inf,
  dprimary = stats::dunif,
  dprimary_args = list(),
  ...
)

qpcens(
  p,
  pdist,
  pwindow = 1,
  D = Inf,
  dprimary = stats::dunif,
  dprimary_args = list(),
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
  [`add_name_attribute()`](https://primarycensored.epinowcast.org/reference/add_name_attribute.md)
  to yield properly tagged functions if they wish to leverage the
  analytical solutions.

- pwindow:

  Primary event window

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
  [`pcd_primary_distributions()`](https://primarycensored.epinowcast.org/reference/pcd_primary_distributions.md)
  for examples. The package can identify base R distributions for
  potential analytical solutions. For non-base R functions, users can
  apply
  [`add_name_attribute()`](https://primarycensored.epinowcast.org/reference/add_name_attribute.md)
  to yield properly tagged functions if they wish to leverage analytical
  solutions.

- dprimary_args:

  List of additional arguments to be passed to dprimary. For example,
  when using `dexpgrowth`, you would pass
  `list(min = 0, max = pwindow, r = 0.2)` to set the minimum, maximum,
  and rate parameters

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
[`new_pcens()`](https://primarycensored.epinowcast.org/reference/new_pcens.md)
and then computes the quantiles using
[`pcens_quantile()`](https://primarycensored.epinowcast.org/reference/pcens_quantile.md).
This approach allows for analytical solutions when available, falling
back to numerical methods when necessary.

For example, if p = 0.5, the function returns the median delay - the
value where 50% of censored events occur by this time and 50% occur
after.

See `methods(pcens_quantile)` for which combinations have analytical
solutions implemented.

## See also

[`new_pcens()`](https://primarycensored.epinowcast.org/reference/new_pcens.md)
and
[`pcens_quantile()`](https://primarycensored.epinowcast.org/reference/pcens_quantile.md)

Primary event censored distribution functions
[`dprimarycensored()`](https://primarycensored.epinowcast.org/reference/dprimarycensored.md),
[`pprimarycensored()`](https://primarycensored.epinowcast.org/reference/pprimarycensored.md),
[`rprimarycensored()`](https://primarycensored.epinowcast.org/reference/rprimarycensored.md)

## Examples

``` r
# Compute delays where 25%, 50%, and 75% of events occur by (quartiles)
# Using lognormal delays with uniform primary events
qprimarycensored(c(0.25, 0.5, 0.75), plnorm, meanlog = 0, sdlog = 1)
#> [1] 1.022948 1.540771 2.498358

# Same quartiles but with exponential growth in primary events
qprimarycensored(
  c(0.25, 0.5, 0.75), plnorm,
  dprimary = dexpgrowth,
  dprimary_args = list(r = 0.2), meanlog = 0, sdlog = 1
)
#> [1] 0.000000 1.557111 2.514701

# Same quartiles but with truncation at 10
qprimarycensored(
  c(0.25, 0.5, 0.75), plnorm,
  dprimary = dexpgrowth,
  dprimary_args = list(r = 0.2), meanlog = 0, sdlog = 1, D = 10
)
#> [1] 0.000000 1.541789 2.459511
```
