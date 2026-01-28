# Compute the primary event censored CDF for delays

This function computes the primary event censored cumulative
distribution function (CDF) for a given set of quantiles. It adjusts the
CDF of the primary event distribution by accounting for the delay
distribution and potential truncation at a maximum delay (D) and minimum
delay (L). The function allows for custom primary event distributions
and delay distributions.

## Usage

``` r
pprimarycensored(
  q,
  pdist,
  pwindow = 1,
  L = 0,
  D = Inf,
  dprimary = stats::dunif,
  dprimary_args = list(),
  ...
)

ppcens(
  q,
  pdist,
  pwindow = 1,
  L = 0,
  D = Inf,
  dprimary = stats::dunif,
  dprimary_args = list(),
  ...
)
```

## Arguments

- q:

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

- ...:

  Additional arguments to be passed to pdist

## Value

Vector of primary event censored CDFs, normalized over \[L, D\] if
truncation is applied

## Details

The primary event censored CDF is computed by integrating the product of
the delay distribution function (CDF) and the primary event distribution
function (PDF) over the primary event window. The integration is
adjusted for truncation if specified.

The primary event censored CDF, \\F\_{\text{cens}}(q)\\, is given by:
\$\$ F\_{\text{cens}}(q) = \int\_{0}^{pwindow} F(q - p) \cdot
f\_{\text{primary}}(p) \\ dp \$\$ where \\F\\ is the CDF of the delay
distribution, \\f\_{\text{primary}}\\ is the PDF of the primary event
times, and \\pwindow\\ is the primary event window.

If truncation is applied (finite D or L \> 0), the CDF is normalized:
\$\$ F\_{\text{cens,norm}}(q) = \frac{F\_{\text{cens}}(q) -
F\_{\text{cens}}(L)}{ F\_{\text{cens}}(D) - F\_{\text{cens}}(L)} \$\$
where \\F\_{\text{cens,norm}}(q)\\ is the normalized CDF. For values \\q
\leq L\\, the function returns 0; for values \\q \geq D\\, it returns 1.

This function creates a `primarycensored` object using
[`new_pcens()`](https://primarycensored.epinowcast.org/dev/reference/new_pcens.md)
and then computes the primary event censored CDF using
[`pcens_cdf()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.md).
This abstraction allows for automatic use of analytical solutions when
available, while seamlessly falling back to numerical integration when
necessary.

See `methods(pcens_cdf)` for which combinations have analytical
solutions implemented.

## See also

[`new_pcens()`](https://primarycensored.epinowcast.org/dev/reference/new_pcens.md)
and
[`pcens_cdf()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.md)

Primary event censored distribution functions
[`dprimarycensored()`](https://primarycensored.epinowcast.org/dev/reference/dprimarycensored.md),
[`qprimarycensored()`](https://primarycensored.epinowcast.org/dev/reference/qprimarycensored.md),
[`rprimarycensored()`](https://primarycensored.epinowcast.org/dev/reference/rprimarycensored.md)

## Examples

``` r
# Example: Lognormal distribution with uniform primary events
pprimarycensored(c(0.1, 0.5, 1), plnorm, meanlog = 0, sdlog = 1)
#> [1] 0.0002753888 0.0475094632 0.2384217081

# Example: Lognormal distribution with exponential growth primary events
pprimarycensored(
  c(0.1, 0.5, 1), plnorm,
  dprimary = dexpgrowth,
  dprimary_args = list(r = 0.2), meanlog = 0, sdlog = 1
)
#> [1] 0.0002496934 0.0440815583 0.2290795695

# Example: Left-truncated distribution (e.g., for generation intervals)
pprimarycensored(
  c(1, 2, 3), plnorm,
  L = 1, D = 10,
  meanlog = 0, sdlog = 1
)
#> [1] 0.0000000 0.5461907 0.7719056
```
