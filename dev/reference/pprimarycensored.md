# Compute the primary event censored CDF for delays

This function computes the primary event censored cumulative
distribution function (CDF) for a given set of quantiles. It adjusts the
CDF of the primary event distribution by accounting for the delay
distribution and potential truncation at a maximum delay (D). The
function allows for custom primary event distributions and delay
distributions.

## Usage

``` r
pprimarycensored(
  q,
  pdist,
  pwindow = 1,
  D = Inf,
  dprimary = stats::dunif,
  dprimary_args = list(),
  pdist_name = lifecycle::deprecated(),
  dprimary_name = lifecycle::deprecated(),
  ...
)

ppcens(
  q,
  pdist,
  pwindow = 1,
  D = Inf,
  dprimary = stats::dunif,
  dprimary_args = list(),
  pdist_name = lifecycle::deprecated(),
  dprimary_name = lifecycle::deprecated(),
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

  Additional arguments to be passed to pdist

## Value

Vector of primary event censored CDFs, normalized by D if finite
(truncation adjustment)

## Details

The primary event censored CDF is computed by integrating the product of
the delay distribution function (CDF) and the primary event distribution
function (PDF) over the primary event window. The integration is
adjusted for truncation if a finite maximum delay (D) is specified.

The primary event censored CDF, \\F\_{\text{cens}}(q)\\, is given by:
\$\$ F\_{\text{cens}}(q) = \int\_{0}^{pwindow} F(q - p) \cdot
f\_{\text{primary}}(p) \\ dp \$\$ where \\F\\ is the CDF of the delay
distribution, \\f\_{\text{primary}}\\ is the PDF of the primary event
times, and \\pwindow\\ is the primary event window.

If the maximum delay \\D\\ is finite, the CDF is normalized by dividing
by \\F\_{\text{cens}}(D)\\: \$\$ F\_{\text{cens,norm}}(q) =
\frac{F\_{\text{cens}}(q)}{F\_{\text{cens}}(D)} \$\$ where
\\F\_{\text{cens,norm}}(q)\\ is the normalized CDF.

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
```
