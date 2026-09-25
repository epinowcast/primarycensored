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
  ...,
  check = TRUE
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
  ...,
  check = TRUE
)
```

## Arguments

- x:

  Vector of quantiles

- pdist:

  Distribution function (CDF). The package can identify base R
  distributions for potential analytical solutions. For non-base R
  functions, users can apply
  [`add_name_attribute()`](https://primarycensored.epinowcast.org/reference/add_name_attribute.md)
  to yield properly tagged functions if they wish to leverage the
  analytical solutions.

- pwindow:

  Primary event window. Use `pwindow = 0` for an exactly observed
  primary event, in which case the delay CDF is used directly.

- swindow:

  Secondary event window (default: 1). Use `swindow = 0` for an exactly
  observed secondary event, in which case a density is returned rather
  than a probability (see Details).

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
  [`pcd_primary_distributions()`](https://primarycensored.epinowcast.org/reference/pcd_primary_distributions.md)
  for examples. The package can identify base R distributions for
  potential analytical solutions. For non-base R functions, users can
  apply
  [`add_name_attribute()`](https://primarycensored.epinowcast.org/reference/add_name_attribute.md)
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

- check:

  Logical; if `TRUE` (the default) `pdist` is validated with
  [`check_pdist()`](https://primarycensored.epinowcast.org/reference/check_pdist.md)
  and `dprimary` with
  [`check_dprimary()`](https://primarycensored.epinowcast.org/reference/check_dprimary.md).
  Set to `FALSE` to skip both when they have already been validated.
  [`check_pdist()`](https://primarycensored.epinowcast.org/reference/check_pdist.md)
  evaluates `pdist` at four points drawn with
  [`stats::runif()`](https://rdrr.io/r/stats/Uniform.html), so skipping
  it avoids that cost and leaves the random number stream untouched.
  Must be given by its full name, as it follows `...`.

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

The function creates a `pcens` object with
[`new_pcens()`](https://primarycensored.epinowcast.org/reference/new_pcens.md)
and computes the PMF with
[`pcens_pmf()`](https://primarycensored.epinowcast.org/reference/pcens_pmf.md).
This evaluates the CDF once for all unique points (including both \\d\\
and \\d + \text{swindow}\\) and reuses these values to calculate the PMF
for each input value.

When the secondary censoring interval extends past the upper truncation
point (\\d + \text{swindow} \> D\\) but the lower endpoint satisfies \\d
\< D\\, the upper endpoint is internally clipped to \\D\\ before
evaluating the CDF. The likelihood for such an observation is \\P(X \in
\[d, \min(d + \text{swindow}, D)\] \mid L \le X \le D)\\, which equals
the usual interval probability when \\d + \text{swindow} \le D\\. This
avoids erroring when an observation's secondary window straddles the
truncation point (relevant for non-parametric delays such as
[`pdiscretestep()`](https://primarycensored.epinowcast.org/reference/pdiscretestep.md)).

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
[`pprimarycensored()`](https://primarycensored.epinowcast.org/reference/pprimarycensored.md).

### Zero-width windows

With `pwindow = 0` the primary event time is known exactly and the
primary event censored CDF is the delay CDF, so the PMF is \\F(d +
\text{swindow}) - F(d)\\.

With `swindow = 0` the secondary event time is known exactly. The
probability of the interval is then zero, so the density of the primary
event censored delay at \\d\\ is returned instead. This is the
derivative of \\F\_{\text{cens}}\\ at \\d\\, the limit of the PMF
divided by `swindow` as `swindow` goes to zero. With `pwindow = 0` as
well it is the delay density. With a uniform primary event distribution
it is \\(F(d) - F(d - \text{pwindow})) / \text{pwindow}\\. Otherwise the
delay density is integrated against the primary event density. The delay
density is found from the name of `pdist` (for example
[`dgamma()`](https://rdrr.io/r/stats/GammaDist.html) for
[`pgamma()`](https://rdrr.io/r/stats/GammaDist.html)) and an error is
raised if it cannot be found. Densities are normalised for truncation in
the same way as probabilities. `swindow` may be a vector, so densities
and probabilities can be mixed in one call.

## See also

Primary event censored distribution functions
[`pprimarycensored()`](https://primarycensored.epinowcast.org/reference/pprimarycensored.md),
[`qprimarycensored()`](https://primarycensored.epinowcast.org/reference/qprimarycensored.md),
[`rprimarycensored()`](https://primarycensored.epinowcast.org/reference/rprimarycensored.md)

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

# Example: exact primary events, and exact secondary events (a density)
dprimarycensored(1:3, pweibull, pwindow = 0, shape = 1.5, scale = 2.0)
#> [1] 0.3343091 0.2086035 0.1001702
dprimarycensored(
  1:3, pweibull,
  pwindow = 1, swindow = 0, shape = 1.5, scale = 2.0
)
#> [1] 0.2978115 0.3343091 0.2086035
```
