# Define a fitdistrplus compatible wrapper around dprimarycensored

Define a fitdistrplus compatible wrapper around dprimarycensored

## Usage

``` r
.dpcens(
  x,
  params,
  pdist,
  dprimary,
  primary_args,
  pprimary = NULL,
  check = TRUE,
  ...
)
```

## Arguments

- x:

  Vector of quantiles

- params:

  A data frame with columns 'swindow', 'pwindow', 'L', and 'D'
  corresponding to the secondary window sizes, primary window sizes,
  upper truncation times, and lower truncation times for each element in
  x.

- pdist:

  Distribution function (CDF). The package can identify base R
  distributions for potential analytical solutions. For non-base R
  functions, users can apply
  [`add_name_attribute()`](https://primarycensored.epinowcast.org/reference/add_name_attribute.md)
  to yield properly tagged functions if they wish to leverage the
  analytical solutions.

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

- ...:

  Additional arguments to be passed to the distribution function
