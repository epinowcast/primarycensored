# Get the pcens object for a likelihood evaluation

Builds a `pcens` object with
[`.build_pcens()`](https://primarycensored.epinowcast.org/dev/reference/dot-build_pcens.md),
or, when `cache` already holds one, updates its delay parameters with
[update()](https://primarycensored.epinowcast.org/dev/reference/update.pcens.md).

## Usage

``` r
.fit_pcens_state(cache, pdist, dprimary, primary_args, pprimary, args)
```

## Arguments

- cache:

  Environment to keep the object in, or `NULL`.

- pdist:

  Distribution function (CDF). The package can identify base R
  distributions for potential analytical solutions. For non-base R
  functions, users can apply
  [`add_name_attribute()`](https://primarycensored.epinowcast.org/dev/reference/add_name_attribute.md)
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

- args:

  Named list of delay distribution parameters.

## Value

An environment with the `pcens` object in `obj`. This is `cache` when it
is not `NULL`.
