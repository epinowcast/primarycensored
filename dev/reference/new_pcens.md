# S3 class for primary event censored distribution computation

S3 class for primary event censored distribution computation

## Usage

``` r
new_pcens(
  pdist,
  dprimary,
  primary_args = NULL,
  pprimary = NULL,
  dprimary_args = NULL,
  ...
)
```

## Arguments

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

  List of additional arguments to be passed to `dprimary` (and the
  looked-up `pprimary`). Replaces the deprecated `dprimary_args`.

- pprimary:

  CDF of the primary event distribution. May be a function or a
  character string naming a primary distribution in
  `pcd_primary_distributions`. When `NULL` (the default), it is looked
  up automatically from the registry using the `"name"` attribute of
  `dprimary`. When both `dprimary` and `pprimary` carry a name, the two
  must agree on everything other than the leading `d`/`p` prefix;
  mismatches such as `dunif` + `pexpgrowth` raise an error.

- dprimary_args:

  \[Deprecated\] Use `primary_args` instead.

- ...:

  Additional arguments to be passed to pdist

## Value

An object with class hierarchy
`c("pcens_{pdist_name}_{dprimary_name}", "pcens_{pdist_name}", "pcens")`.
It is a list with the fields:

- `pdist`:

  The delay distribution CDF.

- `dprimary`:

  The primary event distribution density.

- `primary_args`:

  A list of arguments passed to `dprimary` and `pprimary`.

- `dprimary_args`:

  A copy of `primary_args`, kept for backward compatibility.

- `pprimary`:

  The primary event CDF, or `NULL` if none is available.

- `args`:

  A named list of the delay distribution parameters passed through
  `...`.

It can be used with
[`pcens_cdf()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.md)
to compute the primary event censored CDF and with
[`pcens_pmf()`](https://primarycensored.epinowcast.org/dev/reference/pcens_pmf.md)
to compute the PMF. Use
[update()](https://primarycensored.epinowcast.org/dev/reference/update.pcens.md)
to change its parameters.

## See also

Low level primary event censored distribution objects and methods
[`pcens_cdf()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.md),
[`pcens_cdf.default()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.default.md),
[`pcens_cdf.pcens_pdiscretehazard()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_pdiscretehazard.md),
[`pcens_cdf.pcens_pdiscretestep()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_pdiscretestep.md),
[`pcens_cdf.pcens_pgamma_dunif()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_pgamma_dunif.md),
[`pcens_cdf.pcens_pgengamma.orig_dunif()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_pgengamma.orig_dunif.md),
[`pcens_cdf.pcens_pgengamma_dunif()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_pgengamma_dunif.md),
[`pcens_cdf.pcens_plnorm_dunif()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_plnorm_dunif.md),
[`pcens_cdf.pcens_pweibull_dunif()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_pweibull_dunif.md),
[`pcens_pmf()`](https://primarycensored.epinowcast.org/dev/reference/pcens_pmf.md),
[`pcens_pmf.default()`](https://primarycensored.epinowcast.org/dev/reference/pcens_pmf.default.md),
[`pcens_quantile()`](https://primarycensored.epinowcast.org/dev/reference/pcens_quantile.md),
[`pcens_quantile.default()`](https://primarycensored.epinowcast.org/dev/reference/pcens_quantile.default.md),
[`update.pcens()`](https://primarycensored.epinowcast.org/dev/reference/update.pcens.md)

## Examples

``` r
new_pcens(
  pdist = pgamma, dprimary = dunif,
  primary_args = list(min = 0, max = 1),
  shape = 1, scale = 1
)
#> $pdist
#> function (q, shape, rate = 1, scale = 1/rate, lower.tail = TRUE, 
#>     log.p = FALSE) 
#> {
#>     if (!missing(rate) && !missing(scale)) {
#>         if (abs(rate * scale - 1) < 1e-15) 
#>             warning("specify 'rate' or 'scale' but not both")
#>         else stop("specify 'rate' or 'scale' but not both")
#>     }
#>     .Call(C_pgamma, q, shape, scale, lower.tail, log.p)
#> }
#> <bytecode: 0x564e667a0a58>
#> <environment: namespace:stats>
#> 
#> $dprimary
#> function (x, min = 0, max = 1, log = FALSE) 
#> .Call(C_dunif, x, min, max, log)
#> <bytecode: 0x564e6476ed50>
#> <environment: namespace:stats>
#> 
#> $primary_args
#> $primary_args$min
#> [1] 0
#> 
#> $primary_args$max
#> [1] 1
#> 
#> 
#> $dprimary_args
#> $dprimary_args$min
#> [1] 0
#> 
#> $dprimary_args$max
#> [1] 1
#> 
#> 
#> $pprimary
#> function (q, min = 0, max = 1, lower.tail = TRUE, log.p = FALSE) 
#> .Call(C_punif, q, min, max, lower.tail, log.p)
#> <bytecode: 0x564e66d05a80>
#> <environment: namespace:stats>
#> 
#> $args
#> $args$shape
#> [1] 1
#> 
#> $args$scale
#> [1] 1
#> 
#> 
#> attr(,"class")
#> [1] "pcens_pgamma_dunif" "pcens_pgamma"       "pcens"             
```
