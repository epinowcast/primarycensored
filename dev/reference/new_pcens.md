# S3 class for primary event censored distribution computation

S3 class for primary event censored distribution computation

## Usage

``` r
new_pcens(
  pdist,
  dprimary,
  dprimary_args,
  pdist_name = lifecycle::deprecated(),
  dprimary_name = lifecycle::deprecated(),
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

An object of class `pcens_{pdist_name}_{dprimary_name}`. This contains
the primary event distribution, the delay distribution, the delay
distribution arguments, and any additional arguments. It can be used
with the
[`pcens_cdf()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.md)
function to compute the primary event censored cdf.

## See also

Low level primary event censored distribution objects and methods
[`pcens_cdf()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.md),
[`pcens_cdf.default()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.default.md),
[`pcens_cdf.pcens_pgamma_dunif()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_pgamma_dunif.md),
[`pcens_cdf.pcens_plnorm_dunif()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_plnorm_dunif.md),
[`pcens_cdf.pcens_pweibull_dunif()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_pweibull_dunif.md),
[`pcens_quantile()`](https://primarycensored.epinowcast.org/dev/reference/pcens_quantile.md),
[`pcens_quantile.default()`](https://primarycensored.epinowcast.org/dev/reference/pcens_quantile.default.md)

## Examples

``` r
new_pcens(
  pdist = pgamma, dprimary = dunif, dprimary_args = list(min = 0, max = 1),
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
#> <bytecode: 0x555f242f6dc0>
#> <environment: namespace:stats>
#> 
#> $dprimary
#> function (x, min = 0, max = 1, log = FALSE) 
#> .Call(C_dunif, x, min, max, log)
#> <bytecode: 0x555f2419c560>
#> <environment: namespace:stats>
#> 
#> $dprimary_args
#> $dprimary_args$min
#> [1] 0
#> 
#> $dprimary_args$max
#> [1] 1
#> 
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
#> [1] "pcens_pgamma_dunif"
```
