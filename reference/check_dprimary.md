# Check if a function is a valid bounded probability density function (PDF)

This function tests whether a given function behaves like a valid PDF by
checking if it integrates to approximately 1 over the specified range
and if it takes the arguments min and max.

## Usage

``` r
check_dprimary(dprimary, pwindow, dprimary_args = list(), tolerance = 0.001)
```

## Arguments

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

- pwindow:

  Primary event window

- dprimary_args:

  List of additional arguments to be passed to dprimary. For example,
  when using `dexpgrowth`, you would pass
  `list(min = 0, max = pwindow, r = 0.2)` to set the minimum, maximum,
  and rate parameters

- tolerance:

  The tolerance for the integral to be considered close to 1

## Value

NULL. The function will stop execution with an error message if dprimary
is not a valid PDF.

## See also

Distribution checking functions
[`check_pdist()`](https://primarycensored.epinowcast.org/reference/check_pdist.md),
[`check_truncation()`](https://primarycensored.epinowcast.org/reference/check_truncation.md)

## Examples

``` r
check_dprimary(dunif, pwindow = 1)
```
