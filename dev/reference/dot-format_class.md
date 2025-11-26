# Extract and Combine Distribution Names

This helper function attempts to determine distribution names and uses
those to establish a class name for potential analytical solutions.

## Usage

``` r
.format_class(pdist, dprimary)
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

## Value

a character string representing the combined distribution class
