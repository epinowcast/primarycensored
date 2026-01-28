# Check if a function is a valid cumulative distribution function (CDF)

This function tests whether a given function behaves like a valid CDF by
checking if it's monotonically increasing and bounded between 0 and 1.

## Usage

``` r
check_pdist(pdist, D = Inf, ...)
```

## Arguments

- pdist:

  Distribution function (CDF). The package can identify base R
  distributions for potential analytical solutions. For non-base R
  functions, users can apply
  [`add_name_attribute()`](https://primarycensored.epinowcast.org/dev/reference/add_name_attribute.md)
  to yield properly tagged functions if they wish to leverage the
  analytical solutions.

- D:

  Maximum delay (upper truncation point). If finite, the distribution is
  truncated at D. If set to Inf, no upper truncation is applied.
  Defaults to Inf.

- ...:

  Additional arguments to be passed to pdist

## Value

NULL. The function will stop execution with an error message if pdist is
not a valid CDF.

## See also

Distribution checking functions
[`check_dprimary()`](https://primarycensored.epinowcast.org/dev/reference/check_dprimary.md),
[`check_truncation()`](https://primarycensored.epinowcast.org/dev/reference/check_truncation.md)

## Examples

``` r
check_pdist(pnorm, D = 10)
```
