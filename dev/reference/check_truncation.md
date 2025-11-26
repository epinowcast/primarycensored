# Check if truncation time is appropriate relative to the maximum delay

This function checks if the truncation time D is appropriate relative to
the maximum delay. If D is much larger than necessary, it suggests
considering setting it to `Inf` for better efficiency with minimal
accuracy cost.

## Usage

``` r
check_truncation(delays, D, multiplier = 2)
```

## Arguments

- delays:

  A numeric vector of delay times

- D:

  The truncation time

- multiplier:

  The multiplier for the maximum delay to compare with D. Default is 2.

## Value

Invisible NULL. Prints a message if the condition is met.

## See also

Distribution checking functions
[`check_dprimary()`](https://primarycensored.epinowcast.org/dev/reference/check_dprimary.md),
[`check_pdist()`](https://primarycensored.epinowcast.org/dev/reference/check_pdist.md)

## Examples

``` r
check_truncation(delays = c(1, 2, 3, 4), D = 10, multiplier = 2)
#> The truncation time D (10) is larger than 2 times the maximum observed delay (4). Consider setting D to Inf for better efficiency with minimal accuracy cost for this case.
```
