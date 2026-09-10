# Convert a PMF to discrete-time hazards

Inverts
[`hazards_to_pmf()`](https://primarycensored.epinowcast.org/reference/hazards_to_pmf.md):
given a PMF, computes the discrete-time conditional hazard at each time
point. \$\$h_i = pmf_i / (1 - \sum\_{j \< i} pmf_j)\$\$

## Usage

``` r
pmf_to_hazards(pmf)
```

## Arguments

- pmf:

  Numeric vector of probabilities. Must be non-negative and sum to
  approximately 1.

## Value

Numeric vector of hazards in \\\[0, 1\]\\.

## Details

The returned vector has the same length as `pmf`, with the last entry
equal to 1.

## See also

Other pdiscretestep:
[`ddiscretestep()`](https://primarycensored.epinowcast.org/reference/ddiscretestep.md),
[`hazards_to_pmf()`](https://primarycensored.epinowcast.org/reference/hazards_to_pmf.md),
[`pdiscretestep()`](https://primarycensored.epinowcast.org/reference/pdiscretestep.md),
[`rdiscretestep()`](https://primarycensored.epinowcast.org/reference/rdiscretestep.md)

## Examples

``` r
pmf_to_hazards(c(0.2, 0.3, 0.5))
#> [1] 0.200 0.375 1.000
```
