# Convert discrete-time hazards to a PMF

Given a vector of discrete-time hazards \\h_1, \ldots, h_K\\, computes
the corresponding PMF: \$\$pmf_i = h_i \prod\_{j \< i} (1 - h_j)\$\$

## Usage

``` r
hazards_to_pmf(hazards)
```

## Arguments

- hazards:

  Numeric vector of hazards in \\\[0, 1\]\\. The last entry must equal 1
  (within \\10^{-8}\\), or a vector of length \\K-1\\ may be supplied
  and the trailing 1 will be appended.

## Value

Numeric vector of PMF values.

## Details

The last hazard must equal 1 (to ensure the PMF sums to 1). If a vector
of length \\K-1\\ is supplied (all hazards except the final exit
hazard), 1 is appended automatically.

## See also

Other pdiscretestep:
[`ddiscretestep()`](https://primarycensored.epinowcast.org/reference/ddiscretestep.md),
[`pdiscretestep()`](https://primarycensored.epinowcast.org/reference/pdiscretestep.md),
[`pmf_to_hazards()`](https://primarycensored.epinowcast.org/reference/pmf_to_hazards.md),
[`rdiscretestep()`](https://primarycensored.epinowcast.org/reference/rdiscretestep.md)

## Examples

``` r
hazards_to_pmf(c(0.2, 0.3, 1))
#> [1] 0.20 0.24 0.56
```
