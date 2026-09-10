# Sample from a step distribution

Draws `n` independent samples from the discrete distribution with
probability mass `pmf[i]` at `boundaries[i+1]`.

## Usage

``` r
rdiscretestep(n, boundaries = NULL, pmf)
```

## Arguments

- n:

  Integer. Number of samples to draw.

- boundaries:

  Numeric vector of length \\K+1\\ defining the bin edges. Must be
  strictly increasing. Defaults to `0:K` (unit-width daily bins) where
  `K` is inferred from `length(pmf)`.

- pmf:

  Numeric vector of length \\K\\ giving the probability mass for each
  bin. Must be non-negative and sum to approximately 1; if either
  condition is violated the function returns a vector of zeros (a soft
  simplex penalty for use inside optimisation).

## Value

Numeric vector of length `n`.

## See also

Other pdiscretestep:
[`ddiscretestep()`](https://primarycensored.epinowcast.org/dev/reference/ddiscretestep.md),
[`hazards_to_pmf()`](https://primarycensored.epinowcast.org/dev/reference/hazards_to_pmf.md),
[`pdiscretestep()`](https://primarycensored.epinowcast.org/dev/reference/pdiscretestep.md),
[`pmf_to_hazards()`](https://primarycensored.epinowcast.org/dev/reference/pmf_to_hazards.md)

## Examples

``` r
set.seed(42)
rdiscretestep(10, boundaries = 0:3, pmf = c(0.2, 0.5, 0.3))
#>  [1] 1 1 2 1 3 3 3 2 3 3
```
