# Step (piecewise-constant) PMF

Returns the probability mass for each value in `x`. Mass `pmf[i]` is
located at `boundaries[i+1]` (the right edge of bin *i*); all other
values return 0.

## Usage

``` r
ddiscretestep(x, boundaries = NULL, pmf)
```

## Arguments

- x:

  Numeric vector of values at which to evaluate the PMF.

- boundaries:

  Numeric vector of length \\K+1\\ defining the bin edges. Must be
  strictly increasing. Defaults to `0:K` (unit-width daily bins) where
  `K` is inferred from `length(pmf)`.

- pmf:

  Numeric vector of length \\K\\ giving the probability mass for each
  bin. Must be non-negative and sum to approximately 1; if either
  condition is violated the function returns a vector of
  `.Machine$double.eps` rather than 0 (a soft simplex penalty that keeps
  log-density finite inside fitting closures).

## Value

Numeric vector of PMF values, the same length as `x`.

## Details

Like
[`pdiscretestep`](https://primarycensored.epinowcast.org/dev/reference/pdiscretestep.md),
this function applies a soft simplex penalty: if `pmf` contains negative
entries or fails to sum to 1 (within \\10^{-8}\\), the function returns
near-zero density (`.Machine$double.eps`) rather than erroring. This
makes it safe to call from inside
[`fitdistrplus::fitdist()`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.html)
closures driven by
[`fitdistdoublecens`](https://primarycensored.epinowcast.org/dev/reference/fitdistdoublecens.md).

## See also

Other pdiscretestep:
[`hazards_to_pmf()`](https://primarycensored.epinowcast.org/dev/reference/hazards_to_pmf.md),
[`pdiscretestep()`](https://primarycensored.epinowcast.org/dev/reference/pdiscretestep.md),
[`pmf_to_hazards()`](https://primarycensored.epinowcast.org/dev/reference/pmf_to_hazards.md),
[`rdiscretestep()`](https://primarycensored.epinowcast.org/dev/reference/rdiscretestep.md)

## Examples

``` r
ddiscretestep(c(0, 1, 2, 3), boundaries = 0:3, pmf = c(0.2, 0.5, 0.3))
#> [1] 0.0 0.2 0.5 0.3
```
