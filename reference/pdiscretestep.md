# Step (piecewise-constant) CDF

Returns the CDF of a discrete distribution whose mass is concentrated at
the right boundary of each bin. The CDF is right-continuous and
piecewise constant with jumps at `boundaries[2], ..., boundaries[K+1]`.

## Usage

``` r
pdiscretestep(q, boundaries = NULL, pmf)
```

## Arguments

- q:

  Numeric vector of quantiles.

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

Numeric vector of CDF values, the same length as `q`.

## Details

Below `boundaries[1]` the function returns 0. At `boundaries[i+1]` (the
right edge of bin *i*), F jumps by `pmf[i]`, so \\F(boundaries\[i+1\]) =
\sum\_{j=1}^{i} pmf_j\\. For `q` in \\\[boundaries\[i\],
boundaries\[i+1\])\\, F equals \\\sum\_{j=1}^{i-1} pmf_j\\. At or above
`boundaries[K+1]` the function returns 1.

### Use with [`fitdistdoublecens()`](https://primarycensored.epinowcast.org/reference/fitdistdoublecens.md)

This function carries the attribute `vector_param = "pmf"` so that
[`fitdistdoublecens`](https://primarycensored.epinowcast.org/reference/fitdistdoublecens.md)
can drive it from a flat list of scalar parameters `p1, ..., p_{K-1}`.
The free parameters are the first \\K-1\\ bin probabilities; the last is
set to `1 - sum(p1, ..., p_{K-1})`. When the implied probabilities
violate the simplex (any negative entry, or sum departing from 1 by more
than \\10^{-8}\\), the function returns 0 (or near-zero density in
[`ddiscretestep`](https://primarycensored.epinowcast.org/reference/ddiscretestep.md))
rather than erroring; this drives the optimiser back to the feasible
region.

## See also

Other pdiscretestep:
[`ddiscretestep()`](https://primarycensored.epinowcast.org/reference/ddiscretestep.md),
[`hazards_to_pmf()`](https://primarycensored.epinowcast.org/reference/hazards_to_pmf.md),
[`pmf_to_hazards()`](https://primarycensored.epinowcast.org/reference/pmf_to_hazards.md),
[`rdiscretestep()`](https://primarycensored.epinowcast.org/reference/rdiscretestep.md)

## Examples

``` r
# Two-bin PMF: mass 0.3 at x=1, mass 0.7 at x=2
pdiscretestep(c(0.5, 1, 1.5, 2), boundaries = 0:2, pmf = c(0.3, 0.7))
#> [1] 0.0 0.3 0.3 1.0
```
