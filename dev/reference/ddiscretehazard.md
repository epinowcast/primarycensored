# Hazard-parameterised piecewise-constant PMF

Returns the probability mass for each value in `x`. Converts `hazards`
to a PMF via
[`hazards_to_pmf()`](https://primarycensored.epinowcast.org/dev/reference/hazards_to_pmf.md)
then delegates to
[`ddiscretestep()`](https://primarycensored.epinowcast.org/dev/reference/ddiscretestep.md).

## Usage

``` r
ddiscretehazard(x, boundaries = NULL, hazards)
```

## Arguments

- x:

  Numeric vector of values at which to evaluate the PMF.

- boundaries:

  Numeric vector of length \\K+1\\ defining the bin edges. Must be
  strictly increasing. Defaults to `0:K`.

- hazards:

  Numeric vector of length \\K\\ (or \\K-1\\, in which case a trailing 1
  is appended) giving the discrete-time conditional hazard for each bin.
  Values must be in \\\[0, 1\]\\.

## Value

Numeric vector of PMF values, the same length as `x`.

## Details

Outside fitting it is a deterministic wrapper around
[`ddiscretestep()`](https://primarycensored.epinowcast.org/dev/reference/ddiscretestep.md).
It earns its keep as a fitting parameterisation in
[`fitdistdoublecens()`](https://primarycensored.epinowcast.org/dev/reference/fitdistdoublecens.md)
because the random walk on the logit hazard smooths the recovered PMF.

## See also

[`ddiscretestep()`](https://primarycensored.epinowcast.org/dev/reference/ddiscretestep.md),
[`hazards_to_pmf()`](https://primarycensored.epinowcast.org/dev/reference/hazards_to_pmf.md),
[`pmf_to_hazards()`](https://primarycensored.epinowcast.org/dev/reference/pmf_to_hazards.md),
[`fitdistdoublecens()`](https://primarycensored.epinowcast.org/dev/reference/fitdistdoublecens.md)

Other pdiscretehazard:
[`discretehazard_start()`](https://primarycensored.epinowcast.org/dev/reference/discretehazard_start.md),
[`pdiscretehazard()`](https://primarycensored.epinowcast.org/dev/reference/pdiscretehazard.md),
[`rdiscretehazard()`](https://primarycensored.epinowcast.org/dev/reference/rdiscretehazard.md)

## Examples

``` r
hazards <- c(0.3, 0.5, 1)
ddiscretehazard(c(1, 2, 3), boundaries = 0:3, hazards = hazards)
#> [1] 0.30 0.35 0.35
```
