# Hazard-parameterised piecewise-constant CDF

Returns the CDF of a discrete distribution specified by its bin-wise
discrete-time hazards. Converts `hazards` to a PMF via
[`hazards_to_pmf()`](https://primarycensored.epinowcast.org/reference/hazards_to_pmf.md)
then delegates to
[`pdiscretestep()`](https://primarycensored.epinowcast.org/reference/pdiscretestep.md).

## Usage

``` r
pdiscretehazard(q, boundaries = NULL, hazards)
```

## Arguments

- q:

  Numeric vector of quantiles.

- boundaries:

  Numeric vector of length \\K+1\\ defining the bin edges. Must be
  strictly increasing. Defaults to `0:K`.

- hazards:

  Numeric vector of length \\K\\ (or \\K-1\\, in which case a trailing 1
  is appended) giving the discrete-time conditional hazard for each bin.
  Values must be in \\\[0, 1\]\\.

## Value

Numeric vector of CDF values, the same length as `q`.

## Details

Outside fitting it is a deterministic wrapper around
[`pdiscretestep()`](https://primarycensored.epinowcast.org/reference/pdiscretestep.md).
It earns its keep as a fitting parameterisation in
[`fitdistdoublecens()`](https://primarycensored.epinowcast.org/reference/fitdistdoublecens.md)
because the random walk on the logit hazard smooths the recovered PMF.

## See also

[`pdiscretestep()`](https://primarycensored.epinowcast.org/reference/pdiscretestep.md),
[`hazards_to_pmf()`](https://primarycensored.epinowcast.org/reference/hazards_to_pmf.md),
[`pmf_to_hazards()`](https://primarycensored.epinowcast.org/reference/pmf_to_hazards.md),
[`fitdistdoublecens()`](https://primarycensored.epinowcast.org/reference/fitdistdoublecens.md)

Other pdiscretehazard:
[`ddiscretehazard()`](https://primarycensored.epinowcast.org/reference/ddiscretehazard.md),
[`discretehazard_start()`](https://primarycensored.epinowcast.org/reference/discretehazard_start.md),
[`rdiscretehazard()`](https://primarycensored.epinowcast.org/reference/rdiscretehazard.md)

## Examples

``` r
hazards <- c(0.3, 0.5, 1)
pdiscretehazard(c(0.5, 1, 2, 3), boundaries = 0:3, hazards = hazards)
#> [1] 0.00 0.30 0.65 1.00
```
