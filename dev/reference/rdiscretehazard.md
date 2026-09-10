# Sample from a hazard-parameterised step distribution

Draws `n` independent samples from the discrete distribution defined by
`hazards`. Converts `hazards` to a PMF via
[`hazards_to_pmf()`](https://primarycensored.epinowcast.org/dev/reference/hazards_to_pmf.md)
then delegates to
[`rdiscretestep()`](https://primarycensored.epinowcast.org/dev/reference/rdiscretestep.md).

## Usage

``` r
rdiscretehazard(n, boundaries = NULL, hazards)
```

## Arguments

- n:

  Integer. Number of samples to draw.

- boundaries:

  Numeric vector of length \\K+1\\ defining the bin edges. Must be
  strictly increasing. Defaults to `0:K`.

- hazards:

  Numeric vector of length \\K\\ (or \\K-1\\, in which case a trailing 1
  is appended) giving the discrete-time conditional hazard for each bin.
  Values must be in \\\[0, 1\]\\.

## Value

Numeric vector of length `n`.

## Details

Outside fitting it is a deterministic wrapper around
[`rdiscretestep()`](https://primarycensored.epinowcast.org/dev/reference/rdiscretestep.md).
It earns its keep as a fitting parameterisation in
[`fitdistdoublecens()`](https://primarycensored.epinowcast.org/dev/reference/fitdistdoublecens.md)
because the random walk on the logit hazard smooths the recovered PMF.

## See also

[`rdiscretestep()`](https://primarycensored.epinowcast.org/dev/reference/rdiscretestep.md),
[`hazards_to_pmf()`](https://primarycensored.epinowcast.org/dev/reference/hazards_to_pmf.md),
[`pmf_to_hazards()`](https://primarycensored.epinowcast.org/dev/reference/pmf_to_hazards.md),
[`fitdistdoublecens()`](https://primarycensored.epinowcast.org/dev/reference/fitdistdoublecens.md)

Other pdiscretehazard:
[`ddiscretehazard()`](https://primarycensored.epinowcast.org/dev/reference/ddiscretehazard.md),
[`discretehazard_start()`](https://primarycensored.epinowcast.org/dev/reference/discretehazard_start.md),
[`pdiscretehazard()`](https://primarycensored.epinowcast.org/dev/reference/pdiscretehazard.md)

## Examples

``` r
set.seed(42)
rdiscretehazard(10, boundaries = 0:3, hazards = c(0.3, 0.5, 1))
#>  [1] 1 1 2 1 3 3 1 2 3 1
```
