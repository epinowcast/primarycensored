# Start values for the logit-hazard parameterisation

Builds a named list of starting values for
[`fitdistdoublecens()`](https://primarycensored.epinowcast.org/dev/reference/fitdistdoublecens.md)
with `distr = "discretehazard"`. The free parameters are the logit
intercept `alpha`, the log random-walk (or random-effect) scale
`log_sigma`, and the `K - 1` innovations `eps_1, ..., eps_{K-1}`. The
final bin hazard is pinned to 1 inside the parameterisation so the
implied PMF sums to 1.

## Usage

``` r
discretehazard_start(K, alpha = -2, log_sigma = log(1), eps = 0)
```

## Arguments

- K:

  Integer, number of bins in the hazard parameterisation.

- alpha:

  Numeric, start value for the logit intercept.

- log_sigma:

  Numeric, start value for the log scale.

- eps:

  Numeric, scalar or length-`K - 1` vector of start values for the
  innovations `eps_1, ..., eps_{K-1}`.

## Value

Named list suitable for the `start` argument of
[`fitdistdoublecens()`](https://primarycensored.epinowcast.org/dev/reference/fitdistdoublecens.md).

## See also

Other pdiscretehazard:
[`ddiscretehazard()`](https://primarycensored.epinowcast.org/dev/reference/ddiscretehazard.md),
[`pdiscretehazard()`](https://primarycensored.epinowcast.org/dev/reference/pdiscretehazard.md),
[`rdiscretehazard()`](https://primarycensored.epinowcast.org/dev/reference/rdiscretehazard.md)

## Examples

``` r
discretehazard_start(K = 5)
#> $alpha
#> [1] -2
#> 
#> $log_sigma
#> [1] 0
#> 
#> $eps_1
#> [1] 0
#> 
#> $eps_2
#> [1] 0
#> 
#> $eps_3
#> [1] 0
#> 
#> $eps_4
#> [1] 0
#> 
```
