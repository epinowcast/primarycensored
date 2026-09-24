# Look up the density matching a delay CDF

Finds the `d` function paired with a `p` function by name, for example
[`dgamma()`](https://rdrr.io/r/stats/GammaDist.html) for
[`pgamma()`](https://rdrr.io/r/stats/GammaDist.html). The name is taken
from the `"name"` attribute of `pdist` or inferred with
[`.dist_name()`](https://primarycensored.epinowcast.org/dev/reference/dot-dist_name.md).
The density is searched for from the environment of `pdist`, then in
`stats` and `primarycensored`.

## Usage

``` r
.lookup_ddist(pdist)
```

## Arguments

- pdist:

  Delay distribution CDF.

## Value

The density function. An error is raised if none is found.
