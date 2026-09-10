# Look up the primary event CDF from the registry

Given a primary event density function `dprimary`, looks up the
corresponding CDF function from `pcd_primary_distributions` using the
`"name"` attribute. Returns `NULL` silently when no match is found so
that callers can fall back to numerical integration.

## Usage

``` r
.lookup_pprimary(dprimary)
```

## Arguments

- dprimary:

  Function. The primary event density function.

## Value

A function (the primary CDF) or `NULL`.
