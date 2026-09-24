# Look up the primary event CDF from the registry

Given a primary event density function `dprimary`, looks up the
corresponding CDF function from `pcd_primary_distributions` using the
`"name"` attribute. Returns `NULL` silently when no match is found so
that callers can fall back to numerical integration.

## Usage

``` r
.lookup_pprimary(dprimary, dprim_name = .dist_name(dprimary))
```

## Arguments

- dprimary:

  Function. The primary event density function.

- dprim_name:

  Name of `dprimary`, as given by
  [`.dist_name()`](https://primarycensored.epinowcast.org/dev/reference/dot-dist_name.md).

## Value

A function (the primary CDF) or `NULL`.
