# Resolve the primary CDF, validating against `dprimary` if both supplied

Returns the primary CDF to use. If the user supplies `pprimary`
explicitly (either a function or a string name), it is returned (after
resolving the string via
[`pcd_dist_name`](https://primarycensored.epinowcast.org/reference/pcd_dist_name.md)).
When both `dprimary` and `pprimary` carry a `"name"` attribute, the
names must agree on everything except the leading `d`/`p`; otherwise we
error to catch typos like `dunif` + `pexpgrowth`. If `pprimary` is not
supplied, falls back to a registry lookup against `dprimary` via
[`.lookup_pprimary`](https://primarycensored.epinowcast.org/reference/dot-lookup_pprimary.md),
which may return `NULL`.

## Usage

``` r
.resolve_pprimary(dprimary, pprimary = NULL)
```

## Arguments

- dprimary:

  The primary density function.

- pprimary:

  Optional user-supplied primary CDF (function or string).

## Value

A primary CDF function, or `NULL` if no match was found.
