# Resolve a registry name or alias to its function name

Resolve a registry name or alias to its function name

## Usage

``` r
.canonical_name(name, registry, column)
```

## Arguments

- name:

  Character string, a distribution name as given by
  [`.dist_name()`](https://primarycensored.epinowcast.org/reference/dot-dist_name.md).

- registry:

  A registry data frame,
  [pcd_distributions](https://primarycensored.epinowcast.org/reference/pcd_distributions.md)
  or
  [pcd_primary_distributions](https://primarycensored.epinowcast.org/reference/pcd_primary_distributions.md).

- column:

  The registry column holding the function names, `"pdist"` or
  `"dprimary"`.

## Value

The function name in `column` for the registry row whose `name` or
`aliases` is `name`, or `name` itself if it is already a function name
or is not in the registry.
