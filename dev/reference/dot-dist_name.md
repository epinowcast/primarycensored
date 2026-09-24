# Get the distribution name of a function

Returns the `"name"` attribute of `func` if set. A `stats` function that
is identical to one named in
[pcd_distributions](https://primarycensored.epinowcast.org/dev/reference/pcd_distributions.md)
or
[pcd_primary_distributions](https://primarycensored.epinowcast.org/dev/reference/pcd_primary_distributions.md)
gets that name. Otherwise the name is found with
[`.extract_function_name()`](https://primarycensored.epinowcast.org/dev/reference/dot-extract_function_name.md),
which deparses the function body and is slower.

## Usage

``` r
.dist_name(func)
```

## Arguments

- func:

  Function, for example the `p`- or `d`- form of a distribution
  function.

## Value

Character string with the name of the function, or `"unknown"`.
