# Resolve a delay distribution function from a name or function

Accepts either a function (returned as-is, with its existing `"name"`
attribute preserved) or a character string that is looked up against
[`pcd_distributions`](https://primarycensored.epinowcast.org/dev/reference/pcd_distributions.md).
When a string is supplied, the corresponding base R `p<name>` function
is returned with the `"name"` attribute attached so analytical solutions
can dispatch.

## Usage

``` r
.resolve_pdist(pdist, type = c("p", "d"))
```

## Arguments

- pdist:

  Either a function or a character string.

- type:

  Character string. `"p"` for CDF lookup, `"d"` for density. Defaults to
  `"p"`.

## Value

A function with a `"name"` attribute.
