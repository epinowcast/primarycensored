# Extract base function name

This helper function extracts the base name of a function, removing
namespace prefixes. Base R distribution functions are identified from
the C routine they call. Functions exported by another package (for
example
[`flexsurv::pgengamma.orig()`](http://chjackson.github.io/flexsurv-dev/reference/GenGamma.orig.md))
are identified by their exported name.

## Usage

``` r
.extract_function_name(func)
```

## Arguments

- func:

  Function, for example the `p`- or `d`- form of a distribution
  function.

## Value

Character string representing the base name of the function, or
`"unknown"` if it cannot be determined.
