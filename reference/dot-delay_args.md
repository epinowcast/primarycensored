# Delay distribution arguments an analytical solution handles

Delay distribution arguments an analytical solution handles

## Usage

``` r
.delay_args(object, allowed, use_numeric = FALSE)
```

## Arguments

- object:

  A `primarycensored` object as created by
  [`new_pcens()`](https://primarycensored.epinowcast.org/reference/new_pcens.md).

- allowed:

  Character vector of the argument names the solution uses.

- use_numeric:

  Logical, if TRUE forces use of numeric integration even for
  distributions with analytical solutions. This is primarily useful for
  testing purposes or for settings where the analytical solution breaks
  down.

## Value

The delay arguments, or `NULL` if the numerical method should be used.
