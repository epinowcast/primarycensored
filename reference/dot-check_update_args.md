# Validate the arguments passed to `update.pcens()`

Validate the arguments passed to
[`update.pcens()`](https://primarycensored.epinowcast.org/reference/update.pcens.md)

## Usage

``` r
.check_update_args(object, new_args, primary_args)
```

## Arguments

- object:

  A `pcens` object as created by
  [`new_pcens()`](https://primarycensored.epinowcast.org/reference/new_pcens.md).

- new_args:

  Named list of delay parameters, from `...`.

- primary_args:

  Optional named list of primary event distribution arguments. These are
  merged into `object$primary_args` in the same way as `...` is merged
  into `object$args`. Defaults to `NULL`, which leaves the primary event
  distribution arguments unchanged.

## Value

`NULL` invisibly. Called for its errors.
