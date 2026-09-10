# Resolve `primary_args` / `dprimary_args` with a deprecation

The resolver distinguishes "user supplied nothing" from "user supplied
an empty list" so the deprecated `dprimary_args` path can be detected.
The returned list is never `NULL`.

## Usage

``` r
.resolve_primary_args(
  primary_args,
  dprimary_args,
  fn,
  env = parent.frame(),
  user_env = parent.frame(2)
)
```

## Arguments

- primary_args:

  The new argument value (or `NULL`).

- dprimary_args:

  The old argument value (or `NULL`).

- fn:

  Character string identifying the calling function (used in the
  deprecation message).

- env:

  Environment of the exported function that owns the deprecated
  argument. Defaults to the caller of this helper.

- user_env:

  Environment the exported function was called from. Defaults to the
  caller of the caller of this helper.

## Value

A list (possibly empty) of primary distribution arguments.

## Details

The deprecation is soft
([`lifecycle::deprecate_soft()`](https://lifecycle.r-lib.org/reference/deprecate_soft.html)):
a warning is shown when the exported function is called from the global
environment or from the package under test, and calls from other
packages stay silent.
