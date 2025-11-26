# Deprecation name helper

Provides lifecycle management consistently across several functions.
Currently uses
[`lifecycle::deprecate_soft()`](https://lifecycle.r-lib.org/reference/deprecate_soft.html) -
i.e. warns only when used directly. In future versions, this will switch
to warning unconditionally with
[`lifecycle::deprecate_warn()`](https://lifecycle.r-lib.org/reference/deprecate_soft.html),
then throwing via
[`lifecycle::deprecate_warn()`](https://lifecycle.r-lib.org/reference/deprecate_soft.html),
and finally be deleted along with the subject arguments.

## Usage

``` r
.name_deprecation(
  pdist_name,
  dprimary_name,
  env = rlang::caller_env(),
  user_env = rlang::caller_env(2)
)
```

## Arguments

- pdist_name:

  the deprecated variable to check

- dprimary_name:

  the deprecated variable to check

- env, user_env:

  Pair of environments that define where `deprecate_*()` was called
  (used to determine the package name) and where the function called the
  deprecating function was called (used to determine if
  `deprecate_soft()` should message).

  These are only needed if you're calling `deprecate_*()` from an
  internal helper, in which case you should forward `env = caller_env()`
  and `user_env = caller_env(2)`.
