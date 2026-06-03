# Validate truncation bounds L and D

Internal function to validate that L (lower truncation) and D (upper
truncation) parameters are valid. L must be less than D. L may be
negative or `-Inf` for delay distributions with support below zero.

## Usage

``` r
.check_truncation_bounds(L, D)
```

## Arguments

- L:

  Lower truncation bound (may be negative or `-Inf`)

- D:

  Upper truncation bound

## Value

Invisible NULL if valid, otherwise stops with an error message.
