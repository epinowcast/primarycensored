# Validate truncation bounds L and D

Internal function to validate that L (lower truncation) and D (upper
truncation) parameters are valid: L must be non-negative and less than
D.

## Usage

``` r
.check_truncation_bounds(L, D)
```

## Arguments

- L:

  Lower truncation bound

- D:

  Upper truncation bound

## Value

Invisible NULL if valid, otherwise stops with an error message.
