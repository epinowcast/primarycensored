# Validate truncation bounds in a data frame

Internal function to validate that L (lower truncation) is less than D
(upper truncation) for all rows in a data frame.

## Usage

``` r
.check_truncation_bounds_df(data, L_col, D_col)
```

## Arguments

- data:

  Data frame containing the L and D columns

- L_col:

  Name of the column containing L values

- D_col:

  Name of the column containing D values

## Value

Invisible NULL if valid, otherwise stops with an error message.
