# Group observations that share censoring and truncation settings

Group observations that share censoring and truncation settings

## Usage

``` r
.param_groups(params, cols)
```

## Arguments

- params:

  A data frame of per-observation settings.

- cols:

  Names of the columns to group by.

## Value

A list with one element per unique combination of `cols`. Each element
is a list of the values of `cols` and a logical `mask` selecting the
rows of `params` with those values.
