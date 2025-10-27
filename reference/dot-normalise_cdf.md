# Normalise a primary event censored CDF

Internal function to normalise a primary event censored CDF when
truncation is applied. The CDF is normalised by dividing by its value at
the truncation point D and setting all values beyond D to 1.

## Usage

``` r
.normalise_cdf(result, q, D, pcens_obj, pwindow)
```

## Arguments

- result:

  Numeric vector of CDF values to normalise.

- q:

  Numeric vector of quantiles at which CDF was evaluated.

- D:

  Numeric truncation point

- pcens_obj:

  A primarycensored object as created by
  [`new_pcens()`](https://primarycensored.epinowcast.org/reference/new_pcens.md).

- pwindow:

  Secondary event window

## Value

Normalised CDF values as a numeric vector
