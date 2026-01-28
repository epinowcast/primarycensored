# Normalise a primary event censored CDF

Internal function to normalise a primary event censored CDF when
truncation is applied. The CDF is normalised using (F(q) - F(L)) /
(F(D) - F(L)) and values outside \[L, D\] are clamped to 0 or 1.

## Usage

``` r
.normalise_cdf(result, q, L, D, pcens_obj, pwindow)
```

## Arguments

- result:

  Numeric vector of CDF values to normalise.

- q:

  Numeric vector of quantiles at which CDF was evaluated.

- L:

  Numeric lower truncation point

- D:

  Numeric upper truncation point

- pcens_obj:

  A primarycensored object as created by
  [`new_pcens()`](https://primarycensored.epinowcast.org/dev/reference/new_pcens.md).

- pwindow:

  Secondary event window

## Value

Normalised CDF values as a numeric vector
