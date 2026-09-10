# Test whether a numeric vector lies on the simplex

Returns `TRUE` when the vector contains no missing values, no negative
entries, and sums to 1 within \\10^{-8}\\. Used by the `discretestep`
family to apply a soft penalty inside fitting closures rather than
erroring.

## Usage

``` r
.is_valid_simplex(p, tol = 1e-08)
```
