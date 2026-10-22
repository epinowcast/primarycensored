# Log of a difference between two points from either tail

With `L(t) + U(t)` constant, `L(q) - L(y)` equals `U(y) - U(q)`. This
uses the form with the smaller ratio of the terms.

## Usage

``` r
.exptilt_tail_diff(lower_q, lower_y, upper_q, upper_y)
```

## Arguments

- lower_q, lower_y, upper_q, upper_y:

  Log lower and upper tail quantities at `q` and `y`.

## Value

Vector of the log differences.
