# Check the delays and windows of a uniform primary CDF

Check the delays and windows of a uniform primary CDF

## Usage

``` r
.check_pwindow(q, pwindow)
```

## Arguments

- q:

  Vector of quantiles

- pwindow:

  Primary event window. Use `pwindow = 0` for an exactly observed
  primary event, in which case the delay CDF is used directly.

## Value

Invisible NULL if valid, otherwise stops with an error message.
