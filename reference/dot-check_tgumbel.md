# Check the parameters and bounds of the truncated Gumbel distribution

Check the parameters and bounds of the truncated Gumbel distribution

## Usage

``` r
.check_tgumbel(min, max, mu, beta)
```

## Arguments

- min:

  Minimum value of the distribution range. Default is 0.

- max:

  Maximum value of the distribution range. Default is 1.

- mu:

  Location of the Gumbel distribution before truncation, on the same
  scale as `x`.

- beta:

  Scale of the Gumbel distribution before truncation, positive.

## Value

`NULL` invisibly. Called for its error.
