# Truncated logistic distribution functions

Density, distribution function, and random generation for the logistic
distribution truncated to the interval \[min, max\]. As a primary event
distribution it describes a cumulative exposure probability that
increases smoothly and then saturates, for example an environmental
source that becomes established.

## Usage

``` r
dtlogis(x, min = 0, max = 1, location = 0, scale = 1, log = FALSE)

ptlogis(
  q,
  min = 0,
  max = 1,
  location = 0,
  scale = 1,
  lower.tail = TRUE,
  log.p = FALSE
)

rtlogis(n, min = 0, max = 1, location = 0, scale = 1)
```

## Arguments

- x, q:

  Vector of quantiles.

- min:

  Minimum value of the distribution range. Default is 0.

- max:

  Maximum value of the distribution range. Default is 1.

- location:

  Location of the logistic distribution before truncation, its midpoint,
  on the same scale as `min` and `max`. Default is 0.

- scale:

  Scale of the logistic distribution before truncation. A larger scale
  gives a slower transition. It must be positive. Default is 1.

- log, log.p:

  Logical; if TRUE, probabilities p are given as log(p).

- lower.tail:

  Logical; if TRUE (default), probabilities are P\[X \<= x\], otherwise,
  P\[X \> x\].

- n:

  Number of observations. If `length(n) > 1`, the length is taken to be
  the number required.

## Value

`dtlogis` gives the density, `ptlogis` gives the distribution function,
and `rtlogis` generates random deviates.

The length of the result is determined by `n` for `rtlogis`, and is the
maximum of the lengths of the numerical arguments for the other
functions.

## Details

With the logistic CDF \\L(z) = 1 / (1 + \exp\\-(z - \text{location}) /
\text{scale}\\)\\ and the mass \\D_L = L(max) - L(min)\\ of the
interval, the probability density function on \[min, max\] is

\$\$f(x) = \frac{L'(x)}{D_L}\$\$

and the cumulative distribution function is

\$\$F(x) = \frac{L(x) - L(min)}{D_L}.\$\$

As a primary event distribution the window is \[0, `pwindow`\], so `min`
is 0 and `max` is the window.

The differences of the logistic CDF are computed on the log scale from
the lower tails where the interval is below `location`, and from the
upper tails where it is above. This keeps \\D_L\\ and the CDF accurate
when the interval is far from `location` and \\L(max) - L(min)\\ would
otherwise round to zero.

For random number generation, we use the inverse transform method with
the same tails. A uniform draw \\u\\ gives \\x = L^{-1}(L(min) + u
D_L)\\, evaluated on the log scale. With a large `scale` relative to the
interval the distribution approaches the uniform distribution.

## See also

Primary event distributions
[`expgrowth`](https://primarycensored.epinowcast.org/reference/expgrowth.md),
[`tgumbel`](https://primarycensored.epinowcast.org/reference/tgumbel.md)

## Examples

``` r
x <- seq(0, 1, by = 0.1)
dens <- dtlogis(x, location = 0.5, scale = 0.2)
cumprobs <- ptlogis(x, location = 0.5, scale = 0.2)
samples <- rtlogis(100, location = 0.5, scale = 0.2)
```
