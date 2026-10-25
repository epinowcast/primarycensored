# Truncated Gumbel distribution functions

Density, distribution function, and random generation for a Gumbel
distribution truncated to `[min, max]`. It is a primary event
distribution for exposure that builds up or decays unevenly within the
window.

## Usage

``` r
dtgumbel(x, min = 0, max = 1, mu, beta, log = FALSE)

ptgumbel(q, min = 0, max = 1, mu, beta, lower.tail = TRUE, log.p = FALSE)

rtgumbel(n, min = 0, max = 1, mu, beta)
```

## Arguments

- x, q:

  Vector of quantiles.

- min:

  Minimum value of the distribution range. Default is 0.

- max:

  Maximum value of the distribution range. Default is 1.

- mu:

  Location of the Gumbel distribution before truncation, on the same
  scale as `x`.

- beta:

  Scale of the Gumbel distribution before truncation, positive.

- log, log.p:

  Logical; if TRUE, probabilities p are given as log(p).

- lower.tail:

  Logical; if TRUE (default), probabilities are P\[X \<= x\], otherwise,
  P\[X \> x\].

- n:

  Number of observations. If `length(n) > 1`, the length is taken to be
  the number required.

## Value

`dtgumbel` gives the density, `ptgumbel` gives the distribution
function, and `rtgumbel` generates random deviates.

The length of the result is determined by `n` for `rtgumbel`, and is the
maximum of the lengths of the numerical arguments for the other
functions.

## Details

With \\G(z) = \exp\[-\exp\\-(z - \mu) / \beta\\\]\\, the CDF of the
Gumbel distribution, the density on \[min, max\] is \$\$f(x) =
\frac{G'(x)}{G(max) - G(min)}\$\$ and the cumulative distribution
function is \$\$F(x) = \frac{G(x) - G(min)}{G(max) - G(min)}.\$\$ The
density peaks at \\x = \mu\\ and falls more slowly above it. If \\\mu\\
is far below `min` it decays with rate \\1 / \beta\\, as for
[`dexpgrowth()`](https://primarycensored.epinowcast.org/reference/expgrowth.md)
with \\r = -1 / \beta\\. If \\\mu\\ is far above `max` it is a narrow
spike at `max`.

Differences of \\s(x) = \exp\\-(x - \mu) / \beta\\\\ are evaluated on
the log scale, so the tails stay precise when \\\mu\\ is far from the
window. Random numbers are drawn by inverting the upper tail.

The Stan equivalents are `tgumbel_lpdf()`, `tgumbel_lcdf()` and
`tgumbel_rng()`, with the primary distribution identifier 4 and
`primary_params = [mu, beta]`. Analytical primary event censored CDFs
are described in
[pcens_cdf_gumbel](https://primarycensored.epinowcast.org/reference/pcens_cdf_gumbel.md).

## See also

Primary event distributions
[`expgrowth`](https://primarycensored.epinowcast.org/reference/expgrowth.md),
[`tlogis`](https://primarycensored.epinowcast.org/reference/tlogis.md)

## Examples

``` r
x <- seq(0, 2, by = 0.25)
dens <- dtgumbel(x, min = 0, max = 2, mu = 0.5, beta = 0.3)
cumprobs <- ptgumbel(x, min = 0, max = 2, mu = 0.5, beta = 0.3)
samples <- rtgumbel(100, min = 0, max = 2, mu = 0.5, beta = 0.3)
```
