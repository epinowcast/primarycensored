# Exponential growth distribution functions

Density, distribution function, and random generation for the
exponential growth distribution.

## Usage

``` r
dexpgrowth(x, min = 0, max = 1, r, log = FALSE)

pexpgrowth(q, min = 0, max = 1, r, lower.tail = TRUE, log.p = FALSE)

rexpgrowth(n, min = 0, max = 1, r)
```

## Arguments

- x, q:

  Vector of quantiles.

- min:

  Minimum value of the distribution range. Default is 0.

- max:

  Maximum value of the distribution range. Default is 1.

- r:

  Rate parameter for the exponential growth.

- log, log.p:

  Logical; if TRUE, probabilities p are given as log(p).

- lower.tail:

  Logical; if TRUE (default), probabilities are P\[X \<= x\], otherwise,
  P\[X \> x\].

- n:

  Number of observations. If `length(n) > 1`, the length is taken to be
  the number required.

## Value

`dexpgrowth` gives the density, `pexpgrowth` gives the distribution
function, and `rexpgrowth` generates random deviates.

The length of the result is determined by `n` for `rexpgrowth`, and is
the maximum of the lengths of the numerical arguments for the other
functions.

## Details

The exponential growth distribution is defined on the interval \[min,
max\] with rate parameter (r). Its probability density function (PDF)
is:

\$\$f(x) = \frac{r \cdot \exp(r \cdot (x - min))}{\exp(r \cdot max) -
\exp(r \cdot min)}\$\$

The cumulative distribution function (CDF) is:

\$\$F(x) = \frac{\exp(r \cdot (x - min)) - \exp(r \cdot min)}{ \exp(r
\cdot max) - \exp(r \cdot min)}\$\$

For random number generation, we use the inverse transform sampling
method:

1.  Generate \\u \sim \text{Uniform}(0,1)\\

2.  Set \\F(x) = u\\ and solve for \\x\\: \$\$ x = min + \frac{1}{r}
    \cdot \log(u \cdot (\exp(r \cdot max) - \exp(r \cdot min)) + \exp(r
    \cdot min)) \$\$

This method works because of the probability integral transform theorem,
which states that if \\X\\ is a continuous random variable with CDF
\\F(x)\\, then \\Y = F(X)\\ follows a \\\text{Uniform}(0,1)\\
distribution. Conversely, if \\U\\ is a \\\text{Uniform}(0,1)\\ random
variable, then \\F^{-1}(U)\\ has the same distribution as \\X\\, where
\\F^{-1}\\ is the inverse of the CDF.

In our case, we generate \\u\\ from \\\text{Uniform}(0,1)\\, then solve
\\F(x) = u\\ for \\x\\ to get a sample from our exponential growth
distribution. The formula for \\x\\ is derived by algebraically solving
the equation:

\$\$ u = \frac{\exp(r \cdot (x - min)) - \exp(r \cdot min)}{\exp(r \cdot
max) - \exp(r \cdot min)} \$\$

When \\r\\ is very close to 0 (\\\|r\| \< 1e-10\\), the distribution
approximates a uniform distribution on \[min, max\], and we use a
simpler method to generate samples directly from this uniform
distribution.

## Examples

``` r
x <- seq(0, 1, by = 0.1)
probs <- dexpgrowth(x, r = 0.2)
cumprobs <- pexpgrowth(x, r = 0.2)
samples <- rexpgrowth(100, r = 0.2)
```
