# Methods for log-logistic delays

Analytical primary event censored CDFs for a log-logistic delay with a
uniform primary event window. They honour `use_numeric`. The numerical
method is used where the closed form does not apply or is ill
conditioned.

## Usage

``` r
# S3 method for class 'pcens_pllogis_dunif'
pcens_cdf(object, q, pwindow, use_numeric = FALSE)
```

## Arguments

- object:

  A `primarycensored` object as created by
  [`new_pcens()`](https://primarycensored.epinowcast.org/reference/new_pcens.md).

- q:

  Vector of quantiles

- pwindow:

  Primary event window. Use `pwindow = 0` for an exactly observed
  primary event, in which case the delay CDF is used directly.

- use_numeric:

  Logical, if TRUE forces use of numeric integration even for
  distributions with analytical solutions. This is primarily useful for
  testing purposes or for settings where the analytical solution breaks
  down.

## Value

Vector of computed primary event censored CDFs

## Details

The delay has \\F(x) = 1 / (1 + (x / \lambda)^{-k})\\, given as in
[`flexsurv::pllogis()`](http://chjackson.github.io/flexsurv-dev/reference/Llogis.md)
and
[`actuar::pllogis()`](https://rdrr.io/pkg/actuar/man/Loglogistic.html)
with `shape` and `scale`, or `rate` for the latter. The Stan functions
use `[scale, shape]`, the order of Stan's `loglogistic_cdf()`.

With \\G_1(t) = \int_0^t F(u) du\\ the CDF is \\\\G_1(q) - G_1(q - w)\\
/ w\\, evaluated with series for a shape of at least 0.01. The numerical
method is used where this difference loses precision, so the CDF can
step by about 1e-8 of the smaller tail at the switch. See
[`vignette("analytic-solutions")`](https://primarycensored.epinowcast.org/articles/analytic-solutions.md)
for the derivation and the switching rule.

## Examples

``` r
# Log-logistic delay with a uniform primary
obj <- new_pcens(
  pdist = flexsurv::pllogis, dprimary = dunif,
  primary_args = list(), shape = 2, scale = 5
)
pcens_cdf(obj, q = c(1, 4, 8), pwindow = 2)
#> [1] 0.0065111 0.2644136 0.6596526
```
