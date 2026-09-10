# Method for generalised gamma delay with uniform primary

Analytical solution for the generalised gamma distribution in the Stacy
parameterisation used by
[`flexsurv::pgengamma.orig()`](http://chjackson.github.io/flexsurv-dev/reference/GenGamma.orig.md),
with parameters `shape`, `scale` and `k`. The delay CDF is \\F_T(t) =
P(k, (t / \theta)^a)\\ with \\P\\ the regularised lower incomplete gamma
function, \\a\\ the `shape` and \\\theta\\ the `scale`. The mean is
\\E\[T\] = \theta \Gamma(k + 1/a) / \Gamma(k)\\ and the partial
expectation distribution is \\\tilde F_T(t) = P(k + 1/a, (t /
\theta)^a)\\, so the solution generalises the gamma (`shape = 1`) and
Weibull (`k = 1`) cases. See
[`vignette("analytic-solutions")`](https://primarycensored.epinowcast.org/dev/articles/analytic-solutions.md)
for the derivation.

## Usage

``` r
# S3 method for class 'pcens_pgengamma.orig_dunif'
pcens_cdf(object, q, pwindow, use_numeric = FALSE)
```

## Arguments

- object:

  A `primarycensored` object as created by
  [`new_pcens()`](https://primarycensored.epinowcast.org/dev/reference/new_pcens.md).

- q:

  Vector of quantiles

- pwindow:

  Primary event window

- use_numeric:

  Logical, if TRUE forces use of numeric integration even for
  distributions with analytical solutions. This is primarily useful for
  testing purposes or for settings where the analytical solution breaks
  down.

## Value

Vector of computed primary event censored CDFs

## See also

Low level primary event censored distribution objects and methods
[`new_pcens()`](https://primarycensored.epinowcast.org/dev/reference/new_pcens.md),
[`pcens_cdf()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.md),
[`pcens_cdf.default()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.default.md),
[`pcens_cdf.pcens_pdiscretehazard()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_pdiscretehazard.md),
[`pcens_cdf.pcens_pdiscretestep()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_pdiscretestep.md),
[`pcens_cdf.pcens_pgamma_dunif()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_pgamma_dunif.md),
[`pcens_cdf.pcens_pgengamma_dunif()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_pgengamma_dunif.md),
[`pcens_cdf.pcens_plnorm_dunif()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_plnorm_dunif.md),
[`pcens_cdf.pcens_pweibull_dunif()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_pweibull_dunif.md),
[`pcens_quantile()`](https://primarycensored.epinowcast.org/dev/reference/pcens_quantile.md),
[`pcens_quantile.default()`](https://primarycensored.epinowcast.org/dev/reference/pcens_quantile.default.md)

## Examples

``` r
pcens_obj <- new_pcens(
  pdist = flexsurv::pgengamma.orig,
  dprimary = dunif,
  dprimary_args = list(min = 0, max = 1),
  shape = 1.5,
  scale = 2,
  k = 0.8
)
#> Warning: The `dprimary_args` argument of `new_pcens()` is deprecated as of
#> primarycensored 1.6.0.
#> ℹ Please use the `primary_args` argument instead.
pcens_cdf(pcens_obj, q = c(1, 4, 8), pwindow = 1)
#> [1] 0.1940835 0.9298869 0.9995685
```
