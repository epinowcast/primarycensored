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
/ w\\, where \\G_1(t) = t F(t) (1 - r\_{1/k}(A))\\, \\A = (t /
\lambda)^k\\ and \$\$r_a(A) = \frac{1 + A}{A^{a + 1}} \int_0^A
\frac{y^a}{(1 + y)^2} dy.\$\$ The ratio is at most \\1 / (a + 1)\\. Base
R has no Gauss hypergeometric function, and the incomplete beta form
holds only for \\a \< 1\\. `.loglogistic_ratio()` therefore uses three
series that do not cancel, for \\A \le 1\\, \\1 \< A \le 3\\ and \\A \>
3\\. Each unique endpoint is evaluated once.

**Fallback to the numerical method.** The solution is used for a shape
of at least 0.01. Within this limit a bound on the rounding error of the
difference of \\G_1\\ is compared with the smaller of the CDF and the
survival. The numerical method is used where it exceeds 1e-8 of that
tail, or 1e-6 of the density of the CDF. The bound is not taken below
1e-15 in the upper tail. The same rule is used in Stan. The numerical
method integrates the delay CDF, or the survival where the delay CDF at
the mean primary event time is above a half, with a relative tolerance
of 1e-13. The switch between methods makes the CDF discontinuous by
about 1e-8 of the smaller tail.

## Examples

``` r
# Log-logistic delay with a uniform primary
pllogis <- add_name_attribute(
  function(q, shape, scale) {
    plogis(shape * (log(pmax(q, 0)) - log(scale)))
  },
  "pllogis"
)
obj <- new_pcens(
  pdist = pllogis, dprimary = dunif,
  primary_args = list(), shape = 2, scale = 5
)
pcens_cdf(obj, q = c(1, 4, 8), pwindow = 2)
#> [1] 0.0065111 0.2644136 0.6596526
```
