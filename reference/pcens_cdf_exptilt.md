# Methods for delays with an exponentially tilted primary

Analytical primary event censored CDFs for exponential, gamma, normal,
lognormal, Weibull and generalised gamma delay distributions with an
exponentially tilted primary event window, the
[`dexpgrowth()`](https://primarycensored.epinowcast.org/reference/expgrowth.md)
primary distribution with `r` equal to the tilt \\\rho\\. They honour
`use_numeric`, and use the numerical method of
[`pcens_cdf.default()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.default.md)
when no closed form applies.

## Usage

``` r
# S3 method for class 'pcens_pexp_dexpgrowth'
pcens_cdf(object, q, pwindow, use_numeric = FALSE)

# S3 method for class 'pcens_pgamma_dexpgrowth'
pcens_cdf(object, q, pwindow, use_numeric = FALSE)

# S3 method for class 'pcens_pnorm_dexpgrowth'
pcens_cdf(object, q, pwindow, use_numeric = FALSE)

# S3 method for class 'pcens_plnorm_dexpgrowth'
pcens_cdf(object, q, pwindow, use_numeric = FALSE)

# S3 method for class 'pcens_pweibull_dexpgrowth'
pcens_cdf(object, q, pwindow, use_numeric = FALSE)

# S3 method for class 'pcens_pgengamma.orig_dexpgrowth'
pcens_cdf(object, q, pwindow, use_numeric = FALSE)

# S3 method for class 'pcens_pgengamma_dexpgrowth'
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

With window width \\w\\, delay CDF \\F\\ and \\J(x) = \int e^{-\rho u}
f(u) du\\ up to \\x\\, the CDF at \\q\\ is \$\$ F\_\rho(q) = F(q - w) +
\frac{e^{\rho q} \\J(q) - J(q - w)\\ - \\F(q) - F(q - w)\\}{e^{\rho w} -
1}. \$\$ Every term depends on one endpoint, \\q\\ or \\q - w\\, so each
endpoint is evaluated once and shared between neighbouring `q`. Terms
are on the log scale, and each difference is taken between the lower or
the upper tail terms, whichever loses less precision.

The exponential and gamma forms need \\\lambda + \rho \> 0\\ for rate
\\\lambda\\, so that the tilted delay distribution exists. Otherwise the
numerical method is used. The numerical method is less accurate in the
lower tail of a gamma with shape below 1, with a relative error of about
1e-5 in R and 1e-4 in Stan. The lognormal transform has no closed form.
It is evaluated by quadrature for \\\rho \> 0\\ and by a series of
partial moments for \\\rho \< 0\\. The numerical method is used where
the series needs too many terms or is slower than numerical integration.

The Weibull and the generalised gamma with power \\a\\, scale \\\theta\\
and \\k\\ (\\k = 1\\ for the Weibull) have the transform \\J(x) = \sum_n
(-\rho \theta)^n / n!\\ \Gamma(k + n / a) / \Gamma(k)\\ P(k + n / a, (x
/ \theta)^a)\\, for the regularised incomplete gamma function \\P\\. It
exists for every tilt, and needs \\Q \> 0\\ for
[`flexsurv::pgengamma()`](http://chjackson.github.io/flexsurv-dev/reference/GenGamma.md).
For \\\rho \> 0\\ its terms alternate in sign and cancel. The numerical
method is used at a `q` where the terms of the series cancel by more
than 1e5, where \\\rho q \ge 400\\, and in the upper tail, where the
series has no tail form and the terms of the difference between the
endpoints exceed the smaller of the CDF and its survival function by
more than 1e5.

The direct form cancels as \\\rho \to 0\\. With \\G_k(t) = \int (t -
u)^k f(u) du\\ up to \\t\\, two series in \\\rho\\ truncated after the
second order replace it.

- For \\\|\rho\| w\\ below \\10^{-2}\\, or \\10^{-3}\\ for delays on the
  reals, and \\\|\rho\| (\|q\| + w) \< 0.1\\, the uniform window limit
  with its corrections, \\\\G_1(q) - G_1(q - w)\\ / w + \rho \\G_2(q) -
  w G_1(q) - G_2(q - w) - w G_1(q - w)\\ / (2 w) + \rho^2 \\G_3(q) -
  G_3(q - w) - 3 w \[G_2(q) + G_2(q - w)\] / 2 + w^2 \[G_1(q) - G_1(q -
  w)\] / 2\\ / (6 w)\\. Its relative error is below \\(\|\rho\| w)^3 /
  160\\.

- Otherwise, for delays on the non-negative reals with \\q \< w\\ and
  \\\|\rho\| q \< 5 \times 10^{-3}\\, \\\rho \\G_1(q) + \rho G_2(q) /
  2 + \rho^2 G_3(q) / 6\\ / (e^{\rho w} - 1)\\. Its relative error is
  below \\(\|\rho\| q)^3 / 24\\.

Both errors are below \\10^{-8}\\ at these limits.

The CDF agrees with numerical integration to a relative 1e-6 in the
package tests.

## See also

Low level primary event censored distribution objects and methods
[`new_pcens()`](https://primarycensored.epinowcast.org/reference/new_pcens.md),
[`pcens_cdf()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.md),
[`pcens_cdf.default()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.default.md),
[`pcens_cdf.pcens_pdiscretehazard()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_pdiscretehazard.md),
[`pcens_cdf.pcens_pdiscretestep()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_pdiscretestep.md),
[`pcens_cdf.pcens_pgamma_dunif()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_pgamma_dunif.md),
[`pcens_cdf.pcens_pgengamma.orig_dunif()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_pgengamma.orig_dunif.md),
[`pcens_cdf.pcens_pgengamma_dunif()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_pgengamma_dunif.md),
[`pcens_cdf.pcens_plnorm_dunif()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_plnorm_dunif.md),
[`pcens_cdf.pcens_pweibull_dunif()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_pweibull_dunif.md),
[`pcens_cdf_gumbel`](https://primarycensored.epinowcast.org/reference/pcens_cdf_gumbel.md),
[`pcens_cdf_tlogis`](https://primarycensored.epinowcast.org/reference/pcens_cdf_tlogis.md),
[`pcens_pmf()`](https://primarycensored.epinowcast.org/reference/pcens_pmf.md),
[`pcens_pmf.default()`](https://primarycensored.epinowcast.org/reference/pcens_pmf.default.md),
[`pcens_quantile()`](https://primarycensored.epinowcast.org/reference/pcens_quantile.md),
[`pcens_quantile.default()`](https://primarycensored.epinowcast.org/reference/pcens_quantile.default.md),
[`update.pcens()`](https://primarycensored.epinowcast.org/reference/update.pcens.md)

## Examples

``` r
# Exponential delay with a growing primary event process
pexp_obj <- new_pcens(
  pdist = pexp, dprimary = dexpgrowth,
  primary_args = list(r = 0.3), rate = 0.5
)
pcens_cdf(pexp_obj, q = c(1, 4, 8), pwindow = 2)
#> [1] 0.08649693 0.75597288 0.96697452

# Gamma delay with a declining process
pgamma_obj <- new_pcens(
  pdist = pgamma, dprimary = dexpgrowth,
  primary_args = list(r = -0.3), shape = 2, rate = 1
)
pcens_cdf(pgamma_obj, q = c(1, 4, 8), pwindow = 2)
#> [1] 0.06357351 0.79940359 0.99243183

# Normal delay, for example a difference of event times
pnorm_obj <- new_pcens(
  pdist = pnorm, dprimary = dexpgrowth,
  primary_args = list(r = 0.1), mean = 3, sd = 2
)
pcens_cdf(pnorm_obj, q = c(-1, 3, 8), pwindow = 2)
#> [1] 0.007780575 0.309880218 0.971725556

# Weibull delay
pweibull_obj <- new_pcens(
  pdist = pweibull, dprimary = dexpgrowth,
  primary_args = list(r = 0.2), shape = 1.5, scale = 5
)
pcens_cdf(pweibull_obj, q = c(1, 4, 8), pwindow = 2)
#> [1] 0.01500387 0.36056251 0.80146095
```
