# Method for a normal delay with a truncated Gumbel primary

Analytical primary event censored CDF for a normal delay distribution
with a truncated Gumbel primary event window, the
[`dtgumbel()`](https://primarycensored.epinowcast.org/reference/tgumbel.md)
primary distribution. It honours `use_numeric`, and uses a numerical
method where the series below is not accurate.

## Usage

``` r
# S3 method for class 'pcens_pnorm_dtgumbel'
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

Write \\G(z) = \exp\[-\exp\\-(z - \mu) / \beta\\\]\\ and \\D = G(w) -
G(0)\\ for a window of width \\w\\, and let \\F\\ be the delay CDF and
\\T_f(\xi; \tau)\\ the transform of
[pcens_cdf_exptilt](https://primarycensored.epinowcast.org/reference/pcens_cdf_exptilt.md).
With \\c(q) = \exp\\-(q - \mu) / \beta\\\\ the reversed window \\G(q -
u) = \exp\\-c(q) e^{u / \beta}\\\\ is a power series in \\e^{u /
\beta}\\, which gives \$\$ F_G(q) = F(q - w) + \frac{1}{D} \Big\[ (1 -
G(0)) \\F(q) - F(q - w)\\ + \sum\_{n \ge 1} \frac{\\-c(q)\\^n}{n!}
\Delta_n(q) \Big\], \$\$ with \\\Delta_n(q) = T_f(n / \beta; q) - T_f(n
/ \beta; q - w)\\. The transforms at each endpoint are evaluated once at
every tilt and reused when several `q` share an endpoint. Each
\\\Delta_n\\ is taken between the lower or the upper tail transforms,
whichever loses less precision, as for
[pcens_cdf_exptilt](https://primarycensored.epinowcast.org/reference/pcens_cdf_exptilt.md).

The terms alternate and are bounded by \\s_0^n / n!\\ times the delay
mass in the window, where \\s_0 = e^{\mu / \beta}\\. They are summed on
the log scale with the odd and even terms separate. The series is
truncated where the bound is below 1e-18.

**Fallback to the numerical method.** The relative error of the series
is estimated for each `q`. Where it is above 1e-8 the numerical method
is used for that `q`. The numerical method is used for every `q` where
\\\mu / \beta\\ is above `log(15)` or where the window is narrow
relative to \\\beta\\.

**Numerical method.** The window density is a spike of width about
\\\beta\\ when \\\mu / \beta\\ is large, which the integration of
[`pcens_cdf.default()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.default.md)
can miss. So
[`pcens_cdf.default()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.default.md)
integrates this primary in a variable in which the integrand is smooth,
over the range that holds the mass of the window. This applies to any
delay with this primary.

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
[`pcens_cdf_exptilt`](https://primarycensored.epinowcast.org/reference/pcens_cdf_exptilt.md),
[`pcens_cdf_tlogis`](https://primarycensored.epinowcast.org/reference/pcens_cdf_tlogis.md),
[`pcens_pmf()`](https://primarycensored.epinowcast.org/reference/pcens_pmf.md),
[`pcens_pmf.default()`](https://primarycensored.epinowcast.org/reference/pcens_pmf.default.md),
[`pcens_quantile()`](https://primarycensored.epinowcast.org/reference/pcens_quantile.md),
[`pcens_quantile.default()`](https://primarycensored.epinowcast.org/reference/pcens_quantile.default.md),
[`update.pcens()`](https://primarycensored.epinowcast.org/reference/update.pcens.md)

## Examples

``` r
pnorm_obj <- new_pcens(
  pdist = pnorm, dprimary = dtgumbel,
  primary_args = list(mu = 0.5, beta = 0.5), mean = 3, sd = 2
)
pcens_cdf(pnorm_obj, q = c(-1, 3, 8), pwindow = 2)
#> [1] 0.01008243 0.35544404 0.98021388
```
