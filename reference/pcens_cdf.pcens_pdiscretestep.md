# Method for step CDF delay with general primary event distribution

Computes the analytic primary event censored CDF for a
piecewise-constant (step) delay distribution and an arbitrary primary
event distribution whose CDF \\F\_{primary}\\ is available via
`object$pprimary`.

## Usage

``` r
# S3 method for class 'pcens_pdiscretestep'
pcens_cdf(object, q, pwindow, use_numeric = FALSE)
```

## Arguments

- object:

  A `primarycensored` object as created by
  [`new_pcens()`](https://primarycensored.epinowcast.org/reference/new_pcens.md).

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

## Details

The observation CDF is \$\$F\_{obs}(q) = \int_0^{pwindow}
F\_{step}(q-p)\\dF\_{primary}(p)\$\$ Because \\F\_{step}\\ is piecewise
constant, the integral reduces to \$\$F\_{obs}(q) = \sum_k c_k
\\\[F\_{primary}(p^{end}\_k) - F\_{primary}(p^{start}\_k)\]\$\$ where
\\c_k\\ is the constant value of \\F\_{step}\\ on the \\k\\-th
sub-interval of the primary event window induced by the step-function
knots. The partition is exact for any bin widths, so bins may be wider
or narrower than `pwindow`, and for boundaries that start below zero.

Falls back to `pcens_cdf.default` when `use_numeric = TRUE` or when no
primary CDF is available on the object.

## See also

Low level primary event censored distribution objects and methods
[`new_pcens()`](https://primarycensored.epinowcast.org/reference/new_pcens.md),
[`pcens_cdf()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.md),
[`pcens_cdf.default()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.default.md),
[`pcens_cdf.pcens_pdiscretehazard()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_pdiscretehazard.md),
[`pcens_cdf.pcens_pgamma_dunif()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_pgamma_dunif.md),
[`pcens_cdf.pcens_pgengamma.orig_dunif()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_pgengamma.orig_dunif.md),
[`pcens_cdf.pcens_pgengamma_dunif()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_pgengamma_dunif.md),
[`pcens_cdf.pcens_plnorm_dunif()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_plnorm_dunif.md),
[`pcens_cdf.pcens_pweibull_dunif()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_pweibull_dunif.md),
[`pcens_quantile()`](https://primarycensored.epinowcast.org/reference/pcens_quantile.md),
[`pcens_quantile.default()`](https://primarycensored.epinowcast.org/reference/pcens_quantile.default.md)
