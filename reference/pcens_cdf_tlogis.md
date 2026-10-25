# Methods for delays with a truncated logistic primary

Analytical primary event censored CDFs for exponential, gamma and normal
delay distributions with a truncated logistic primary event window, the
[`dtlogis()`](https://primarycensored.epinowcast.org/reference/tlogis.md)
primary distribution with `location` and `scale`. They honour
`use_numeric`, and use the numerical method of
[`pcens_cdf.default()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.default.md)
when no closed form applies.

## Usage

``` r
# S3 method for class 'pcens_pexp_dtlogis'
pcens_cdf(object, q, pwindow, use_numeric = FALSE)

# S3 method for class 'pcens_pgamma_dtlogis'
pcens_cdf(object, q, pwindow, use_numeric = FALSE)

# S3 method for class 'pcens_pnorm_dtlogis'
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

With the window \\\[0, w\]\\ and the logistic CDF \\L\\ with location
\\m\\ and scale \\s\\, the window density is \\L'(z) / D_L\\ with \\D_L
= L(w) - L(0)\\. For a delay with density \\f\\ and CDF \\F\\, the
primary event censored CDF at \\q\\ is \$\$F_L(q) = F(q - w) +
\frac{1}{D_L} \int\_{q - w}^{q} \\L(q - u) - L(0)\\ f(u) du.\$\$
Expanding \\L\\ as a geometric series writes the integral as sums of the
transforms \\T_f(\xi; \tau) = \int\_{-\infty}^{\tau} e^{\xi u} f(u) du\\
at the tilts \\\xi = \pm n / s\\. The positive tilts are for the part of
the window above the location and the negative tilts for the part below
it, so the window is split at the location when it is inside it. Each
transform depends on one endpoint, so each endpoint is evaluated once
and reused when several `q` share it, as for the integer delays of
[`pcens_pmf()`](https://primarycensored.epinowcast.org/reference/pcens_pmf.md).
Terms that hold \\L(0)\\ are written as differences so nothing cancels
when the window is far from the location.

**Truncation.** Near the location the series converge slowly. Each
series is summed with binomial taper weights (the Euler transform after
\\n_0\\ terms), which bounds the error for each term by \\x^{n_0} ((1 -
x) / 2)^M / (1 + x)\\. The fewest terms, at most 64, with a total error
below \\10^{-10} D_L\\ are used. The numerical method is used where
there is no such rule.

**Admissibility.** The exponential and gamma forms need the largest
positive tilt of the series to be below the rate, which holds for a
location after the window or a large rate. Otherwise, and for a window
that is not a single positive finite number,
[`pcens_cdf.default()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.default.md)
is used. The normal form has no restriction.

**Precision.** The CDF agrees with a reference integral to a relative
difference of about 1e-8 in general. In the far tails the error grows
with the cancellation in the integral, to about 1e-7 for a scale much
larger than the window.

A new delay distribution is supported by methods for the internal tilt
generics in `R/tilt_transform.R` and a `pcens_cdf` method for its class
that calls `.pcens_cdf_tlogis()`.

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
[`pcens_cdf_gumbel`](https://primarycensored.epinowcast.org/reference/pcens_cdf_gumbel.md),
[`pcens_pmf()`](https://primarycensored.epinowcast.org/reference/pcens_pmf.md),
[`pcens_pmf.default()`](https://primarycensored.epinowcast.org/reference/pcens_pmf.default.md),
[`pcens_quantile()`](https://primarycensored.epinowcast.org/reference/pcens_quantile.md),
[`pcens_quantile.default()`](https://primarycensored.epinowcast.org/reference/pcens_quantile.default.md),
[`update.pcens()`](https://primarycensored.epinowcast.org/reference/update.pcens.md)

## Examples

``` r
# Normal delay, for example a difference of event times, with a logistic
# window that saturates inside the primary event window
pnorm_obj <- new_pcens(
  pdist = pnorm, dprimary = dtlogis,
  primary_args = list(location = 0.5, scale = 0.2), mean = 3, sd = 2
)
pcens_cdf(pnorm_obj, q = c(-1, 3, 8), pwindow = 2)
#> [1] 0.01209716 0.39163747 0.98590408

# Gamma delay with the midpoint after the window, where the analytic
# solution applies for any rate
pgamma_obj <- new_pcens(
  pdist = pgamma, dprimary = dtlogis,
  primary_args = list(location = 3, scale = 0.5), shape = 2, rate = 1
)
pcens_cdf(pgamma_obj, q = c(1, 4, 8), pwindow = 2)
#> [1] 0.008255513 0.700264706 0.987951915
```
