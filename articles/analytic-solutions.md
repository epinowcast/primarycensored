# Analytic solutions for censored delay distributions

## 1 Introduction

For a non-uniform primary event window each analytic solution for a
continuous delay in the package writes the primary event censored
distribution function as a weighted sum of truncated exponential-moment
transforms of the delay
density^(\[[1](#ref-TransformPaperPlaceholder)\]). The uniform window is
the limit of this sum as the tilt goes to zero, and needs only the delay
CDF and its partial moment. Where no tractable form applies the package
integrates numerically. In R it uses
[`pcens_cdf.default()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.default.md)
and in Stan an ODE solver, or for some pairs a numerical method specific
to the pair.

This vignette gives the definitions, the result for each window and the
transform for each delay family. It then lists which pairs of window and
delay the package solves analytically. The derivations are in the
transform paper^(\[[1](#ref-TransformPaperPlaceholder)\]). The double
interval censored model is in Brand et
al.^(\[[2](#ref-BrandInPrepPrimaryCensored)\]).

## 2 Notation

The notation follows the transform
paper^(\[[1](#ref-TransformPaperPlaceholder)\]).

| Symbol | Meaning | In the package |
|:---|:---|:---|
| \\F\_{\mathrm{IP}}\\, \\f\_{\mathrm{IP}}\\ | CDF and density of the delay | `pdist` and its arguments |
| \\w\\, \\w_P\\ | width of the primary event window | `pwindow` |
| \\z \in \[0, w\]\\ | primary event time, measured from the window start |  |
| \\q\\ | evaluation time, measured from the window start | `q` of [`pcens_cdf()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.md) |
| \\f\_{\mathrm{win}}(z; w)\\ | density of the primary event time in the window | `dprimary` and `primary_args` |
| \\H_w(z)\\ | CDF of the primary event time in the window | `pprimary` |
| \\F\_{\mathrm{cens}}(q; w)\\ | primary event censored CDF | [`pcens_cdf()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.md) |
| \\\rho\\ | exponential tilt of the window | `r` of [`dexpgrowth()`](https://primarycensored.epinowcast.org/reference/expgrowth.md) |

The general results write the width as \\w\\ and the window sections
write it as \\w_P\\. The window sections add the window parameters to
\\f\_{\mathrm{win}}\\, for example \\f\_{\mathrm{win}}(z; w_P, \rho)\\.
The subscript IP stands for incubation period. \\F\_{\mathrm{IP}}\\ may
be any delay. A delay on the non-negative reals has
\\F\_{\mathrm{IP}}(x) = f\_{\mathrm{IP}}(x) = 0\\ for \\x \< 0\\. The
main text of the transform paper takes the delay to be non-negative. The
normal delay has support on the reals. For the normal delay the
integrals over the delay \\u\\ in equations [(3.3)](#eq:kernel),
[(3.4)](#eq:transform), [(3.7)](#eq:partial-moment) and [(3.8)](#eq:kn)
have lower limit \\-\infty\\ instead of 0, as in the normal appendix of
the transform paper. \\F\_{\mathrm{IP}}(q - w)\\ is then not zero for
\\q \< w\\.

## 3 The transform reduction

The primary event censored CDF is

\\ F\_{\mathrm{cens}}(q; w) = \int_0^w F\_{\mathrm{IP}}(q - z)
f\_{\mathrm{win}}(z; w)\\dz. \tag{3.1} \\

Let

\\ H_w(z) = \int_0^z f\_{\mathrm{win}}(s; w)\\ds, \tag{3.2} \\

with \\H_w(z) = 0\\ for \\z \leq 0\\ and \\H_w(z) = 1\\ for \\z \geq
w\\. Exchanging the order of integration in equation [(3.1)](#eq:cens)
gives

\\ F\_{\mathrm{cens}}(q; w) = \int_0^q H_w(q - u)
f\_{\mathrm{IP}}(u)\\du. \tag{3.3} \\

The truncated exponential-moment transform of the delay density is

\\ T_f(\xi; \tau) = \int_0^\tau e^{\xi u} f\_{\mathrm{IP}}(u)\\du,
\tag{3.4} \\

with \\T_f(\xi; \tau) = 0\\ for \\\tau \leq 0\\ when the delay is
non-negative. The delay CDF is \\F\_{\mathrm{IP}}(\tau) = T_f(0;
\tau)\\. Suppose the reversed kernel \\H_w(q - u)\\ has the expansion

\\ H_w(q - u) = \sum_k a_k(q, w) e^{\xi_k u} \tag{3.5} \\

for every \\u\\ in the support of the delay up to \\q\\. Substituting it
into equation [(3.3)](#eq:kernel) gives

\\ F\_{\mathrm{cens}}(q; w) = \sum_k a_k(q, w) T_f(\xi_k; q). \tag{3.6}
\\

The window fixes the coefficients \\a_k(q, w)\\ and the tilts \\\xi_k\\.
The delay fixes whether \\T_f\\ is tractable at those tilts.

The reversed kernel is 1 for \\u \leq q - w\\. This saturated part
contributes \\F\_{\mathrm{IP}}(q - w)\\. The expansion is then needed
only on \\\[q - w, q\]\\, and it may differ between parts of that
interval. Each result below for the exponentially tilted, logistic and
Gumbel windows has the form \\F\_{\mathrm{IP}}(q - w)\\ plus weighted
differences \\T_f(\xi_k; \tau^+) - T_f(\xi_k; \tau^-)\\. The endpoints
\\\tau^\pm\\ are \\q\\ and \\q - w\\, and for the logistic window also a
split point \\u^\star\\.

The uniform window uses the partial moment

\\ M(\tau) = \left. \frac{\partial T_f(\xi; \tau)}{\partial \xi}
\right\|\_{\xi = 0} = \int_0^\tau u f\_{\mathrm{IP}}(u)\\du \tag{3.7} \\

and, for \\n \geq 1\\,

\\ K_n(t) = \int_0^t (t - u)^n f\_{\mathrm{IP}}(u)\\du. \tag{3.8} \\

\\K_1(t)\\ is the integral of \\F\_{\mathrm{IP}}\\ up to \\t\\, and

\\ K_1(t) = t\\T_f(0; t) - M(t). \tag{3.9} \\

## 4 Primary event windows

| Window | R primary | Stan primary id | \\f\_{\mathrm{win}}\\ | Tilts \\\xi_k\\ | Result |
|:---|:---|:---|:---|:---|:---|
| Uniform | `dunif` | 1 | \\1 / w_P\\ | none, uses \\M\\ | [(4.2)](#eq:uniform) |
| Exponential tilt | [`dexpgrowth()`](https://primarycensored.epinowcast.org/reference/expgrowth.md), `r` \\= \rho\\ | 2 | \\\rho e^{\rho z} / (e^{\rho w_P} - 1)\\ | \\0\\ and \\-\rho\\ | [(4.1)](#eq:exptilt) |
| Truncated logistic | `dtlogis()`, `location` \\= m\\, `scale` \\= s\\ | 3 | \\L'(z) / D_L\\ | \\n / s\\ for \\n \geq 0\\, \\-n / s\\ for \\n \geq 1\\ | [(4.3)](#eq:logistic-left) to [(4.5)](#eq:logistic-split) |
| Truncated Gumbel | `dtgumbel()`, `mu` \\= \mu_G\\, `beta` \\= \beta\\ | 4 | \\G'(z) / D_G\\ | \\n / \beta\\ for \\n \geq 0\\ | [(4.6)](#eq:gumbel) |

### 4.1 Exponentially tilted window

The window density is \\f\_{\mathrm{win}}(z; w_P, \rho) = \rho e^{\rho
z} / (e^{\rho w_P} - 1)\\ and the window CDF is

\\ H\_{\rho, w_P}(z) = \frac{e^{\rho z} - 1}{e^{\rho w_P} - 1}, \qquad 0
\leq z \leq w_P. \\

A tilt \\\rho \> 0\\ places more mass late in the window, as for a
growing primary event
process^(\[[3](#ref-Park2024),[4](#ref-Zhao2021BETS)\]). A tilt \\\rho
\< 0\\ places more mass early in the window. Let \\J\_\rho(x) =
T_f(-\rho; x)\\. The CDF is

\\ F\_\rho(q; w_P) = F\_{\mathrm{IP}}(q - w_P) + \frac{e^{\rho q}
\\J\_\rho(q) - J\_\rho(q - w_P)\\ - \\F\_{\mathrm{IP}}(q) -
F\_{\mathrm{IP}}(q - w_P)\\}{e^{\rho w_P} - 1}. \tag{4.1} \\

Collecting the \\F\_{\mathrm{IP}}(q - w_P)\\ terms gives the main result
of the transform paper^(\[[1](#ref-TransformPaperPlaceholder)\]). In the
collected form the numerator loses precision through cancellation.
Equation [(4.1)](#eq:exptilt) keeps \\F\_{\mathrm{IP}}(q - w_P)\\ out of
the numerator, and the package evaluates it on the log scale. Writing
\\e^{\rho q} = e^{\rho w_P} e^{\rho (q - w_P)}\\ makes each term a
function of \\q\\ or of \\q - w_P\\ alone. For an integer `pwindow` of
at least 1, \\q - w_P\\ is an integer delay too.
`primarycensored_analytical_lcdf_vectorized()` in Stan therefore
computes the terms once per integer delay \\t\\ and uses them for the
CDF at \\t\\ and at \\t + w_P\\.

### 4.2 Uniform window as the limit \\\rho \to 0\\

As \\\rho \to 0\\ the window density tends to \\1 / w_P\\ and equation
[(4.1)](#eq:exptilt) tends to

\\ F_0(q; w_P) = F\_{\mathrm{IP}}(q - w_P) + \frac{1}{w_P} \int\_{q -
w_P}^q (q - u) f\_{\mathrm{IP}}(u)\\du = \frac{K_1(q) - K_1(q -
w_P)}{w_P}. \tag{4.2} \\

By equation [(3.9)](#eq:k1) this needs only \\F\_{\mathrm{IP}}\\ and the
partial moment \\M\\. It is closed form whenever \\M\\ is. More families
have a closed form \\M\\ than a closed form tilted transform, for
example the lognormal. Brand et al. give \\M\\ for the gamma, lognormal
and Weibull delays^(\[[2](#ref-BrandInPrepPrimaryCensored)\]). For a
gamma delay with daily windows the resulting probability mass function
matches that of Cori et al.^(\[[5](#ref-cori2013new)\]).

### 4.3 Truncated logistic window

With midpoint \\m\\ and scale \\s \> 0\\ let

\\ L(z) = \frac{1}{1 + e^{-(z - m) / s}}, \qquad D_L = L(w_P) - L(0),
\qquad A(q) = e^{-(q - m) / s}, \qquad u^\star = q - m. \\

The window density is \\f\_{\mathrm{win}}(z; w_P, m, s) = L'(z) / D_L\\
and the window CDF is \\H\_{L, w_P}(z) = \\L(z) - L(0)\\ / D_L\\. The
reversed function \\L(q - u) = 1 / \\1 + A(q) e^{u / s}\\\\ is a
geometric series in \\e^{u / s}\\ for \\u \< u^\star\\ and in \\e^{-u /
s}\\ for \\u \> u^\star\\. For \\m \< 0\\,

\\ F_L(q; w_P) = F\_{\mathrm{IP}}(q - w_P) + \frac{1}{D_L} \Big\[
\sum\_{n = 0}^\infty (-1)^n A(q)^n \\T_f(n / s; q) - T_f(n / s; q -
w_P)\\ - L(0) \\F\_{\mathrm{IP}}(q) - F\_{\mathrm{IP}}(q - w_P)\\
\Big\]. \tag{4.3} \\

For \\m \geq w_P\\,

\\ F_L(q; w_P) = F\_{\mathrm{IP}}(q - w_P) + \frac{1}{D_L} \Big\[
\sum\_{n = 1}^\infty (-1)^{n - 1} A(q)^{-n} \\T_f(-n / s; q) - T_f(-n /
s; q - w_P)\\ - L(0) \\F\_{\mathrm{IP}}(q) - F\_{\mathrm{IP}}(q - w_P)\\
\Big\]. \tag{4.4} \\

For \\0 \leq m \< w_P\\ the integral is split at \\u^\star\\,

\\ \begin{aligned} F_L(q; w_P) = F\_{\mathrm{IP}}(q - w_P) +
\frac{1}{D_L} \Big\[ &\sum\_{n = 0}^\infty (-1)^n A(q)^n \\T_f(n / s;
u^\star) - T_f(n / s; q - w_P)\\ \\ + &\sum\_{n = 1}^\infty (-1)^{n - 1}
A(q)^{-n} \\T_f(-n / s; q) - T_f(-n / s; u^\star)\\ - L(0)
\\F\_{\mathrm{IP}}(q) - F\_{\mathrm{IP}}(q - w_P)\\ \Big\].
\end{aligned} \tag{4.5} \\

### 4.4 Truncated Gumbel window

With location \\\mu_G\\ and scale \\\beta \> 0\\ let

\\ G(z) = \exp\left\[-e^{-(z - \mu_G) / \beta}\right\], \qquad D_G =
G(w_P) - G(0), \qquad c(q) = e^{-(q - \mu_G) / \beta}. \\

The window density is \\f\_{\mathrm{win}}(z; w_P, \mu_G, \beta) = G'(z)
/ D_G\\ and the window CDF is \\H\_{G, w_P}(z) = \\G(z) - G(0)\\ /
D_G\\. The reversed function \\G(q - u) = \exp\\-c(q) e^{u / \beta}\\\\
is a power series in \\e^{u / \beta}\\, which gives

\\ F_G(q; w_P) = F\_{\mathrm{IP}}(q - w_P) + \frac{1}{D_G} \Big\[
\sum\_{n = 0}^\infty \frac{\\-c(q)\\^n}{n!} \\T_f(n / \beta; q) - T_f(n
/ \beta; q - w_P)\\ - G(0) \\F\_{\mathrm{IP}}(q) - F\_{\mathrm{IP}}(q -
w_P)\\ \Big\]. \tag{4.6} \\

The series converges for every \\q\\ and needs no split.

## 5 Delay transforms

The table gives \\T_f\\ and \\M\\ for each delay family with a transform
in the package. \\\gamma(\alpha, x)\\ is the lower incomplete gamma
function and \\P(\alpha, x) = \gamma(\alpha, x) / \Gamma(\alpha)\\.
\\\Phi\\ and \\\phi\\ are the standard normal CDF and density. The
classes are those of the transform
paper^(\[[1](#ref-TransformPaperPlaceholder)\]). A transform is analytic
if it is a closed form in elementary or standard functions such as \\P\\
and \\\Phi\\, semi-analytic if it needs a convergent series of special
functions, and numerical if it needs quadrature. The numerical class is
distinct from the numerical method, which integrates equation
[(3.1)](#eq:cens) directly.

| Family (R; Stan id) | \\T_f(\xi; \tau)\\ | \\M(\tau)\\ | Class | Transform used when |
|:---|:---|:---|:---|:---|
| Exponential, rate \\\lambda\\ (`pexp`; 4) | \\\frac{\lambda}{\lambda - \xi} \\1 - e^{-(\lambda - \xi) \tau}\\\\ | \\P(2, \lambda \tau) / \lambda\\ | analytic | \\\xi \< \lambda\\ |
| Gamma, shape \\\alpha\\, rate \\\lambda\\ (`pgamma`; 2) | \\\left(\frac{\lambda}{\lambda - \xi}\right)^\alpha P(\alpha, (\lambda - \xi) \tau)\\ | \\\frac{\alpha}{\lambda} P(\alpha + 1, \lambda \tau)\\ | analytic | \\\xi \< \lambda\\ |
| Normal, mean \\\mu\\, standard deviation \\\sigma\\ (`pnorm`; 18) | \\e^{\xi \mu + \xi^2 \sigma^2 / 2}\\ \Phi\left(\frac{\tau - \mu - \xi \sigma^2}{\sigma}\right)\\ | \\\mu \Phi(z\_\tau) - \sigma \phi(z\_\tau)\\, \\z\_\tau = \frac{\tau - \mu}{\sigma}\\ | analytic | always |
| Weibull, shape \\k\\, scale \\\lambda\\ (`pweibull`; 3) | \\\sum\_{n \geq 0} \frac{(\xi \lambda)^n}{n!} \gamma\left(1 + \frac{n}{k}, (\tau / \lambda)^k\right)\\ | \\\lambda \gamma\left(1 + \frac{1}{k}, (\tau / \lambda)^k\right)\\ | semi-analytic | \\\|\xi\| \tau \< 400\\ and the checks under Numerical safeguards pass |
| Generalised gamma, \\(a, b, c)\\ (`pgengamma.orig`; 5, and `pgengamma` in R) | \\\sum\_{n \geq 0} \frac{\xi^n}{n!} b^{-n / c} \frac{\gamma(a + n / c, b \tau^c)}{\Gamma(a)}\\ | \\b^{-1 / c} \frac{\gamma(a + 1 / c, b \tau^c)}{\Gamma(a)}\\ | semi-analytic | \\\|\xi\| \tau \< 400\\ and the checks under Numerical safeguards pass |
| Lognormal, \\\mu\\, \\\sigma\\ (`plnorm`; 1) | \\\int\_{-\infty}^{\log \tau} e^{\xi e^x} \frac{1}{\sigma} \phi\left(\frac{x - \mu}{\sigma}\right) dx\\ | \\e^{\mu + \sigma^2 / 2}\\ \Phi\left(\frac{\log \tau - \mu}{\sigma} - \sigma\right)\\ | numerical for \\\xi \< 0\\, semi-analytic for \\\xi \> 0\\ | \\\xi \< 0\\ with \\\log(-\xi) + 2 \log \sigma + \mu \< 690\\, or \\\xi \> 0\\ within the series limits below |

The generalised gamma row uses the density \\c b^a u^{a c - 1} e^{-b
u^c} / \Gamma(a)\\ of the transform paper. The Stacy form of
[`flexsurv::pgengamma.orig()`](http://chjackson.github.io/flexsurv-dev/reference/GenGamma.orig.md)^(\[[6](#ref-stacy1962generalization)\])
has `shape` \\= c\\, `scale` \\= b^{-1 / c}\\ and `k` \\= a\\. The
Weibull is the case \\a = 1\\. The Prentice form of
[`flexsurv::pgengamma()`](http://chjackson.github.io/flexsurv-dev/reference/GenGamma.md)
with \\Q \> 0\\ is mapped in R to the Stacy form with `shape` \\= Q /
\sigma\\, `scale` \\= e^\mu Q^{2 \sigma / Q}\\ and `k` \\= Q^{-2}\\. For
\\Q \leq 0\\ the numerical method is used.

The lognormal transform has no closed form. The package uses
Gauss-Legendre quadrature for \\\xi \< 0\\. For \\\xi \> 0\\ it uses the
series \\\sum\_{j \geq 0} \xi^j e^{j \mu + j^2 \sigma^2 / 2}
\Phi\left(\frac{\log \tau - \mu}{\sigma} - j \sigma\right) / j!\\ of
positive terms. In R this series is used where \\\xi \tau \leq 200\\ or
\\\xi w_P \> 2\\. In Stan it is used where \\\xi \tau \leq 60\\ or \\\xi
w_P \> 2\\. Both use at most 20000 terms. Elsewhere the numerical method
is used for that \\q\\.

Some families have a closed form \\M\\ or \\K_1\\ but no tilted
transform in the package. They have analytic solutions for the uniform
window only. Here \\I_x(\cdot, \cdot)\\ is the regularised incomplete
beta function, \\A\_\tau = (\tau / \lambda)^k\\ and \\r_a\\ is defined
below.

| Family (Stan id) | \\M(\tau)\\ or \\K_1(t)\\ | Needs | Available in |
|:---|:---|:---|:---|
| Log-logistic, shape \\k\\, scale \\\lambda\\ (31) | \\M(\tau) = \tau F\_{\mathrm{IP}}(\tau) r\_{1 / k}(A\_\tau)\\ | \\k \geq 0.01\\ | R and Stan |
| Beta, shapes \\\alpha_1\\, \\\alpha_2\\ (9) | \\M(\tau) = \frac{\alpha_1}{\alpha_1 + \alpha_2} I\_{\min(\tau, 1)}(\alpha_1 + 1, \alpha_2)\\ | non-centrality 0 | R and Stan |
| Chi-square, \\\nu\\ degrees of freedom (13) | the gamma with \\\alpha = \nu / 2\\ and \\\lambda = 1 / 2\\ | non-centrality 0 | R and Stan |
| Inverse gamma, shape \\\alpha\\, scale \\\theta\\ (16) | \\M(\tau) = \frac{\theta}{\alpha - 1} \\1 - P(\alpha - 1, \theta / \tau)\\\\ | \\\alpha \> 1\\ | Stan |
| Inverse chi-square, \\\nu\\ (19) | the inverse gamma with \\\alpha = \nu / 2\\ and \\\theta = 1 / 2\\ | \\\nu \> 2\\ | Stan |
| Scaled inverse chi-square, \\\nu\\ and scale \\\varsigma\\ (22) | the inverse gamma with \\\alpha = \nu / 2\\ and \\\theta = \nu \varsigma^2 / 2\\ | \\\nu \> 2\\ | Stan |
| Pareto, shape \\\alpha\\, minimum \\y_m\\ (21) | \\K_1(t) = y_m \left\\e^{\ell} - 1 + \frac{e^{-(\alpha - 1) \ell} - 1}{\alpha - 1}\right\\\\, \\\ell = \log(t / y_m)\\, for \\t \geq y_m\\ |  | Stan |

For the Pareto with \\\alpha = 1\\ the second term takes its limit
\\-\ell\\. The non-central beta and chi-square use the numerical method.

The log-logistic ratio is

\\ r_a(A\_\tau) = \frac{1 + A\_\tau}{A\_\tau^{a + 1}} \int_0^{A\_\tau}
\frac{y^a}{(1 + y)^2}\\dy. \\

For a delay \\U\\ with CDF \\F\_{\mathrm{IP}}\\ it is \\E\[(U / \tau)^{a
k} \mid U \leq \tau\]\\. The integral is a Gauss hypergeometric
function, which base R does not have, so the package evaluates \\r_a\\
by one of three series in \\A\_\tau\\. In Stan the log-logistic
parameters are ordered `[scale, shape]`.

The step and hazard delays of
[`pdiscretestep()`](https://primarycensored.epinowcast.org/reference/pdiscretestep.md)
and
[`pdiscretehazard()`](https://primarycensored.epinowcast.org/reference/pdiscretehazard.md)
(Stan ids 26, 27 and 28) do not use a transform. Their CDF is piecewise
constant, so equation [(3.1)](#eq:cens) is a finite sum of differences
of the window CDF. In R
[`pcens_cdf.pcens_pdiscretestep()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_pdiscretestep.md)
and
[`pcens_cdf.pcens_pdiscretehazard()`](https://primarycensored.epinowcast.org/reference/pcens_cdf.pcens_pdiscretehazard.md)
compute this exactly for any window with a `pprimary` function.

## 6 Coverage by window and delay

The table lists the pairs of window and delay that the package solves.
Each entry is the class of the quantities the window needs, “exact” for
a finite sum, or “numerical method”. The uniform window needs \\M\\ or
\\K_1\\, and the other windows need \\T_f\\. The log-logistic \\M\\ is
semi-analytic because it is evaluated by series in \\A\_\tau\\. The
logistic and Gumbel windows add a truncated series over the tilts. The R
method for a pair is `pcens_cdf.pcens_p<delay>_d<primary>`, for example
[`pcens_cdf.pcens_pgamma_dexpgrowth()`](https://primarycensored.epinowcast.org/reference/pcens_cdf_exptilt.md).

| Delay | Uniform | Exponential tilt | Truncated logistic | Truncated Gumbel |
|:---|:---|:---|:---|:---|
| Exponential | analytic | analytic if \\\lambda \> -\rho\\ | analytic if \\m \geq w_P\\ or \\\lambda\\ exceeds the largest tilt | numerical method |
| Gamma | analytic | analytic if \\\lambda \> -\rho\\ | analytic if \\m \geq w_P\\ or \\\lambda\\ exceeds the largest tilt | numerical method |
| Normal | analytic | analytic | analytic | analytic |
| Weibull | analytic | semi-analytic | numerical method | numerical method |
| Generalised gamma | analytic | semi-analytic | numerical method | numerical method |
| Lognormal | analytic | numerical for \\\rho \> 0\\, semi-analytic for \\\rho \< 0\\ | numerical method | numerical method |
| Log-logistic | semi-analytic | numerical method | numerical method | numerical method |
| Beta, chi-square | analytic | numerical method | numerical method | numerical method |
| Inverse gamma, inverse chi-square, scaled inverse chi-square, Pareto | analytic in Stan, numerical method in R | numerical method | numerical method | numerical method |
| Step and hazard | exact | exact | exact | exact in R, numerical method in Stan |

The logistic window needs the tilts \\n / s\\ and \\-n / s\\ up to the
number of series terms. For \\m \geq w_P\\ it needs no positive tilt.
The Gumbel window needs the tilts \\n / \beta\\ for \\n = 0, \ldots,
N\\, with \\N\\ at least 2. The normal delay has a transform at every
tilt, and is the only delay the package pairs with the Gumbel window.
The Gumbel series checks under Numerical safeguards still apply to it.
The Weibull, generalised gamma and lognormal delays have transforms, but
no method pairs them with the logistic or Gumbel windows.

In Stan `primarycensored_cdf()` and `primarycensored_lcdf()` choose
between the analytical and the numerical path with
`check_for_analytical_params()`. It takes `pwindow`, as the
admissibility of some windows, such as the truncated logistic, depends
on it. `primarycensored_analytical_lcdf()` rejects parameters that fail
these checks.

| Window | Stan functions | Parameter check |
|:---|:---|:---|
| Uniform | `primarycensored_<dist>_uniform_terms()` and `primarycensored_uniform_lcdf_from_terms()` | `check_for_uniform_terms()` and `check_uniform_terms_params()` |
| Exponential tilt | `primarycensored_exptilt_lcdf()` and `primarycensored_analytical_lcdf_vectorized()`, with `log_tilt_transform_pair()` | `check_for_tilt_transform()`, with `check_for_tilt_transform_at()` for the lognormal |
| Truncated logistic | `primarycensored_tlogis_lcdf()` and `primarycensored_tlogis_lcdf_vectorized()` | `check_for_tlogis_params()` |
| Truncated Gumbel | `primarycensored_gumbel_lcdf()`, and `primarycensored_gumbel_terms()` with `primarycensored_gumbel_lcdf_from_terms()` in `primarycensored_analytical_lcdf_vectorized()` | `check_for_gumbel_params()` |

For the Gumbel window, with any delay, and for the log-logistic delay
with the uniform window the numerical method is specific to the pair. In
R these are `.gumbel_numeric()` and `.loglogistic_numeric_cdf()`. In
Stan they are `primarycensored_gumbel_numeric_lcdf()` and
`loglogistic_numeric_lcdf()`. The numerical method is used in the
following cases.

| Condition | R | Stan |
|:---|:---|:---|
| The transform is not available at the parameters | numerical method for every `q` | numerical method |
| `pwindow` is not a single positive finite number | numerical method for every `q` | not applicable |
| A series or error check fails at one `q` | numerical method for that `q` | numerical method for that `q` |
| Lognormal delay, exponential tilt, fewer than 10 `q` and \\\|\rho\| w_P \leq 1\\ | numerical method for every `q` | not applicable |
| `use_numeric = TRUE` | numerical method for every `q` | not applicable |

## 7 Numerical safeguards

### 7.1 Small tilt

The numerator of equation [(4.1)](#eq:exptilt) loses precision through
cancellation as \\\rho \to 0\\. For a gamma delay the loss grows with
the shape. For \\\|\rho\| w_P \< 10^{-2}\\, or \\10^{-3}\\ for delays on
the reals, and \\\|\rho\| (\|q\| + w_P) \< 0.1\\, the package replaces
equation [(4.1)](#eq:exptilt) with its Taylor expansion in \\\rho\\
about the uniform solution, written in the \\K_n\\ of equation
[(3.8)](#eq:kn). To first order this is \\F_0(q; w_P) + \rho \\K_2(q) -
w_P K_1(q) - K_2(q - w_P) - w_P K_1(q - w_P)\\ / (2 w_P)\\. The package
also keeps the second order term, \\\rho^2 \[2 \\K_3(q) - K_3(q -
w_P)\\ - 3 w_P \\K_2(q) + K_2(q - w_P)\\ + w_P^2 \\K_1(q) - K_1(q -
w_P)\\\] / (12 w_P)\\, so the error is third order, \\O\big((\rho
w_P)^3\big)\\. The relative truncation error is below \\(\|\rho\| w_P)^3
/ 160\\.

This form subtracts terms of the size of \\\|q\| + w_P\\, so its
rounding error grows with the distance from the origin. Equation
[(4.1)](#eq:exptilt) has a rounding error of about \\10^{-15} /
(\|\rho\| w_P)\\ at any \\q\\. The expansion has the smaller error while
\\\|\rho\| (\|q\| + w_P)\\ is below about 0.1, so equation
[(4.1)](#eq:exptilt) is used beyond that. In R
[`.exptilt_cdf_finite()`](https://primarycensored.epinowcast.org/reference/dot-exptilt_cdf_finite.md)
and in Stan `primarycensored_exptilt_form()` choose the form for each
delay.

### 7.2 Small delay

For a delay on the non-negative reals with \\q \< w_P\\ the exponential
tilt CDF is \\\\e^{\rho q} J\_\rho(q) - F\_{\mathrm{IP}}(q)\\ / (e^{\rho
w_P} - 1)\\. For \\\|\rho\| q \< 5 \times 10^{-3}\\ the package replaces
the numerator with its Taylor expansion in \\\rho\\, \\\rho K_1(q) +
\rho^2 K_2(q) / 2 + \rho^3 K_3(q) / 6\\. The relative truncation error
is below \\(\|\rho\| q)^3 / 24\\. The errors of both small tilt forms
are below \\10^{-8}\\ at these limits. For delays on the non-negative
reals the \\K_n\\ are computed from the partial moments \\\int_0^t u^j
f\_{\mathrm{IP}}(u)\\du\\, \\j = 0, \ldots, 3\\, by
[`.moments_from_partial()`](https://primarycensored.epinowcast.org/reference/dot-moments_from_partial.md)
in R and `primarycensored_moments_from_partial()` in Stan. The logistic
window uses \\\\L'(0) K_1(q) + L''(0) K_2(q) / 2\\ / D_L\\ for \\q \<
w_P\\ and \\q \< 10^{-4} s\\.

### 7.3 Tail choice

A transform difference is taken from the lower transforms or from the
upper transforms \\T_f(\xi; \infty) - T_f(\xi; \tau)\\, whichever loses
less precision. The exponential, gamma and normal have upper transforms
for every admissible tilt, the lognormal for \\\xi \leq 0\\ and the
other families only at \\\xi = 0\\.

### 7.4 Weibull and generalised gamma series

The series is used only where \\\|\xi\| \tau \< 400\\ and it converges
within 400 terms. For \\\xi \< 0\\ its terms alternate, and it is used
only where the sum of the absolute terms is below \\10^5\\ times the
transform. R and Stan apply this check in the same way.

These series have no upper tail form, so a second check is made at each
\\q\\ on the numerator of equation [(4.1)](#eq:exptilt). It compares the
size of the terms of the numerator, divided by \\e^{\rho w_P} - 1\\,
with a reference. Where the ratio is above \\10^5\\ that \\q\\ uses the
numerical method. In R the reference is the smaller of the CDF and its
survival function, so the check bounds the relative error in both tails.
In Stan the reference is 1, so the check bounds the absolute error of
the CDF at about \\10^{-11}\\. Upper tail probability masses of about
this size or smaller can therefore lose relative accuracy in Stan.

### 7.5 Log-logistic partial moment

The uniform solution differences \\K_1(t) = t F\_{\mathrm{IP}}(t) \\1 -
r\_{1 / k}(A_t)\\\\ at \\q\\ and \\q - w_P\\. Each \\q\\ is screened
with a bound on the rounding error of this difference. That \\q\\ uses
the numerical method where the bound exceeds \\10^{-8}\\ times the
smaller of the CDF and its survival function, or \\10^{-6}\\ times the
density of \\F\_{\mathrm{cens}}\\. R and Stan use the same rule.

### 7.6 Logistic series

The geometric series converge slowly near \\u^\star\\. They are summed
with a binomial taper of the partial sums (Euler summation) to an error
below \\10^{-10} D_L\\, with at most 64 terms.

### 7.7 Gumbel series

The terms alternate. The sum of their absolute values is up to about
\\e^{c(0)}\\ times the delay mass in the window, so the relative
precision lost grows as \\e^{c(0)}\\. The series is used only for \\c(0)
\leq 15\\ and where its estimated relative error is below \\10^{-8}\\.
The estimate adds the rounding error, from the sum of the absolute
terms, and a bound on the truncation error. Where it is above
\\10^{-8}\\ that \\q\\ uses the numerical method.

## 8 Departures from the paper

1.  The exponential and gamma transforms are finite for every \\\xi\\,
    but the package has no stable form for \\\xi \geq \lambda\\.
2.  The lognormal transform is semi-analytic for \\\xi \> 0\\, where the
    transform paper evaluates it by quadrature.
3.  The transform paper states the logistic regimes with strict
    inequalities. The package includes the boundaries as above.
4.  The logistic and Gumbel series coefficients are rearranged so that
    the \\L(0)\\ and \\G(0)\\ terms do not cancel when \\D_L\\ or
    \\D_G\\ is small. The result is unchanged.

## 9 Extending the package

A new delay needs its transform. In R it needs methods for
[`.pcens_tilt_transform()`](https://primarycensored.epinowcast.org/reference/pcens_tilt.md),
[`.pcens_tilt_available()`](https://primarycensored.epinowcast.org/reference/pcens_tilt.md),
[`.pcens_tilt_lower()`](https://primarycensored.epinowcast.org/reference/pcens_tilt.md)
and
[`.pcens_tilt_moments()`](https://primarycensored.epinowcast.org/reference/pcens_tilt.md).
In Stan it needs an entry in `log_tilt_transform_pair()`,
`primarycensored_tilt_moments()`, `check_for_tilt_transform()` and
`check_for_exptilt()`. A delay whose transform holds only at some `q`
also needs the per-point checks that the lognormal uses,
`.pcens_tilt_fits()` in R and `check_for_tilt_transform_at()` in Stan. A
pair of window and delay then needs its `pcens_cdf` method in R and its
dispatch in Stan.

## 10 Learning more

- For more mathematical background on the analytic solutions see the
  [`vignette("why-it-works")`](https://primarycensored.epinowcast.org/articles/why-it-works.md).
- For a more introductory explanation of the primary event censored
  distribution see the
  [`vignette("primarycensored")`](https://primarycensored.epinowcast.org/articles/primarycensored.md).

## References

1\.

PLACEHOLDER: authors to be confirmed. (n.d.). *A transform-based
marginalisation approach for non-uniform primary-censored
epidemiological delays* \[Unpublished manuscript\].

2\.

Brand, S. P. C., Nemcova, B., Pearson, C. A. B., Jamieson, N., Charniga,
K., Azam, J. M., Howes, A., Ward, T., Meyers, L. A., Funk, S., Overton,
C. E., Park, S. W., & Abbott, S. (2026). *A scalable marginalisation
approach for double interval censored epidemiological delays*. In press.

3\.

Park, S. W., Akhmetzhanov, A. R., Charniga, K., Cori, A., Davies, N. G.,
Dushoff, J., Funk, S., Gostic, K., Grenfell, B., Linton, N. M.,
Lipsitch, M., Lison, A., Overton, C. E., Ward, T., & Abbott, S. (2024).
Estimating epidemiological delay distributions for infectious diseases.
*bioRxiv*. <https://doi.org/10.1101/2024.01.12.24301247>

4\.

Zhao, Q., Ju, N., Bacallado, S., & Shah, R. D. (2021). BETS: The dangers
of selection bias in early analyses of the coronavirus disease
(COVID-19) pandemic. *The Annals of Applied Statistics*, *15*(1),
363–390. <https://doi.org/10.1214/20-AOAS1401>

5\.

Cori, A., Ferguson, N. M., Fraser, C., & Cauchemez, S. (2013). A new
framework and software to estimate time-varying reproduction numbers
during epidemics. *American Journal of Epidemiology*, *178*(9),
1505–1512.

6\.

Stacy, E. W. (1962). A generalization of the gamma distribution. *Ann.
Math. Stat.*, *33*(3), 1187–1192.
<https://doi.org/10.1214/aoms/1177704481>
