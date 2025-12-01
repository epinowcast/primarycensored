# Fit a distribution to doubly censored data

This function wraps the custom approach for fitting distributions to
doubly censored data using fitdistrplus and primarycensored. It handles
primary censoring (when the primary event time is not known exactly),
secondary censoring (when the secondary event time is
interval-censored), and right truncation (when events are only observed
up to a maximum delay).

## Usage

``` r
fitdistdoublecens(
  censdata,
  distr,
  left = "left",
  right = "right",
  pwindow = "pwindow",
  D = "D",
  dprimary = stats::dunif,
  dprimary_name = lifecycle::deprecated(),
  dprimary_args = list(),
  truncation_check_multiplier = 2,
  ...
)
```

## Arguments

- censdata:

  A data frame with columns 'left' and 'right' representing the lower
  and upper bounds of the censored observations. Unlike
  [`fitdistrplus::fitdistcens()`](https://lbbe-software.github.io/fitdistrplus/reference/fitdistcens.html)
  `NA` is not supported for either the upper or lower bounds.

- distr:

  A character string naming the distribution to be fitted. This should
  be the base name of a distribution with corresponding `d` (density)
  and `p` (cumulative distribution) functions available. For example,
  use `"gamma"` (which will use `dgamma` and `pgamma`), `"lnorm"` (for
  `dlnorm` and `plnorm`), `"weibull"`, `"norm"`, etc. Custom
  distributions can also be used as long as the corresponding
  `d<distr>()` and `p<distr>()` functions are defined and loaded.

- left:

  Column name for lower bound of observed values (default: "left").

- right:

  Column name for upper bound of observed values (default: "right").

- pwindow:

  Column name for primary window (default: "pwindow").

- D:

  Column name for maximum delay (truncation point). If finite, the
  distribution is truncated at D. If set to Inf, no truncation is
  applied. (default: "D").

- dprimary:

  Function to generate the probability density function (PDF) of primary
  event times. This function should take a value `x` and a `pwindow`
  parameter, and return a probability density. It should be normalized
  to integrate to 1 over \[0, pwindow\]. Defaults to a uniform
  distribution over \[0, pwindow\]. Users can provide custom functions
  or use helper functions like `dexpgrowth` for an exponential growth
  distribution. See
  [`pcd_primary_distributions()`](https://primarycensored.epinowcast.org/reference/pcd_primary_distributions.md)
  for examples. The package can identify base R distributions for
  potential analytical solutions. For non-base R functions, users can
  apply
  [`add_name_attribute()`](https://primarycensored.epinowcast.org/reference/add_name_attribute.md)
  to yield properly tagged functions if they wish to leverage analytical
  solutions.

- dprimary_name:

  **\[deprecated\]** this argument will be ignored in future versions;
  use
  [`add_name_attribute()`](https://primarycensored.epinowcast.org/reference/add_name_attribute.md)
  on `dprimary` instead

- dprimary_args:

  List of additional arguments to be passed to dprimary. For example,
  when using `dexpgrowth`, you would pass
  `list(min = 0, max = pwindow, r = 0.2)` to set the minimum, maximum,
  and rate parameters

- truncation_check_multiplier:

  Numeric multiplier to use for checking if the truncation time D is
  appropriate relative to the maximum delay. Set to NULL to skip the
  check. Default is 2.

- ...:

  Additional arguments to be passed to
  [`fitdistrplus::fitdist()`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.html).

## Value

An object of class "fitdist" as returned by fitdistrplus::fitdist.

## Details

### How distribution functions are resolved

The `distr` parameter specifies the base name of the distribution. The
function automatically looks up the corresponding density (`d`) and
cumulative distribution (`p`) functions by prepending these prefixes to
the distribution name. For example:

- `distr = "gamma"` uses
  [`dgamma()`](https://rdrr.io/r/stats/GammaDist.html) and
  [`pgamma()`](https://rdrr.io/r/stats/GammaDist.html)

- `distr = "lnorm"` uses
  [`dlnorm()`](https://rdrr.io/r/stats/Lognormal.html) and
  [`plnorm()`](https://rdrr.io/r/stats/Lognormal.html)

- `distr = "weibull"` uses
  [`dweibull()`](https://rdrr.io/r/stats/Weibull.html) and
  [`pweibull()`](https://rdrr.io/r/stats/Weibull.html)

Any distribution available in base R or loaded packages can be used, as
long as the corresponding `d<distr>` and `p<distr>` functions exist and
follow standard R distribution function conventions (first argument is
`x` for density, `q` for CDF).

### What this function does internally

This function creates custom density and CDF functions that account for
primary censoring, secondary censoring, and truncation using
[`dprimarycensored()`](https://primarycensored.epinowcast.org/reference/dprimarycensored.md)
and
[`pprimarycensored()`](https://primarycensored.epinowcast.org/reference/pprimarycensored.md).
These custom functions are then passed to
[`fitdistrplus::fitdist()`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.html)
for maximum likelihood estimation.

The function handles varying observation windows across observations,
making it suitable for real-world data where truncation times or
censoring windows may differ between observations.

## See also

Modelling wrappers for external fitting packages
[`pcd_as_stan_data()`](https://primarycensored.epinowcast.org/reference/pcd_as_stan_data.md),
[`pcd_cmdstan_model()`](https://primarycensored.epinowcast.org/reference/pcd_cmdstan_model.md)

## Examples

``` r
# Example with normal distribution
set.seed(123)
n <- 1000
true_mean <- 5
true_sd <- 2
pwindow <- 2
swindow <- 2
D <- 10
samples <- rprimarycensored(
  n, rnorm,
  mean = true_mean, sd = true_sd,
  pwindow = pwindow, swindow = swindow, D = D
)

delay_data <- data.frame(
  left = samples,
  right = samples + swindow,
  pwindow = rep(pwindow, n),
  D = rep(D, n)
)

fit_norm <- fitdistdoublecens(
  delay_data,
  distr = "norm",
  start = list(mean = 0, sd = 1)
)

summary(fit_norm)
#> Fitting of the distribution ' pcens_dist ' by maximum likelihood 
#> Parameters : 
#>      estimate Std. Error
#> mean 5.007126 0.07883554
#> sd   2.020160 0.06962184
#> Loglikelihood:  -1398.874   AIC:  2801.747   BIC:  2811.563 
#> Correlation matrix:
#>           mean        sd
#> mean 1.0000000 0.3248076
#> sd   0.3248076 1.0000000
#> 
```
