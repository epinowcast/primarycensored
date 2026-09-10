# Fit a distribution to doubly censored data

This function wraps the custom approach for fitting distributions to
doubly censored data using fitdistrplus and primarycensored. It handles
primary censoring (when the primary event time is not known exactly),
secondary censoring (when the secondary event time is
interval-censored), and truncation (when events are only observed within
a delay range \[L, D\]).

## Usage

``` r
fitdistdoublecens(
  censdata,
  distr,
  left = "left",
  right = "right",
  pwindow = "pwindow",
  L = "L",
  D = "D",
  dprimary = dunif,
  primary_args = NULL,
  pprimary = NULL,
  dprimary_args = NULL,
  truncation_check_multiplier = 2,
  prior = NULL,
  hazard_model = c("rw", "re"),
  check = TRUE,
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

  A character string naming the distribution to be fitted. Special
  values `"discretestep"` and `"discretehazard"` select the
  non-parametric step-distribution fitting; see Details.

- left:

  Column name for lower bound of observed values (default: "left").

- right:

  Column name for upper bound of observed values (default: "right").

- pwindow:

  Column name for primary window (default: "pwindow").

- L:

  Column name for minimum delay (lower truncation point). For any finite
  L the distribution is left-truncated at L; use `L = -Inf` for no left
  truncation. This is useful for modelling generation intervals where
  day 0 is excluded, particularly when used in renewal models. (default:
  "L"). If the column is not present in censdata, L = -Inf is assumed.

- D:

  Column name for maximum delay (upper truncation point). If finite, the
  distribution is truncated at D. If set to Inf, no upper truncation is
  applied. (default: "D"). Observations whose secondary censoring
  interval straddles `D` (`left < D <= right`) are accepted: the upper
  endpoint is internally clipped to `D` and the likelihood becomes
  `P(X in [left, min(right, D)] | L <= X <= D)`. This is a no-op for the
  standard parametric case where `right <= D`. Observations with
  `left >= D` are rejected because under truncation at `D` no event with
  latent value `>= D` is observable.

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

- primary_args:

  List of additional arguments to be passed to dprimary (and the
  matching primary CDF). For example, when using `dexpgrowth`, you would
  pass `list(min = 0, max = pwindow, r = 0.2)` to set the minimum,
  maximum, and rate parameters. Replaces the deprecated `dprimary_args`;
  defaults to `NULL`.

- pprimary:

  Optional CDF for the primary event distribution. May be a function or
  a character string naming a primary distribution in
  `pcd_primary_distributions`. Defaults to `NULL`, in which case the
  primary CDF is looked up automatically from the registry using the
  `"name"` attribute of `dprimary`. When both `dprimary` and `pprimary`
  carry a `"name"` attribute (or are base R functions whose name can be
  inferred), the two names must agree on everything other than the
  leading `d`/`p` prefix; mismatches such as `dunif` + `pexpgrowth`
  raise an error. Supplying `pprimary` explicitly is mainly useful when
  using a custom primary distribution whose CDF is not in the registry.

- dprimary_args:

  \[Deprecated\] Use `primary_args` instead.

- truncation_check_multiplier:

  Numeric multiplier to use for checking if the truncation time D is
  appropriate relative to the maximum delay. Set to NULL to skip the
  check. Default is 2.

- prior:

  Optional list of prior settings used by the dist function's
  `fit_penalty` attribute (currently only `"discretehazard"`). Each
  element is itself a list with `mean` and `sd` entries. Defaults are
  used for any component not supplied. See
  [`pdiscretehazard()`](https://primarycensored.epinowcast.org/reference/pdiscretehazard.md)
  for the default values.

- hazard_model:

  One of `"rw"` (default) or `"re"`. Only consulted when
  `distr = "discretehazard"`. `"rw"` selects the random-walk transform
  `logit(h_i) = alpha + sigma * cumsum(eps)`; `"re"` selects the IID
  logit random-effect transform `logit(h_i) = alpha + sigma * eps_i`.
  See Details.

- check:

  Logical; if `TRUE` (the default) `pdist` is validated with
  [`check_pdist()`](https://primarycensored.epinowcast.org/reference/check_pdist.md)
  and `dprimary` with
  [`check_dprimary()`](https://primarycensored.epinowcast.org/reference/check_dprimary.md).
  Neither changes across a fit, so validation runs on the first
  likelihood evaluation only rather than on every one. Set to `FALSE` to
  skip it entirely. For non-parametric distributions, `start` is
  required and determines the number of bins; pass `boundaries` here to
  override the default `0:K` unit-width bins.

- ...:

  Additional arguments to be passed to
  [`fitdistrplus::fitdist()`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.html).

## Value

An object of class "fitdist" as returned by fitdistrplus::fitdist.

## Details

### How distribution functions are resolved

The `distr` argument names a distribution. The function looks up the
density and CDF functions by prepending `d` and `p` to the name (e.g.
`distr = "gamma"` resolves to
[`dgamma()`](https://rdrr.io/r/stats/GammaDist.html) and
[`pgamma()`](https://rdrr.io/r/stats/GammaDist.html)). Custom
distributions can be used as long as the corresponding `d<distr>()` and
`p<distr>()` functions are defined.

### Non-parametric distributions

Two non-parametric distributions are supported. They share a common
fitting machinery: the dist function carries a `vector_param` attribute
(`"pmf"` for
[`pdiscretestep()`](https://primarycensored.epinowcast.org/reference/pdiscretestep.md)/[`ddiscretestep()`](https://primarycensored.epinowcast.org/reference/ddiscretestep.md),
`"hazards"` for
[`pdiscretehazard()`](https://primarycensored.epinowcast.org/reference/pdiscretehazard.md)/[`ddiscretehazard()`](https://primarycensored.epinowcast.org/reference/ddiscretehazard.md))
that drives this function to build a closure mapping flat scalar
parameters into the underlying vector argument.

- `distr = "discretestep"`: free parameters `p1, ..., p_{K-1}` (in
  `[0, 1]`); the last bin probability is `1 - sum(p1, ..., p_{K-1})`.
  See
  [`pdiscretestep()`](https://primarycensored.epinowcast.org/reference/pdiscretestep.md)
  for parameterisation details and the soft simplex penalty applied when
  probabilities are infeasible.

- `distr = "discretehazard"`: free parameters `alpha`, `log_sigma`,
  `eps_1, ..., eps_{K-1}`. The hazard form parameterises the same family
  of step distributions as `"discretestep"`, but its free parameters
  drive either a Gaussian random walk on the logit hazard
  (`hazard_model = "rw"`, the default,
  `logit(h_i) = alpha + sigma * cumsum(eps)`) or an IID logit
  random-effect transform (`hazard_model = "re"`,
  `logit(h_i) = alpha + sigma * eps_i` with `eps_i ~ N(0, 1)`). The
  smoothing of the random walk regularises the recovered PMF against
  over-fitting in sparse data and replaces the simplex constraint with
  an unconstrained optimisation; the random-effect variant models
  hazards as independent draws around `alpha` rather than a smoothed
  trajectory. See
  [`pdiscretehazard()`](https://primarycensored.epinowcast.org/reference/pdiscretehazard.md)
  for full parameterisation details and the MAP-equivalent prior penalty
  applied during fitting; pass `prior` to override the default prior
  settings.

For non-parametric distributions `K` is implied by `length(start)`:
`K = length(start) + 1` for `"discretestep"` and `K = length(start) - 1`
for `"discretehazard"`. `start` is therefore required.

## See also

[`pdiscretestep()`](https://primarycensored.epinowcast.org/reference/pdiscretestep.md)
[`pdiscretehazard()`](https://primarycensored.epinowcast.org/reference/pdiscretehazard.md)

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
#> mean 5.003695 0.07759293
#> sd   1.997836 0.06707002
#> Loglikelihood:  -1401.056   AIC:  2806.112   BIC:  2815.927 
#> Correlation matrix:
#>           mean        sd
#> mean 1.0000000 0.3174891
#> sd   0.3174891 1.0000000
#> 

# \donttest{
# Example with discretestep (non-parametric PMF) distribution
set.seed(42)
true_pmf <- c(0.1, 0.3, 0.4, 0.15, 0.05)
step_samples <- rprimarycensored(
  500, rdiscretestep,
  boundaries = 0:5, pmf = true_pmf,
  pwindow = 1, swindow = 1, D = 6
)
step_data <- data.frame(
  left = step_samples,
  right = step_samples + 1,
  pwindow = rep(1, 500),
  D = rep(6, 500)
)
fit_step <- fitdistdoublecens(
  step_data,
  distr = "discretestep",
  boundaries = 0:5,
  start = as.list(setNames(rep(0.2, 4), paste0("p", 1:4)))
)

# Example with discretehazard (logit-hazard random walk) distribution
fit_haz <- fitdistdoublecens(
  step_data,
  distr = "discretehazard",
  boundaries = 0:5,
  start = c(
    list(alpha = -2, log_sigma = log(0.5)),
    as.list(setNames(rep(0, 4), paste0("eps_", 1:4)))
  )
)
#> Warning: diag(V) had non-positive or NA entries; the non-finite result may be dubious
#> Warning: NaNs produced
# }
```
