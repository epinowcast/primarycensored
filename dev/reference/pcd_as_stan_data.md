# Prepare data for primarycensored Stan model

This function takes in delay data and prepares it for use with the
primarycensored Stan model.

## Usage

``` r
pcd_as_stan_data(
  data,
  delay = "delay",
  delay_upper = "delay_upper",
  n = "n",
  pwindow = "pwindow",
  start_relative_obs_time = "start_relative_obs_time",
  relative_obs_time = "relative_obs_time",
  dist_id,
  primary_id,
  param_bounds,
  primary_param_bounds,
  priors,
  primary_priors,
  compute_log_lik = FALSE,
  use_reduce_sum = FALSE,
  truncation_check_multiplier = 2,
  dist_options = NULL
)
```

## Arguments

- data:

  A data frame containing the delay data.

- delay:

  Column name for observed delays (default: "delay")

- delay_upper:

  Column name for upper bound of delays (default: "delay_upper")

- n:

  Column name for count of observations (default: "n")

- pwindow:

  Column name for primary window (default: "pwindow")

- start_relative_obs_time:

  Column name for start of relative observation time, used as the lower
  truncation point L. Values may be any finite real number (including
  negatives) or `-Inf` to indicate no lower truncation. If the column is
  not present in data, L = `-Inf` is assumed for all observations.
  (default: "start_relative_obs_time")

- relative_obs_time:

  Column name for relative observation time, used as the upper
  truncation point D. Values may be any finite real number (including
  negatives, paired with a smaller `start_relative_obs_time`) or `+Inf`
  to indicate no upper truncation. (default: "relative_obs_time")

- dist_id:

  Integer identifying the delay distribution: You can use
  [`pcd_stan_dist_id()`](https://primarycensored.epinowcast.org/dev/reference/pcd_stan_dist_id.md)
  to get the dist ID for a distribution or look at the
  [pcd_distributions](https://primarycensored.epinowcast.org/dev/reference/pcd_distributions.md)
  data set.

- primary_id:

  Integer identifying the primary distribution: You can use
  [`pcd_stan_dist_id()`](https://primarycensored.epinowcast.org/dev/reference/pcd_stan_dist_id.md)
  to get the primary dist ID for a distribution (make sure to select the
  "primary" type) or look at the
  [pcd_primary_distributions](https://primarycensored.epinowcast.org/dev/reference/pcd_primary_distributions.md)
  data set.

- param_bounds:

  A list with elements `lower` and `upper`, each a numeric vector
  specifying bounds for the delay distribution parameters.

- primary_param_bounds:

  A list with elements `lower` and `upper`, each a numeric vector
  specifying bounds for the primary distribution parameters.

- priors:

  A list with elements `location` and `scale`, each a numeric vector
  specifying priors for the delay distribution parameters.

- primary_priors:

  A list with elements `location` and `scale`, each a numeric vector
  specifying priors for the primary distribution parameters.

- compute_log_lik:

  Logical; compute log likelihood? (default: FALSE)

- use_reduce_sum:

  Logical; use reduce_sum for performance? (default: FALSE)

- truncation_check_multiplier:

  Numeric multiplier to use for checking if the truncation time D is
  appropriate relative to the maximum delay for each unique D value. Set
  to NULL to skip the check. Default is 2.

- dist_options:

  Optional list carrying the shape of a non-parametric delay. When
  `dist_id` is one of `26` (step CDF, Dirichlet prior on the PMF), `27`
  (step CDF, random walk on the logit hazards) or `28` (step CDF, IID
  logit random effects on the hazards), this list **must** be supplied
  and must contain:

  - `K`: integer, the number of bins.

  - `boundaries`: numeric vector of length `K + 1`.

  Priors and `param_bounds` follow the same channel as the parametric
  path: pass them through `priors` and `param_bounds`. The semantics of
  `priors` depend on `dist_id`:

  - `dist_id = 26`: `priors$scale` is the length-`K` Dirichlet
    concentration vector; `priors$location` is unused.

  - `dist_id = 27` or `28`:
    `priors$location = c(alpha_mean, log_sigma_mean)` and
    `priors$scale = c(alpha_sd, log_sigma_sd)`.

  `param_bounds` is unused for the non-parametric paths; pass
  `list(lower = numeric(0), upper = numeric(0))`. If `priors` is empty
  when a non-parametric `dist_id` is given, sensible defaults are
  applied (`Dirichlet(1, ..., 1)` for `dist_id = 26`; `N(0, 5)` on
  `alpha` and `N(0, 1)` on `log_sigma` for 27 and 28).

## Value

A list containing the data formatted for use with
[`pcd_cmdstan_model()`](https://primarycensored.epinowcast.org/dev/reference/pcd_cmdstan_model.md)

## See also

Modelling wrappers for external fitting packages
[`fitdistdoublecens()`](https://primarycensored.epinowcast.org/dev/reference/fitdistdoublecens.md),
[`pcd_cmdstan_model()`](https://primarycensored.epinowcast.org/dev/reference/pcd_cmdstan_model.md)

## Examples

``` r
data <- data.frame(
  delay = c(1, 2, 3),
  delay_upper = c(2, 3, 4),
  n = c(10, 20, 15),
  pwindow = c(1, 1, 2),
  relative_obs_time = c(10, 10, 10)
)
stan_data <- pcd_as_stan_data(
  data,
  dist_id = 1,
  primary_id = 1,
  param_bounds = list(lower = c(0, 0), upper = c(10, 10)),
  primary_param_bounds = list(lower = numeric(0), upper = numeric(0)),
  priors = list(location = c(1, 1), scale = c(1, 1)),
  primary_priors = list(location = numeric(0), scale = numeric(0))
)
#> The truncation time D (10) is larger than 2 times the maximum observed delay (3). Consider setting D to Inf for better efficiency with minimal accuracy cost for this case.
```
