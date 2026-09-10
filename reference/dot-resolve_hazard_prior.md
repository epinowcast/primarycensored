# Resolve the prior settings for the hazard fit_penalty

Defaults: `alpha ~ N(0, 5)`, `log_sigma ~ N(log(0.5), 1)`. User
overrides are merged in: each component of `prior_settings` may itself
be a list with `mean` and `sd` entries.

## Usage

``` r
.resolve_hazard_prior(prior_settings = NULL)
```
