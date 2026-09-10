# Shared MAP penalty for the logit-hazard families

Both RW and RE variants share priors on `alpha`, `log_sigma`, and
`eps_*`; only how the hazards are built from those scalars differs (see
[`.make_hazard_transform()`](https://primarycensored.epinowcast.org/reference/dot-make_hazard_transform.md)).

## Usage

``` r
.hazard_fit_penalty(par_named, N, prior_settings = NULL)
```
