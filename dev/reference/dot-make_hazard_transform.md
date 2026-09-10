# Build a hazard-vector transform for a given model

Returns a function that maps a named list of scalar parameters (`alpha`,
`log_sigma`, `eps_1`, ..., `eps_{K-1}`) to a length-\\K\\ hazard vector
with the final entry pinned to 1. `model = "rw"` uses a Gaussian random
walk on the logit hazard; `model = "re"` treats the innovations as IID
logit random effects around the intercept.

## Usage

``` r
.make_hazard_transform(model = c("rw", "re"))
```
