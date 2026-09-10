# Build the (dpcens_dist, ppcens_dist) closure pair for fitdistdoublecens

Single path: scalar named parameters are gathered from the closure's
environment, optionally folded into the vector argument named by
`vector_param` via the distribution's `param_transform`, and dispatched
through `.dpcens`/`.ppcens`. Formals on the closures are derived from
the supplied `start` list (parametric) or from a `vector_param`-aware
naming convention (non-parametric).

## Usage

``` r
.build_pcens_closures(
  pdist,
  ddist,
  params,
  dprimary,
  primary_args,
  pprimary = NULL,
  vector_param,
  param_transform = NULL,
  fit_penalty,
  prior,
  N,
  start,
  pdist_extras = list()
)
```
