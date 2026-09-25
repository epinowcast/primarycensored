# Build the (dpcens_dist, ppcens_dist) closure pair for fitdistdoublecens

Single path: scalar named parameters are gathered from the closure's
environment, optionally folded into the vector argument named by
`vector_param` via the distribution's `param_transform`, and dispatched
through `.dpcens`/`.ppcens`. Formals on the closures are derived from
the supplied `start` list (parametric) or from a `vector_param`-aware
naming convention (non-parametric). Only the supplied parameters become
closure arguments, so a distribution with redundant parameterisations
(e.g. gamma with `rate` and `scale`) is fitted in whichever one `start`
uses.

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
  fix_names = NULL,
  pdist_extras = list(),
  check_once = function() FALSE
)
```

## Arguments

- fix_names:

  Character vector of parameter names held fixed through the `fix.arg`
  argument of
  [`fitdistrplus::fitdist()`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.html),
  or `NULL`. For parametric distributions these are added to the closure
  arguments so fixed parameters reach `pdist`.

- check_once:

  A function returning `TRUE` the first time it is called and `FALSE`
  afterwards, used to validate `pdist` and `dprimary` on the first
  likelihood evaluation only.
