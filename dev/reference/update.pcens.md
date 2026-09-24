# Update the parameters of a pcens object

Replaces the delay distribution parameters, and optionally the primary
event distribution arguments, of an existing `pcens` object. The delay
and primary distribution functions, the primary event CDF and the class
are kept as they are. The updated object therefore dispatches to the
same
[`pcens_cdf()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.md)
method without looking up names or rebuilding the class. This is cheaper
than calling
[`new_pcens()`](https://primarycensored.epinowcast.org/dev/reference/new_pcens.md)
again when one distribution is evaluated for many parameter sets, for
example posterior draws.

## Usage

``` r
# S3 method for class 'pcens'
update(object, ..., primary_args = NULL)
```

## Arguments

- object:

  A `pcens` object as created by
  [`new_pcens()`](https://primarycensored.epinowcast.org/dev/reference/new_pcens.md).

- ...:

  Named delay distribution parameters. Each one replaces the entry of
  the same name in `object$args`, or is added if not present. Parameters
  that are not given keep their current values.

- primary_args:

  Optional named list of primary event distribution arguments. These are
  merged into `object$primary_args` in the same way as `...` is merged
  into `object$args`. Defaults to `NULL`, which leaves the primary event
  distribution arguments unchanged.

## Value

A `pcens` object with the same class as `object` and updated `args`,
`primary_args` and `dprimary_args` fields. See
[`new_pcens()`](https://primarycensored.epinowcast.org/dev/reference/new_pcens.md)
for the fields of a `pcens` object.

## Details

Parameters are merged rather than replaced as a whole, so
`update(object, scale = 3)` changes `scale` and keeps all other
parameters. Parameters cannot be removed; use
[`new_pcens()`](https://primarycensored.epinowcast.org/dev/reference/new_pcens.md)
for that.

A name in `...` that is not already in `object$args` must be an argument
of `object$pdist` other than its first, unless `pdist` takes `...`.
Otherwise an error is raised. Names in `primary_args` are not checked
against `dprimary`.

## See also

Low level primary event censored distribution objects and methods
[`new_pcens()`](https://primarycensored.epinowcast.org/dev/reference/new_pcens.md),
[`pcens_cdf()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.md),
[`pcens_cdf.default()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.default.md),
[`pcens_cdf.pcens_pdiscretehazard()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_pdiscretehazard.md),
[`pcens_cdf.pcens_pdiscretestep()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_pdiscretestep.md),
[`pcens_cdf.pcens_pgamma_dunif()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_pgamma_dunif.md),
[`pcens_cdf.pcens_pgengamma.orig_dunif()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_pgengamma.orig_dunif.md),
[`pcens_cdf.pcens_pgengamma_dunif()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_pgengamma_dunif.md),
[`pcens_cdf.pcens_plnorm_dunif()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_plnorm_dunif.md),
[`pcens_cdf.pcens_pweibull_dunif()`](https://primarycensored.epinowcast.org/dev/reference/pcens_cdf.pcens_pweibull_dunif.md),
[`pcens_pmf()`](https://primarycensored.epinowcast.org/dev/reference/pcens_pmf.md),
[`pcens_pmf.default()`](https://primarycensored.epinowcast.org/dev/reference/pcens_pmf.default.md),
[`pcens_quantile()`](https://primarycensored.epinowcast.org/dev/reference/pcens_quantile.md),
[`pcens_quantile.default()`](https://primarycensored.epinowcast.org/dev/reference/pcens_quantile.default.md)

## Examples

``` r
obj <- new_pcens(
  pdist = pgamma, dprimary = dunif,
  primary_args = list(min = 0, max = 1),
  shape = 1, scale = 1
)
obj <- update(obj, shape = 2, scale = 3)
pcens_cdf(obj, q = c(1, 5, 10), pwindow = 1)
#> [1] 0.01571917 0.44166025 0.82397787

# Update the primary event distribution arguments
obj <- new_pcens(
  pdist = pgamma, dprimary = dexpgrowth,
  primary_args = list(r = 0.2),
  shape = 2, scale = 3
)
obj <- update(obj, primary_args = list(r = 0.5))
pcens_pmf(obj, x = 0:5, pwindow = 1)
#> [1] 0.01385490 0.07361513 0.11112116 0.12144831 0.11699152 0.10530254
```
