## Submission

This release follows 1.5.2 quickly because 1.5.2 introduced a regression
where `fitdistdoublecens()` no longer accepted `fix.arg`, which is fixed
here. It also makes the main R functions and `fitdistdoublecens()`
substantially faster, adds support for zero-width censoring windows, and
adds `update()` and `pcens_pmf()` for `pcens` objects.

## R CMD check results

0 errors | 0 warnings | 1 note

## Reverse dependencies

Checked 2 reverse dependencies (distspec, EpiNow2) with r-devel/recheck,
comparing against 1.5.2. No regressions. The one NOTE on EpiNow2
(checking compiled code) is raised by both versions and is unrelated to
this release.

## Comments

- NOTE about cmdstanr availability: cmdstanr is listed in Suggests and is
  available from the stan-dev r-universe repository specified in
  Additional_repositories. This is an optional dependency for Stan-based
  functionality and the package works without it.
