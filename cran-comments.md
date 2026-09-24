## Submission

This release follows 1.5.2 (on CRAN since 2026-09-11) quickly because 1.5.2
introduced a regression where `fitdistdoublecens()` no longer accepted
`fix.arg`. This release fixes it, alongside other bug fixes and new features.

It also makes the main R functions and `fitdistdoublecens()` substantially
faster, adds support for zero-width primary and secondary censoring
windows, and adds an `update()` method and a `pcens_pmf()` generic for
`pcens` objects.

The one change in behaviour is that `dprimarycensored()` with
`swindow = 0` now returns the primary event censored density. It
previously returned 0.

## Test environments

- Local macOS (aarch64), R 4.6.1, `R CMD check --as-cran`
- GitHub Actions: Ubuntu (R release and oldrel-1), macOS and Windows
  (R release), and an `--as-cran` check on Ubuntu

## R CMD check results

0 errors | 0 warnings | 1 note

## Reverse dependencies

Checked 2 reverse dependencies (distspec, EpiNow2) with r-devel/recheck,
comparing against 1.5.2. No changes between the two versions. distspec is
OK. EpiNow2 raises one NOTE (checking compiled code) under both versions,
so it is unrelated to this release.

We also ran the tests of epidist (not on CRAN) against both versions, with
the same results.

## Comments

- NOTE about cmdstanr availability: cmdstanr is listed in Suggests and is
  available from the stan-dev r-universe repository specified in
  Additional_repositories. This is an optional dependency for Stan-based
  functionality and the package works without it.
