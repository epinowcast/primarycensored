## Submission

This release adds non-parametric delay distributions, both a direct PMF
over fixed bins and a discrete-time hazard parameterisation, and analytical
primary event censored solutions for the generalised gamma delay
distribution. It also fixes a gradient bug in the Stan likelihood, where
the log CDF returned a finite value with a non-finite gradient deep in the
lower tail of a narrow lognormal delay.

All additions are backwards compatible.

## R CMD check results

0 errors | 0 warnings | 1 note

## Reverse dependencies

Checked 2 reverse dependencies (distspec, EpiNow2) with r-devel/recheck,
comparing against 1.5.1. No regressions. The one NOTE on EpiNow2
(checking compiled code) is raised by both versions and is unrelated to
this release.

## Comments

- NOTE about cmdstanr availability: cmdstanr is listed in Suggests and is
  available from the stan-dev r-universe repository specified in
  Additional_repositories. This is an optional dependency for Stan-based
  functionality and the package works without it.
