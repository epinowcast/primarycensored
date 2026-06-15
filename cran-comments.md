## Submission

This is a patch release fixing a performance regression in the Stan
likelihood that we identified after the 1.5.0 release. A change in 1.5.0
caused positive-support delay distributions to evaluate a redundant
truncation-normalisation term on every likelihood call, slowing the Stan
model by around 30%. 1.5.1 removes the redundant computation. No
user-facing API has changed and results are unchanged. This explains the
short interval since the 1.5.0 release.

## R CMD check results

0 errors | 0 warnings | 1 note

## Reverse dependencies

Checked 1 reverse dependency (EpiNow2) using
r-devel/recheck. No issues found.

## Comments

- NOTE about cmdstanr availability: cmdstanr is listed in Suggests and is
  available from the stan-dev r-universe repository specified in
  Additional_repositories. This is an optional dependency for Stan-based
  functionality and the package works without it.
