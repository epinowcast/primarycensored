# Generate random samples from a primary event censored distribution

This function generates random samples from a primary event censored
distribution. It adjusts the distribution by accounting for the primary
event distribution and potential truncation at a maximum delay (D). The
function allows for custom primary event distributions and delay
distributions.

## Usage

``` r
rprimarycensored(
  n,
  rdist,
  pwindow = 1,
  swindow = 1,
  D = Inf,
  rprimary = stats::runif,
  rprimary_args = list(),
  oversampling_factor = 1.2,
  ...
)

rpcens(
  n,
  rdist,
  pwindow = 1,
  swindow = 1,
  D = Inf,
  rprimary = stats::runif,
  rprimary_args = list(),
  oversampling_factor = 1.2,
  ...
)
```

## Arguments

- n:

  Number of random samples to generate.

- rdist:

  Function to generate random samples from the delay distribution for
  example [`stats::rlnorm()`](https://rdrr.io/r/stats/Lognormal.html)
  for lognormal distribution.

- pwindow:

  Primary event window

- swindow:

  Integer specifying the window size for rounding the delay (default is
  1). If `swindow = 0` then no rounding is applied.

- D:

  Maximum delay (truncation point). If finite, the distribution is
  truncated at D. If set to Inf, no truncation is applied. Defaults to
  Inf.

- rprimary:

  Function to generate random samples from the primary distribution
  (default is [`stats::runif()`](https://rdrr.io/r/stats/Uniform.html)).

- rprimary_args:

  List of additional arguments to be passed to rprimary.

- oversampling_factor:

  Factor by which to oversample the number of samples to account for
  truncation (default is 1.2).

- ...:

  Additional arguments to be passed to the distribution function.

## Value

Vector of random samples from the primary event censored distribution
censored by the secondary event window.

## Details

The mathematical formulation for generating random samples from a
primary event censored distribution is as follows:

1.  Generate primary event times (p) from the specified primary event
    distribution (f_p) with parameters phi, defined between 0 and the
    primary event window (pwindow): \$\$p \sim f_p(\phi), \quad p \in
    \[0, pwindow\]\$\$

2.  Generate delays (d) from the specified delay distribution (f_d) with
    parameters theta: \$\$d \sim f_d(\theta)\$\$

3.  Calculate the total delays (t) by adding the primary event times and
    the delays: \$\$t = p + d\$\$

4.  Apply truncation (i.e. remove any delays that fall outside the
    observation window) to ensure that the delays are within the
    specified range \[0, D\], where D is the maximum observable delay:
    \$\$t\_{truncated} = \\t \mid 0 \leq t \< D\\\$\$

5.  Round the truncated delays to the nearest secondary event window
    (swindow): \$\$t\_{valid} = \lfloor \frac{t\_{truncated}}{swindow}
    \rfloor \times swindow\$\$

The function oversamples to account for potential truncation and
generates additional samples if needed to reach the desired number of
valid samples.

## See also

Primary event censored distribution functions
[`dprimarycensored()`](https://primarycensored.epinowcast.org/dev/reference/dprimarycensored.md),
[`pprimarycensored()`](https://primarycensored.epinowcast.org/dev/reference/pprimarycensored.md),
[`qprimarycensored()`](https://primarycensored.epinowcast.org/dev/reference/qprimarycensored.md)

## Examples

``` r
# Example: Lognormal distribution with uniform primary events
rprimarycensored(10, rlnorm, meanlog = 0, sdlog = 1)
#>  [1] 3 0 2 1 0 0 1 2 1 0

# Example: Lognormal distribution with exponential growth primary events
rprimarycensored(
  10, rlnorm,
  rprimary = rexpgrowth, rprimary_args = list(r = 0.2),
  meanlog = 0, sdlog = 1
)
#>  [1] 0 2 3 0 0 0 1 1 0 1
```
