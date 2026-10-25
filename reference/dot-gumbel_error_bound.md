# Estimated relative error of the truncated Gumbel series

The machine precision times the sum of the absolute terms and the size
of the logs added, plus the bound \\s_0^{N + 1} / (N + 1)!\\ on the
truncation times the delay mass.

## Usage

``` r
.gumbel_error_bound(
  log_pos,
  log_neg,
  log_bracket,
  log_mass,
  log_c,
  log_s0,
  n_terms,
  scale_t,
  n
)
```

## Arguments

- log_pos, log_neg:

  Logs of the sums of the positive and negative terms.

- log_bracket:

  Log of the bracket, `log_pos` minus `log_neg`.

- log_mass:

  Log of the delay mass in the window.

- log_c:

  Log of \\c(q) = \exp\\-(q - \mu) / \beta\\\\.

- log_s0:

  Log of \\s_0 = e^{\mu / \beta}\\, that is `mu / beta`.

- n_terms:

  Number of terms, \\N\\.

- scale_t:

  Matrix, shaped as `log_delta`, of the size of the log transforms, for
  the rounding estimate.

- n:

  The term indices, `seq_len(n_terms)`.

## Value

Vector of estimated relative errors.
