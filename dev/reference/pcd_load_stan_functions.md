# Load Stan functions as a string

Load Stan functions as a string

## Usage

``` r
pcd_load_stan_functions(
  functions = NULL,
  stan_path = primarycensored::pcd_stan_path(),
  wrap_in_block = FALSE,
  write_to_file = FALSE,
  output_file = "pcd_functions.stan"
)
```

## Arguments

- functions:

  Character vector of function names to load. Defaults to all functions.

- stan_path:

  Character string, the path to the Stan code. Defaults to the path to
  the Stan code in the primarycensored package.

- wrap_in_block:

  Logical, whether to wrap the functions in a `functions{}` block.
  Default is FALSE.

- write_to_file:

  Logical, whether to write the output to a file. Default is FALSE.

- output_file:

  Character string, the path to write the output file if write_to_file
  is TRUE. Defaults to "pcd_functions.stan".

## Value

A character string containing the requested Stan functions

## See also

Tools for working with package Stan functions
[`pcd_stan_dist_id()`](https://primarycensored.epinowcast.org/dev/reference/pcd_stan_dist_id.md),
[`pcd_stan_files()`](https://primarycensored.epinowcast.org/dev/reference/pcd_stan_files.md),
[`pcd_stan_functions()`](https://primarycensored.epinowcast.org/dev/reference/pcd_stan_functions.md),
[`pcd_stan_path()`](https://primarycensored.epinowcast.org/dev/reference/pcd_stan_path.md)
