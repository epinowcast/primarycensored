# Get Stan function names from Stan files

This function reads all Stan files in the specified directory and
extracts the names of all functions defined in those files.

## Usage

``` r
pcd_stan_functions(stan_path = primarycensored::pcd_stan_path())
```

## Arguments

- stan_path:

  Character string specifying the path to the directory containing Stan
  files. Defaults to the Stan path of the primarycensored package.

## Value

A character vector containing unique names of all functions found in the
Stan files.

## See also

Tools for working with package Stan functions
[`pcd_load_stan_functions()`](https://primarycensored.epinowcast.org/dev/reference/pcd_load_stan_functions.md),
[`pcd_stan_dist_id()`](https://primarycensored.epinowcast.org/dev/reference/pcd_stan_dist_id.md),
[`pcd_stan_files()`](https://primarycensored.epinowcast.org/dev/reference/pcd_stan_files.md),
[`pcd_stan_path()`](https://primarycensored.epinowcast.org/dev/reference/pcd_stan_path.md)
