# Extract function names or content from Stan code

Extract function names or content from Stan code

## Usage

``` r
.extract_stan_functions(content, names_only = FALSE, functions = NULL)
```

## Arguments

- content:

  Character vector containing Stan code

- names_only:

  Logical, if TRUE extract function names, otherwise extract function
  content.

- functions:

  Optional, character vector of function names to extract content for.

## Value

Character vector of function names or content
