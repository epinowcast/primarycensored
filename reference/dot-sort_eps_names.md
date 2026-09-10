# Sort eps\_\* parameter names by trailing numeric suffix

Lexicographic sorting of `c("eps_1", ..., "eps_10")` returns `eps_10`
before `eps_2`, which would scramble the random-walk innovations when
\\K \> 10\\. This helper extracts the trailing integer and orders by it.

## Usage

``` r
.sort_eps_names(eps_names)
```
