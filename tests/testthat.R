library(testthat)
library(primarycensored)

test_results <- test_check("primarycensored")

if (any(as.data.frame(test_results)$warning > 0)) {
  stop("tests failed with warnings")
}
