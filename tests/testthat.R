library(testthat)
library(primarycensoreddist)

test_results <- test_check("primarycensoreddist")

if (any(as.data.frame(test_results)$warning > 0)) {
  stop("tests failed with warnings")
}
