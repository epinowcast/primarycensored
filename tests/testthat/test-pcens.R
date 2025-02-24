test_that("new_pcens creates object with correct structure", {
  pdist <- pgamma
  dprimary <- dunif
  shape <- 2
  rate <- 1

  obj <- new_pcens(
    pdist,
    dprimary,
    list(),
    shape = shape,
    rate = rate
  )

  expect_s3_class(obj, "pcens_pgamma_dunif")
  expect_identical(body(obj$pdist), body(pgamma))
  expect_identical(formals(obj$pdist), formals(pgamma))
  expect_identical(body(obj$dprimary), body(dunif))
  expect_identical(formals(obj$dprimary), formals(dunif))
  expect_identical(obj$args, list(shape = shape, rate = rate))

  new_obj <- new_pcens(
    pgamma,
    dunif,
    list(),
    shape = shape,
    rate = rate
  )
  expect_identical(obj, new_obj)
})
