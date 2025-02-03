test_that("pcd_distributions returns expected structure", {
  dists <- pcd_distributions()

  expect_s3_class(dists, "data.frame")
  expect_named(dists, c("name", "params", "aliases", "stan_id"))
  expect_true(all(c("lnorm", "gamma", "weibull", "exp") %in% dists$name))
  expect_gte(nrow(dists), 4)
})

test_that("pcd_primary_distributions returns expected structure", {
  prim_dists <- pcd_primary_distributions()

  expect_s3_class(prim_dists, "data.frame")
  expect_named(prim_dists, c("name", "params", "aliases", "stan_id"))
  expect_true(all(c("unif", "expgrowth") %in% prim_dists$name))
})

test_that("pcd_dist_id works for valid distributions", {
  # Test delay distributions
  expect_identical(pcd_dist_id("lnorm", "delay"), 1L)
  expect_identical(pcd_dist_id("lognormal", "delay"), 1L)
  expect_identical(pcd_dist_id("gamma", "delay"), 2L)

  # Test primary distributions
  expect_identical(pcd_dist_id("unif", "primary"), 1L)
  expect_identical(pcd_dist_id("uniform", "primary"), 1L)
  expect_identical(pcd_dist_id("expgrowth", "primary"), 2L)
})

test_that("pcd_dist_id gives informative errors", {
  expect_error(
    pcd_dist_id("nonexistent", "delay"),
    "No delay distribution found matching: nonexistent."
  )

  expect_error(
    pcd_dist_id("badprim", "primary"),
    "No primary distribution found matching: badprim."
  )
})

test_that("Distribution IDs match Stan model definitions", {
  # This ensures consistency with the Stan code's dist_id numbering
  delay_dists <- pcd_distributions()
  expect_identical(delay_dists$stan_id, seq_len(nrow(delay_dists)))

  prim_dists <- pcd_primary_distributions()
  expect_identical(prim_dists$stan_id, seq_len(nrow(prim_dists)))
})

test_that("name suggestion helper works", {
  # Delay distributions
  expect_match(
    .suggest_dist_name("lnorm", "delay"),
    "Did you mean: lnorm, lognormal"
  )
  expect_match(
    .suggest_dist_name("gama", "delay"),
    "Did you mean: gamma"
  )

  # Primary distributions
  expect_match(
    .suggest_dist_name("expgrowth", "primary"),
    "Did you mean: expgrowth, exponential_growth"
  )
  expect_match(
    .suggest_dist_name("uniff", "primary"),
    "Did you mean: unif, uniform"
  )
})
