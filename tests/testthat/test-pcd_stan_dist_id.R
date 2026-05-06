test_that("pcd_stan_dist_id works for valid distributions", {
  # Test delay distributions
  expect_identical(pcd_stan_dist_id("lnorm", "delay"), 1L)
  expect_identical(pcd_stan_dist_id("lognormal", "delay"), 1L)
  expect_identical(pcd_stan_dist_id("gamma", "delay"), 2L)

  # Test primary distributions
  expect_identical(pcd_stan_dist_id("unif", "primary"), 1L)
  expect_identical(pcd_stan_dist_id("uniform", "primary"), 1L)
  expect_identical(pcd_stan_dist_id("expgrowth", "primary"), 2L)
})

test_that("pcd_stan_dist_id gives informative errors", {
  expect_error(
    pcd_stan_dist_id("nonexistent", "delay"),
    "No delay distribution found matching: nonexistent"
  )

  expect_error(
    pcd_stan_dist_id("badprim", "primary"),
    "No primary distribution found matching: badprim"
  )
})

test_that("Distribution IDs match Stan model definitions", {
  # This ensures consistency with the Stan code's dist_id numbering
  delay_dists <- pcd_distributions
  expect_identical(delay_dists$stan_id, seq_len(nrow(delay_dists)))

  prim_dists <- pcd_primary_distributions
  expect_identical(prim_dists$stan_id, seq_len(nrow(prim_dists)))
})

test_that("pcd_stan_dist_id returns 26L for discretestep distribution", {
  expect_identical(pcd_stan_dist_id("discretestep"), 26L)
  expect_identical(pcd_stan_dist_id("nonparametric"), 26L)
})

test_that("pcd_stan_dist_id returns 27L for discretehazard distribution", {
  expect_identical(pcd_stan_dist_id("discretehazard"), 27L)
  expect_identical(pcd_stan_dist_id("hazard"), 27L)
})
