test_that("pcd_stan_dist_id works for valid distributions", {
  # Test delay distributions
  expect_identical(pcd_stan_dist_id("lnorm", "delay"), 1L)
  expect_identical(pcd_stan_dist_id("lognormal", "delay"), 1L)
  expect_identical(pcd_stan_dist_id("gamma", "delay"), 2L)
  expect_identical(pcd_stan_dist_id("gengamma", "delay"), 5L)
  expect_identical(pcd_stan_dist_id("generalized gamma", "delay"), 5L)

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
  # Identifiers are unique and increasing. Some are reserved for delays
  # added in other releases.
  expect_false(anyDuplicated(delay_dists$stan_id) > 0L)
  expect_false(is.unsorted(delay_dists$stan_id, strictly = TRUE))
  expect_identical(delay_dists$stan_id[1:28], 1:28)

  prim_dists <- pcd_primary_distributions
  expect_identical(prim_dists$stan_id, seq_len(nrow(prim_dists)))
})

test_that("pcd_stan_dist_id returns 26L for discretestep distribution", {
  expect_identical(pcd_stan_dist_id("discretestep"), 26L)
  expect_identical(pcd_stan_dist_id("nonparametric"), 26L)
})

test_that("pcd_stan_dist_id returns 27L and 28L for the two hazard variants", {
  expect_identical(pcd_stan_dist_id("discretehazard_rw"), 27L)
  expect_identical(pcd_stan_dist_id("hazard random walk"), 27L)
  expect_identical(pcd_stan_dist_id("discretehazard_re"), 28L)
  expect_identical(pcd_stan_dist_id("hazard random effect"), 28L)
})

test_that("pcd_stan_dist_id returns 31L for the log-logistic distribution", {
  expect_identical(pcd_stan_dist_id("loglogistic"), 31L)
  expect_identical(pcd_stan_dist_id("log-logistic"), 31L)
})
