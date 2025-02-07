test_that("pcd_dist_name works for valid distributions", {
  # Test delay distributions
  expect_identical(pcd_dist_name("lnorm", "delay"), "plnorm")
  expect_identical(pcd_dist_name("lognormal", "delay"), "plnorm")
  expect_identical(pcd_dist_name("gamma", "delay"), "pgamma")

  # Test primary distributions
  expect_identical(pcd_dist_name("unif", "primary"), "dunif")
  expect_identical(pcd_dist_name("uniform", "primary"), "dunif")
  expect_identical(pcd_dist_name("expgrowth", "primary"), "dexpgrowth")
})

test_that("pcd_dist_name gives informative errors", {
  expect_error(
    pcd_dist_name("nonexistent", "delay"),
    "No delay distribution found matching: nonexistent"
  )

  expect_error(
    pcd_dist_name("badprim", "primary"),
    "No primary distribution found matching: badprim"
  )
})
