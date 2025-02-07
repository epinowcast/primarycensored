test_that(".suggest_dist_name works for delay distributions", {
  # Test close match suggestion
  expect_identical(
    .suggest_dist_name("gamm"),
    "Did you mean: gamma?"
  )

  # Test multiple close matches
  expect_identical(
    .suggest_dist_name("exp"),
    "Did you mean: exp?"
  )

  # Test no close matches shows all distributions
  expect_identical(
    .suggest_dist_name("notadist"),
    paste0(
      "Available distributions:",
      toString(unique(primarycensored::pcd_distributions$name))
    )
  )
})

test_that(".suggest_dist_name works for primary distributions", {
  # Test close match suggestion
  expect_identical(
    .suggest_dist_name("unif", type = "primary"),
    "Did you mean: unif?"
  )

  # Test no close matches shows all distributions
  expect_identical(
    .suggest_dist_name("notadist", type = "primary"),
    paste0(
      "Available distributions:",
      toString(unique(primarycensored::pcd_primary_distributions$name))
    )
  )
})
