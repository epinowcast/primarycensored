test_that(".extract_function_name returns 'unknown' for closures", {
  fn <- function(x) x + 1 # no .Call line in body
  expect_identical(
    primarycensored:::.extract_function_name(fn), "unknown"
  )
})

test_that(".lookup_pprimary returns NULL when name attribute is unknown", {
  fn <- function(x) x
  expect_null(primarycensored:::.lookup_pprimary(fn))
})

test_that(".lookup_pprimary returns NULL for unregistered named functions", {
  fn <- primarycensored::add_name_attribute(function(x) x, "ddoesnotexist")
  expect_null(primarycensored:::.lookup_pprimary(fn))
})

test_that(".resolve_pdist errors for non-character non-function input", {
  expect_error(
    primarycensored:::.resolve_pdist(123L),
    "function or a single character string"
  )
  expect_error(
    primarycensored:::.resolve_pdist(c("a", "b")),
    "function or a single character string"
  )
})

test_that(".resolve_pdist errors when no matching function exists", {
  expect_error(
    primarycensored:::.resolve_pdist("definitelynotadistribution"),
    "No distribution found"
  )
})

test_that(".resolve_pdist resolves base R p<name> via stats namespace", {
  fn <- primarycensored:::.resolve_pdist("norm")
  expect_true(is.function(fn))
  # Identical to stats::pnorm (we should not pick up an unrelated user
  # object that happens to share the name).
  expect_identical(body(fn), body(stats::pnorm))
})

test_that(".resolve_pdist falls back to base lookup for distributions ", {
  # `geom` is not in `pcd_distributions` but `pgeom` exists in `stats`,
  # so this path exercises the registry-miss fallback.
  fn <- primarycensored:::.resolve_pdist("geom")
  expect_true(is.function(fn))
  expect_identical(attr(fn, "name"), "pgeom")
})

test_that(".resolve_pdist errors when registry entry has no base R impl", {
  # `gengamma` is registered (so it bypasses the registry-miss path) but
  # has `pdist = NA`, indicating no base R implementation.
  expect_error(
    primarycensored:::.resolve_pdist("gengamma"),
    "no base R implementation"
  )
})

test_that(".resolve_primary_args errors when both new and old are supplied", {
  expect_error(
    primarycensored:::.resolve_primary_args(
      list(), list(), "fn"
    ),
    "Supply only one"
  )
})

test_that(".resolve_primary_args returns the dprimary_args value", {
  expect_warning(
    out <- primarycensored:::.resolve_primary_args(
      NULL, list(r = 0.2), "fn"
    ),
    "deprecated"
  )
  expect_identical(out, list(r = 0.2))
})

test_that(".resolve_pprimary errors on non-scalar character", {
  expect_error(
    primarycensored:::.resolve_pprimary(
      stats::dunif, c("unif", "expgrowth")
    ),
    "function or a single character string"
  )
})

test_that(".resolve_pprimary errors on non-function non-string input", {
  expect_error(
    primarycensored:::.resolve_pprimary(stats::dunif, 42L),
    "function or a single character string"
  )
})

test_that(".resolve_pprimary resolves a registry name to a function", {
  fn <- primarycensored:::.resolve_pprimary(
    primarycensored::dexpgrowth, "expgrowth"
  )
  expect_true(is.function(fn))
  expect_identical(attr(fn, "name"), "pexpgrowth")
})

test_that(".sort_eps_names orders by trailing index, not lexicographically", {
  nm <- c("eps_10", "eps_2", "eps_1")
  expect_identical(
    primarycensored:::.sort_eps_names(nm),
    c("eps_1", "eps_2", "eps_10")
  )
  # Single-element / empty input is a no-op
  expect_identical(
    primarycensored:::.sort_eps_names("eps_1"), "eps_1"
  )
  expect_identical(
    primarycensored:::.sort_eps_names(character(0)), character(0)
  )
})
