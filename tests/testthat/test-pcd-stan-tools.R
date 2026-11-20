test_that("pcd_stan_path returns correct path", {
  path <- pcd_stan_path()
  expect_true(file.exists(path))
  expect_true(grepl(file.path("stan", "functions$"), path))
  expect_true(dir.exists(path))
})

test_that("pcd_stan_functions returns unique function names", {
  functions <- pcd_stan_functions()
  expect_type(functions, "character")
  expect_gt(length(functions), 0)
  expect_length(functions, length(unique(functions)))

  # Check if specific Stan functions are included
  expected_functions <- c(
    "expgrowth_pdf", "expgrowth_lpdf", "expgrowth_cdf", "expgrowth_lcdf",
    "expgrowth_rng", "primarycensored_ode", "dist_lcdf",
    "primarycensored_lcdf", "primarycensored_lpmf"
  )
  for (func in expected_functions) {
    expect_true(
      func %in% functions,
      info = paste("Function", func, "not found in pcd_stan_functions output")
    )
  }
})

test_that("pcd_load_stan_functions loads specific functions", {
  specific_functions <- c(
    "expgrowth_pdf", "expgrowth_lpdf", "primarycensored_ode"
  )
  stan_code <- pcd_load_stan_functions(functions = specific_functions)
  expect_type(stan_code, "character")
  expect_gt(nchar(stan_code), 0)

  for (func in specific_functions) {
    expect_true(
      grepl(func, stan_code),
      info = paste("Function", func, "not found in loaded Stan code")
    )
  }

  # Check that other functions are not included
  expect_false(grepl("real expgrowth_cdf(", stan_code, fixed = TRUE))
  expect_false(grepl("real dist_lcdf(", stan_code, fixed = TRUE))
})

test_that("pcd_load_stan_functions loads all functions as string", {
  all_functions <- pcd_stan_functions()
  stan_code <- pcd_load_stan_functions(all_functions)
  expect_type(stan_code, "character")
  expect_gt(nchar(stan_code), 0)

  for (func in all_functions) {
    expect_true(
      grepl(func, stan_code),
      info = paste("Function", func, "not found in loaded Stan code")
    )
  }
})

test_that("pcd_load_stan_functions wraps functions in block when specified", {
  stan_code <- pcd_load_stan_functions(wrap_in_block = TRUE)
  expect_true(startsWith(stan_code, "functions {"))
  expect_true(endsWith(stan_code, "}"))
})

test_that("pcd_load_stan_functions writes to file when specified", {
  output_file <- tempfile(fileext = ".stan")
  suppressMessages(pcd_load_stan_functions(
    write_to_file = TRUE,
    output_file = output_file
  ))
  expect_true(file.exists(output_file))
  file_content <- readLines(output_file)
  expect_gt(length(file_content), 0)
  expect_true(
    grepl("Stan functions from primarycensored version", file_content[1],
      fixed = TRUE
    )
  )
  unlink(output_file)
})

test_that("pcd_load_stan_functions handles non-existent functions gracefully", {
  non_existent_function <- "non_existent_function"
  stan_code <- pcd_load_stan_functions(functions = non_existent_function)
  expect_length(stan_code, 1)
})

test_that("pcd_load_stan_functions loads functions from specific files", {
  expgrowth_functions <- c(
    "expgrowth_pdf", "expgrowth_lpdf", "expgrowth_cdf",
    "expgrowth_lcdf", "expgrowth_rng"
  )
  stan_code <- pcd_load_stan_functions(functions = expgrowth_functions)
  for (func in expgrowth_functions) {
    expect_true(
      grepl(func, stan_code),
      info = paste("Function", func, "not found in loaded Stan code")
    )
  }
  expect_false(grepl("primarycensored_ode", stan_code, fixed = TRUE))

  primarycensored_functions <- c(
    "primarycensored_ode", "dist_lcdf", "primarycensored_lcdf"
  )
  stan_code <- pcd_load_stan_functions(functions = primarycensored_functions)
  for (func in primarycensored_functions) {
    expect_true(
      grepl(func, stan_code),
      info = paste("Function", func, "not found in loaded Stan code")
    )
  }
  expect_false(grepl("expgrowth_pdf", stan_code, fixed = TRUE))
})

test_that("pcd_load_stan_functions with dependencies includes all deps", {
  result <- pcd_load_stan_functions(
    "primarycensored_pmf",
    dependencies = TRUE
  )

  # Should contain the requested function
  expect_true(
    grepl("primarycensored_pmf", result, fixed = TRUE),
    info = "Result should contain primarycensored_pmf"
  )

  # Should contain direct dependency
  expect_true(
    grepl("primarycensored_lpmf", result, fixed = TRUE),
    info = "Result should contain primarycensored_lpmf"
  )

  # Should contain transitive dependency
  expect_true(
    grepl("primarycensored_lcdf", result, fixed = TRUE),
    info = "Result should contain primarycensored_lcdf"
  )
})

test_that("pcd_load_stan_functions with dependencies works for ODE", {
  result <- pcd_load_stan_functions(
    "primarycensored_ode",
    dependencies = TRUE
  )

  # Should contain dependencies
  expect_true(grepl("dist_lcdf", result, fixed = TRUE))
  expect_true(grepl("primary_lpdf", result, fixed = TRUE))
  expect_true(grepl("expgrowth_lpdf", result, fixed = TRUE))
})

test_that("pcd_stan_function_deps finds functions passed as arguments", {
  # primarycensored_lcdf passes primarycensored_ode to ode_rk45()
  expect_true(
    "primarycensored_ode" %in% pcd_stan_function_deps("primarycensored_lcdf")
  )
})

test_that("pcd_load_stan_functions dependencies orders correctly", {
  result <- pcd_load_stan_functions(
    "primarycensored_pmf",
    dependencies = TRUE
  )

  # Dependencies should appear before the functions that use them
  pos_lpmf <- regexpr("real primarycensored_lpmf(", result, fixed = TRUE)[1]
  pos_pmf <- regexpr("real primarycensored_pmf(", result, fixed = TRUE)[1]

  expect_lt(pos_lpmf, pos_pmf)
})

test_that("pcd_load_stan_functions with deps handles no-dep function", {
  result <- pcd_load_stan_functions("expgrowth_pdf", dependencies = TRUE)
  expect_type(result, "character")
  expect_true(grepl("expgrowth_pdf", result, fixed = TRUE))
})

test_that("pcd_load_stan_functions with deps errors on missing function", {
  expect_error(
    pcd_load_stan_functions("non_existent_function", dependencies = TRUE),
    "not found"
  )
})

test_that("pcd_stan_function_deps returns dependencies in order", {
  deps <- pcd_stan_function_deps("primarycensored_pmf")

  expect_type(deps, "character")
  expect_gt(length(deps), 1)

  # The requested function should be last

  expect_identical(deps[length(deps)], "primarycensored_pmf")

  # Should include direct dependency
  expect_true("primarycensored_lpmf" %in% deps)
})

test_that("pcd_stan_function_deps handles function with no deps", {
  deps <- pcd_stan_function_deps("expgrowth_pdf")

  expect_type(deps, "character")
  expect_length(deps, 1)
  expect_identical(deps, "expgrowth_pdf")
})

test_that("pcd_stan_function_deps errors on missing function", {
  expect_error(
    pcd_stan_function_deps("non_existent_function"),
    "not found"
  )
})

test_that("pcd_stan_files returns correct files", {
  # Test with no functions specified
  all_files <- pcd_stan_files()
  expect_type(all_files, "character")
  expect_gt(length(all_files), 0)
  expect_true(all(grepl(".stan", all_files, fixed = TRUE)))

  # Test with specific functions
  expgrowth_functions <- c("expgrowth_pdf", "expgrowth_lpdf")
  expgrowth_files <- pcd_stan_files(functions = expgrowth_functions)
  expect_type(expgrowth_files, "character")
  expect_gt(length(expgrowth_files), 0)
  expect_true(all(grepl("expgrowth", expgrowth_files, fixed = TRUE)))

  # Test with functions from different files
  mixed_functions <- c("expgrowth_pdf", "primarycensored_ode")
  mixed_files <- pcd_stan_files(functions = mixed_functions)
  expect_type(mixed_files, "character")
  expect_gt(length(mixed_files), 1)
  expect_true(any(grepl("expgrowth", mixed_files, fixed = TRUE)))
  expect_true(any(grepl("primarycensored", mixed_files, fixed = TRUE)))

  # Test with non-existent function
  non_existent <- pcd_stan_files(functions = "non_existent_function")
  expect_type(non_existent, "character")
  expect_length(non_existent, 0)
})

overload_code <- c(
  "real f(real x) {",
  "  return x;",
  "}",
  "real g(real x) { return x; }",
  "real f(vector x,",
  "       real y) {",
  "  if (y > 0) {",
  "    return sum(x);",
  "  }",
  "  return y;",
  "}"
)

test_that(".stan_function_defs locates each definition", {
  defs <- .stan_function_defs(overload_code)
  expect_identical(defs$name, c("f", "g", "f"))
  expect_identical(defs$start, c(1L, 4L, 5L))
  expect_identical(defs$end, c(3L, 4L, 11L))
  expect_identical(
    defs$signature, c("real f(real x)", "real g(real x)",
                      "real f(vector x, real y)")
  )
})

test_that(".select_stan_overloads keeps the first overload by default", {
  defs <- .stan_function_defs(overload_code)
  expect_identical(.select_stan_overloads(defs), c(TRUE, TRUE, FALSE))
})

test_that(".select_stan_overloads chooses by index or signature", {
  defs <- .stan_function_defs(overload_code)
  expect_identical(
    .select_stan_overloads(defs, list(f = 2)), c(FALSE, TRUE, TRUE)
  )
  expect_identical(
    .select_stan_overloads(defs, list(f = "vector x")), c(FALSE, TRUE, TRUE)
  )
  expect_identical(
    .select_stan_overloads(defs, list(f = 1:2)), c(TRUE, TRUE, TRUE)
  )
  # Functions with one definition are always kept
  expect_identical(
    .select_stan_overloads(defs, list(g = 1)), c(TRUE, TRUE, FALSE)
  )
})

test_that(".select_stan_overloads errors on choices that match nothing", {
  defs <- .stan_function_defs(overload_code)
  expect_error(
    .select_stan_overloads(defs, list(f = 3)),
    "No definition of f matches"
  )
  expect_error(
    .select_stan_overloads(defs, list(f = "matrix x")),
    "No definition of f matches"
  )
  expect_error(
    .select_stan_overloads(defs, list(h = 1)), "Function\\(s\\) not found: h"
  )
})

test_that("pcd_load_stan_functions keeps one primarycensored_lpmf overload", {
  count_defs <- function(code) {
    lengths(regmatches(code, gregexpr("real primarycensored_lpmf(", code,
                                      fixed = TRUE)))
  }
  for (functions in list(NULL, "primarycensored_lpmf")) {
    per_row <- pcd_load_stan_functions(functions)
    expect_identical(count_defs(per_row), 1L)
    expect_true(grepl("primarycensored_lpmf(data int d,", per_row,
                      fixed = TRUE))
    grouped <- pcd_load_stan_functions(
      functions, overloads = list(primarycensored_lpmf = "array[] int n")
    )
    expect_identical(count_defs(grouped), 1L)
    expect_true(grepl("primarycensored_lpmf(data array[] int n,", grouped,
                      fixed = TRUE))
    both <- pcd_load_stan_functions(
      functions, overloads = list(primarycensored_lpmf = 1:2)
    )
    expect_identical(count_defs(both), 2L)
  }
})

test_that("pcd_stan_functions lists overloaded functions once", {
  expect_identical(sum(pcd_stan_functions() == "primarycensored_lpmf"), 1L)
})
