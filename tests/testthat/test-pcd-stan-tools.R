test_that("pcd_stan_path returns correct path", {
  path <- pcd_stan_path()
  expect_true(file.exists(path))
  expect_true(grepl("stan$", path)) # nolint
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
    "expgrowth_rng", "primary_censored_integrand", "dist_lcdf",
    "primary_censored_dist_lcdf", "primary_censored_dist_lpmf"
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
    "expgrowth_pdf", "expgrowth_lpdf", "primary_censored_integrand"
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
  expect_false(grepl("expgrowth_cdf", stan_code, fixed = TRUE))
  expect_false(grepl("dist_lcdf", stan_code, fixed = TRUE))
})

test_that("pcd_load_stan_functions loads all functions as string", {
  all_functions <- pcd_stan_functions()
  stan_code <- pcd_load_stan_functions(functions = all_functions)
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
  pcd_load_stan_functions(
    write_to_file = TRUE,
    output_file = output_file
  )
  expect_true(file.exists(output_file))
  file_content <- readLines(output_file)
  expect_gt(length(file_content), 0)
  expect_true(
    grepl("Stan functions from primarycensoreddist version", file_content[1],
      fixed = TRUE
    )
  )
  unlink(output_file)
})

test_that("pcd_load_stan_functions handles non-existent functions gracefully", {
  non_existent_function <- "non_existent_function"
  stan_code <- pcd_load_stan_functions(functions = non_existent_function)
  expect_identical(stan_code, "")
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
  expect_false(grepl("primary_censored_integrand", stan_code, fixed = TRUE))

  primary_censored_functions <- c(
    "primary_censored_integrand", "dist_lcdf", "primary_censored_dist_lcdf"
  )
  stan_code <- pcd_load_stan_functions(functions = primary_censored_functions)
  for (func in primary_censored_functions) {
    expect_true(
      grepl(func, stan_code),
      info = paste("Function", func, "not found in loaded Stan code")
    )
  }
  expect_false(grepl("expgrowth_pdf", stan_code, fixed = TRUE))
})
