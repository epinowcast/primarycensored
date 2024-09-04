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
    grepl("Stan functions from primarycensoreddist version", file_content[1],
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

test_that("pcd_stan_files returns correct files", {
  # Test with no functions specified
  all_files <- pcd_stan_files()
  expect_type(all_files, "character")
  expect_gt(length(all_files), 0)
  expect_true(all(grepl("\\.stan$", all_files)))

  # Test with specific functions
  expgrowth_functions <- c("expgrowth_pdf", "expgrowth_lpdf")
  expgrowth_files <- pcd_stan_files(functions = expgrowth_functions)
  expect_type(expgrowth_files, "character")
  expect_gt(length(expgrowth_files), 0)
  expect_true(all(grepl("expgrowth", expgrowth_files)))

  # Test with functions from different files
  mixed_functions <- c("expgrowth_pdf", "primary_censored_integrand")
  mixed_files <- pcd_stan_files(functions = mixed_functions)
  expect_type(mixed_files, "character")
  expect_gt(length(mixed_files), 1)
  expect_true(any(grepl("expgrowth", mixed_files)))
  expect_true(any(grepl("primary_censored", mixed_files)))

  # Test with non-existent function
  non_existent <- pcd_stan_files(functions = "non_existent_function")
  expect_type(non_existent, "character")
  expect_length(non_existent, 0)
})
