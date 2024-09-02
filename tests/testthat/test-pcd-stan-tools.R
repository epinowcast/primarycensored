test_that("pcd_stan_path returns correct path", {
  path <- pcd_stan_path()
  expect_true(file.exists(path))
  expect_true(grepl("stan$", path))
})

test_that("pcd_stan_functions returns unique function names", {
  functions <- pcd_stan_functions()
  expect_type(functions, "character")
  expect_true(length(functions) > 0)
  expect_equal(length(functions), length(unique(functions)))
})

test_that("pcd_load_stan_functions loads all functions as string", {
  all_functions <- pcd_stan_functions()
  stan_code <- pcd_load_stan_functions(functions = all_functions)
  expect_type(stan_code, "character")
  expect_true(nchar(stan_code) > 0)
  expect_true(grepl("functions \\{", stan_code))
})

test_that("pcd_load_stan_functions writes to file when specified", {
  output_file <- tempfile(fileext = ".stan")
  pcd_load_stan_functions(write_to_file = TRUE, output_file = output_file)
  expect_true(file.exists(output_file))
  file_content <- readLines(output_file)
  expect_true(length(file_content) > 0)
  unlink(output_file)
})
