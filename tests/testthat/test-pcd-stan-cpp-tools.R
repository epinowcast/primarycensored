cpp_function <- "primarycensored_uniform_lcdf_vectorized"
# nolint next: nonportable_path_linter.
cpp_signature <- paste0("vector\\s+", cpp_function, "\\s*\\([^;{]*\\)")

test_that("pcd_stan_header returns the header", {
  header <- pcd_stan_header()
  expect_true(file.exists(header))
  expect_identical(basename(header), "primarycensored.hpp")
})

test_that("pcd_copy_stan_header copies the header and what it includes", {
  path <- file.path(tempdir(), "pcd_header_copy")
  on.exit(unlink(path, recursive = TRUE))
  files <- pcd_copy_stan_header(path)
  expect_true(all(file.exists(files)))
  expect_true(file.exists(file.path(path, "primarycensored.hpp")))
  includes <- grep(
    "#include \"", readLines(file.path(path, "primarycensored.hpp")),
    fixed = TRUE,
    value = TRUE
  )
  included <- sub('.*#include "([^"]+)".*', "\\1", includes)
  expect_true(all(file.exists(file.path(path, included))))
  # Copying again overwrites rather than failing
  expect_no_error(pcd_copy_stan_header(path))
})

test_that("pcd_load_stan_functions keeps Stan bodies by default", {
  code <- pcd_load_stan_functions(cpp_function)
  expect_match(code, paste0(cpp_signature, "\\s*\\{"))
})

test_that("pcd_load_stan_functions declares C++ functions with cpp = TRUE", {
  code <- pcd_load_stan_functions(cpp_function, cpp = TRUE)
  expect_match(code, paste0(cpp_signature, ";"))
  expect_no_match(code, paste0(cpp_signature, "\\s*\\{"))
})

test_that("pcd_load_stan_functions with cpp = TRUE only changes the C++
  functions", {
  all_stan <- pcd_load_stan_functions()
  all_cpp <- pcd_load_stan_functions(cpp = TRUE)
  # Every function is still there, with the same count
  names_of <- function(code) {
    sort(.extract_stan_functions(strsplit(code, "\n", fixed = TRUE)[[1]],
      names_only = TRUE
    ))
  }
  expect_identical(names_of(all_cpp), names_of(all_stan))
  # The C++ functions are declared exactly as expected
  expect_true(grepl(
    paste0(
      "real primarycensored_uniform_lcdf(data real d, data int dist_id,\n",
      "                                  array[] real params, ",
      "data real pwindow);"
    ),
    all_cpp,
    fixed = TRUE
  ))
  # and only they lose their bodies: the code shrinks by their bodies alone
  bodies <- vapply(
    c(
      "primarycensored_uniform_lcdf",
      "primarycensored_uniform_lcdf_vectorized"
    ),
    function(f) {
      def <- .extract_stan_functions(
        strsplit(all_stan, "\n", fixed = TRUE)[[1]],
        functions = f
      )
      nchar(def) - nchar(sub("\\s*\\{.*$", "", def)) - 1L
    },
    numeric(1)
  )
  expect_identical(nchar(all_stan) - nchar(all_cpp), as.integer(sum(bodies)))
})

test_that("pcd_load_stan_functions with cpp = TRUE leaves out what only the
  C++ functions need", {
  deps_stan <- pcd_stan_function_deps(cpp_function)
  expect_gt(length(deps_stan), 1)
  code <- pcd_load_stan_functions(
    cpp_function,
    dependencies = TRUE, cpp = TRUE
  )
  expect_match(code, paste0(cpp_signature, ";"))
  expect_no_match(
    code, "primarycensored_uniform_terms_lcdf_vectorized\\s*\\("
  )
  # Functions still needed by other requested functions are kept
  code <- pcd_load_stan_functions(
    "primarycensored_sone_pmf_vectorized",
    dependencies = TRUE, cpp = TRUE
  )
  expect_match(code, paste0(cpp_signature, ";"))
  expect_match(
    code, "vector primarycensored_uniform_terms_lcdf_vectorized\\s*\\("
  )
})

test_that("pcd_cmdstan_model compiles with the header and rewritten files
  with cpp = TRUE", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)))
  model <- suppressMessages(suppressWarnings(
    pcd_cmdstan_model(cpp = TRUE)
  ))
  expect_true(file.exists(model$exe_file()))
  expect_false(pcd_stan_path() %in% model$include_paths())
  cpp_path <- model$include_paths()
  expect_match(
    paste(readLines(
      file.path(cpp_path, "primarycensored_analytical_cdf.stan")
    ), collapse = "\n"),
    paste0(cpp_signature, ";")
  )
})
