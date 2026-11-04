cpp_function <- "primarycensored_uniform_lcdf_vectorized"
# nolint next: nonportable_path_linter.
cpp_signature <- paste0("vector\\s+", cpp_function, "\\s*\\([^;{]*\\)")

test_that("pcd_stan_header returns the header", {
  header <- pcd_stan_header()
  expect_true(file.exists(header))
  expect_identical(basename(header), "primarycensored.hpp")
})

test_that("pcd_load_stan_functions declares C++ functions by default", {
  code <- pcd_load_stan_functions(cpp_function)
  expect_match(code, paste0(cpp_signature, ";"))
  expect_no_match(code, paste0(cpp_signature, "\\s*\\{"))
})

test_that("pcd_load_stan_functions keeps Stan bodies with cpp = FALSE", {
  code <- pcd_load_stan_functions(cpp_function, cpp = FALSE)
  expect_match(code, paste0(cpp_signature, "\\s*\\{"))
})

test_that("pcd_load_stan_functions copies the header next to the file when
  it declares C++ functions", {
  path <- file.path(tempdir(), "pcd_header_copy")
  on.exit(unlink(path, recursive = TRUE))
  dir.create(path)
  output_file <- file.path(path, "functions.stan")
  expect_message(
    pcd_load_stan_functions(
      cpp_function,
      write_to_file = TRUE, output_file = output_file
    ),
    "C++ header written to",
    fixed = TRUE
  )
  header <- file.path(path, "primarycensored.hpp")
  expect_true(file.exists(header))
  includes <- grep("#include \"", readLines(header), fixed = TRUE, value = TRUE)
  included <- sub('.*#include "([^"]+)".*', "\\1", includes)
  expect_true(all(file.exists(file.path(path, included))))
  # Writing again overwrites rather than failing
  expect_no_error(suppressMessages(pcd_load_stan_functions(
    cpp_function,
    write_to_file = TRUE, output_file = output_file
  )))
})

test_that("pcd_load_stan_functions does not copy the header when it is not
  needed", {
  for (args in list(
    list(functions = "expgrowth_pdf"),
    list(functions = cpp_function, cpp = FALSE)
  )) {
    path <- file.path(tempdir(), "pcd_no_header")
    dir.create(path)
    output_file <- file.path(path, "functions.stan")
    suppressMessages(do.call(pcd_load_stan_functions, c(
      args,
      list(write_to_file = TRUE, output_file = output_file)
    )))
    expect_identical(list.files(path), "functions.stan")
    unlink(path, recursive = TRUE)
  }
})

test_that(".pcd_declare_cpp_functions cuts a signature over several lines at
  the first brace", {
  code <- paste(
    "// A comment with { a brace",
    "real other(real x) { return x; }",
    "vector primarycensored_uniform_lcdf_vectorized(data int start,",
    "                                               data int n)",
    "{",
    "  if (n > 0) { return rep_vector(0, n); }",
    "  return rep_vector(1, n);",
    "}",
    sep = "\n"
  )
  expect_identical(
    .pcd_declare_cpp_functions(code),
    paste(
      "// A comment with { a brace",
      "real other(real x) { return x; }",
      "vector primarycensored_uniform_lcdf_vectorized(data int start,",
      "                                               data int n);",
      sep = "\n"
    )
  )
})

test_that("pcd_load_stan_functions with cpp = TRUE only changes the C++
  functions", {
  all_stan <- pcd_load_stan_functions(cpp = FALSE)
  all_cpp <- pcd_load_stan_functions()
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

test_that("pcd_cmdstan_model passes the cpp settings to cmdstanr", {
  skip_if_not_installed("cmdstanr")
  local_mocked_bindings(
    cmdstan_model = function(...) list(...),
    .package = "cmdstanr"
  )
  stan_path <- pcd_stan_path()
  args <- pcd_cmdstan_model()
  cpp_path <- file.path(tempdir(), "primarycensored_cpp_stan")
  expect_identical(args$include_paths, cpp_path)
  expect_identical(args$user_header, pcd_stan_header())
  expect_identical(args$dir, cpp_path)
  expect_match(
    paste(readLines(
      file.path(cpp_path, "primarycensored_analytical_cdf.stan")
    ), collapse = "\n"),
    paste0(cpp_signature, ";")
  )
  # Unchanged files are not rewritten, so a compiled model is reused
  mtime <- file.mtime(list.files(cpp_path, full.names = TRUE))
  Sys.sleep(1)
  pcd_cmdstan_model()
  expect_identical(
    file.mtime(list.files(cpp_path, full.names = TRUE)), mtime
  )
  # Other include paths are kept, and user arguments are respected
  extra <- tempdir()
  args <- pcd_cmdstan_model(
    include_paths = c(stan_path, extra),
    user_header = "my.hpp", dir = extra
  )
  expect_identical(args$include_paths, c(cpp_path, extra))
  expect_identical(args$user_header, "my.hpp")
  expect_identical(args$dir, extra)
  # cpp = FALSE leaves everything alone
  args <- pcd_cmdstan_model(cpp = FALSE)
  expect_identical(args$include_paths, stan_path)
  expect_null(args$user_header)
  expect_null(args$dir)
  # Without the package Stan files nothing is changed and no warning is
  # given for paths that do not exist
  missing_path <- file.path(tempdir(), "does-not-exist")
  expect_no_warning(args <- pcd_cmdstan_model(include_paths = missing_path))
  expect_identical(args$include_paths, missing_path)
  expect_null(args$user_header)
})
