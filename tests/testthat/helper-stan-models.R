# Compiled Stan models shared across test files, so each is compiled once
# per test run. Models are compiled with model methods so tests can call
# `log_prob()` and `grad_log_prob()`.
.stan_test_models <- new.env(parent = emptyenv())

# Compiles `code` once under `name`. Keep `name` free of underscores, as
# cmdstanr's expose_functions() wrongly flags them as reserved names.
compile_stan_test_model <- function(name, code) {
  if (is.null(.stan_test_models[[name]])) {
    path <- file.path(tempdir(), paste0(name, ".stan"))
    writeLines(code, path)
    .stan_test_models[[name]] <- suppressMessages(suppressWarnings(
      cmdstanr::cmdstan_model(
        path,
        include_paths = pcd_stan_path(), compile_model_methods = TRUE
      )
    ))
  }
  .stan_test_models[[name]]
}

# A model whose functions block holds all the package Stan functions,
# followed by `code` for the other blocks.
compile_stan_functions_model <- function(name, code) {
  compile_stan_test_model(name, paste0(
    pcd_load_stan_functions(wrap_in_block = TRUE, write_to_file = FALSE),
    "\n", code
  ))
}

# The package model. With `grouped = FALSE` no rows are grouped, so every
# row uses primarycensored_lpmf(), as a reference for the grouped model.
compile_pcens_model <- function(grouped = TRUE) {
  code <- paste(
    readLines(system.file("stan", "pcens_model.stan",
      package = "primarycensored"
    )),
    collapse = "\n"
  )
  if (!grouped) {
    pattern <- "primarycensored_sone_groups\\([^;]*\\)"
    stopifnot(grepl(pattern, code))
    code <- sub(pattern, "rep_array(0, N)", code)
  }
  compile_stan_test_model(
    if (grouped) "pcensgrouped" else "pcensperrow", code
  )
}
