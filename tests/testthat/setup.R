# Set data.table print options for compatibility
options(datatable.print.class = FALSE)
options(datatable.print.keys = FALSE)

if (
  !on_ci() || not_on_cran()
) {
  if (requireNamespace("cmdstanr", quietly = TRUE)) { # nolint
    if (!is.null(cmdstanr::cmdstan_version())) { # nolint
      library(cmdstanr) # nolint
      temp_path <- file.path(tempdir(), "pcd_stan_functions.stan")
      stan_functions <- pcd_load_stan_functions(
        wrap_in_block = TRUE,
        write_to_file = TRUE,
        output_file = temp_path
      )
      model <- suppressMessages(suppressWarnings(
        cmdstanr::cmdstan_model( # nolint
          temp_path
        )
      ))
      model$expose_functions(global = TRUE)
    }
  }
}
