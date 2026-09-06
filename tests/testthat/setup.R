# Set data.table print options for compatibility
options(datatable.print.class = FALSE)
options(datatable.print.keys = FALSE)

if (!on_ci() || not_on_cran()) {
  # nolint start
  if (requireNamespace("cmdstanr", quietly = TRUE)) {
    if (!is.null(cmdstanr::cmdstan_version())) {
      # nolint end
      library(cmdstanr) # nolint
      temp_path <- file.path(tempdir(), "pcd_stan_functions.stan")
      # Use local source Stan files if available (e.g. during
      # development), otherwise fall back to the installed package
      local_stan <- file.path("inst", "stan", "functions")
      stan_path <- if (dir.exists(local_stan)) {
        local_stan
      } else {
        pcd_stan_path()
      }
      stan_functions <- pcd_load_stan_functions(
        wrap_in_block = TRUE,
        write_to_file = TRUE,
        output_file = temp_path,
        stan_path = stan_path
      )
      model <- suppressMessages(
        suppressWarnings(
          cmdstanr::cmdstan_model(
            # nolint
            temp_path
          )
        )
      )
      model$expose_functions(global = TRUE)
      # The Stan helper `hazards_to_pmf` has the same name as the R-side
      # function but a different contract (Stan takes the full K-vector
      # of hazards with the final entry already pinned to 1; R appends
      # the trailing 1 if missing and validates inputs). Globally
      # exposing the Stan version shadows the R one in tests, so we
      # remove it from the global env to restore the R lookup.
      if (exists("hazards_to_pmf", envir = .GlobalEnv, inherits = FALSE)) {
        rm("hazards_to_pmf", envir = .GlobalEnv)
      }
    }
  }
}
