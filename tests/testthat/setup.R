# Set data.table print options for compatibility
options(datatable.print.class = FALSE)
options(datatable.print.keys = FALSE)

if (
  !on_ci() || (on_ci() && Sys.info()["sysname"] == "Linux" && not_on_cran())
) {
  library(cmdstanr)
  stan_functions <- pcd_load_stan_functions(
    wrap_in_block = TRUE,
    write_to_file = TRUE,
    output_file = file.path(tempdir(), "pcd_stan_functions.stan")
  )
  model <- suppressMessages(suppressWarnings(cmdstanr::cmdstan_model(
    file.path(tempdir(), "pcd_stan_functions.stan")
  )))
  model$expose_functions(global = TRUE)
}
