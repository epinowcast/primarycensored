# Set data.table print options for compatibility
options(datatable.print.class = FALSE)
options(datatable.print.keys = FALSE)

if (on_ci() && Sys.info()["sysname"] == "Linux" && not_on_cran()) {
  # we only expose stan functions on linux CI
  # because we only test these functions on linux
  suppressMessages(suppressWarnings())
}
