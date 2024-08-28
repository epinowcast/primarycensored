#' Get the path to the Stan code
#'
#' @return A character string with the path to the Stan code
#' @export
#' @aliases pcd_stan_path
pcd_stan_path <- function() {
  system.file("stan", package = "primarycensoreddist")
}

#' List available Stan functions
#'
#' @inheritParams pcd_load_stan_functions
#' @return A character vector of available Stan function names
#' @export
#' @aliases pcd_stan_functions
pcd_stan_functions <- function(
    stan_path = primarycensoreddist::pcd_stan_path()) {
  stan_files <- list.files(
    stan_path,
    pattern = "\\.stan$", full.names = TRUE,
    recursive = TRUE
  )
  functions <- character(0)
  for (file in stan_files) {
    content <- readLines(file)
    func_lines <- grep(
      "^(real|vector|matrix|void)\\s+\\w+\\s*\\(", content,
      value = TRUE
    )
    functions <- c(
      functions, gsub("^.*?\\s+(\\w+)\\s*\\(.*$", "\\1", func_lines)
    )
  }
  unique(functions)
}

#' Load Stan functions as a string
#'
#' @param functions Character vector of function names to load. Defaults to all
#' functions.
#'
#' @param stan_path Character string, the path to the Stan code. Defaults to the
#' path to the Stan code in the primarycensoreddist package.
#'
#' @param wrap_in_block Logical, whether to wrap the functions in a
#' `functions{}` block. Default is FALSE.
#'
#' @param write_to_file Logical, whether to write the output to a file. Default
#' is FALSE.
#'
#' @param output_file Character string, the path to write the output file if
#' write_to_file is TRUE. Defaults to "pcd_stan_functions.stan".
#'
#' @return A character string containing the requested Stan functions
#' @export
#' @aliases pcd_load_stan_functions
pcd_load_stan_functions <- function(
    functions = NULL, stan_path = primarycensoreddist::pcd_stan_path(),
    wrap_in_block = FALSE, write_to_file = FALSE,
    output_file = "pcd_stan_functions.stan") {
  stan_files <- list.files(
    stan_paths,
    pattern = "\\.stan$", full.names = TRUE,
    recursive = TRUE
  )
  all_content <- character(0)

  for (file in stan_files) {
    content <- readLines(file)
    if (is.null(functions)) {
      all_content <- c(all_content, content)
    } else {
      for (func in functions) {
        start_line <- grep(
          paste0("^(real|vector|matrix|void)\\s+", func, "\\s*\\("), content
        )
        if (length(start_line) > 0) {
          end_line <- which(
            cumsum(
              grepl("^\\s*\\{", content[start_line:length(content)])
            ) ==
              cumsum(
                grepl("^\\s*\\}", content[start_line:length(content)])
              )
          )[1] + start_line - 1
          all_content <- c(all_content, content[start_line:end_line])
        }
      }
    }
  }

  # Add version comment
  version_comment <- paste(
    "// Stan functions from primarycensoreddist version",
    utils::packageVersion("primarycensoreddist")
  )
  all_content <- c(version_comment, all_content)

  if (wrap_in_block) {
    all_content <- c("functions {", all_content, "}")
  }

  result <- paste(all_content, collapse = "\n")

  if (write_to_file) {
    writeLines(result, output_file)
    message("Stan functions written to:", output_file, "\n")
  }

  return(result)
}
