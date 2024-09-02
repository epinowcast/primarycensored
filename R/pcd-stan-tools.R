#' Get the path to the Stan code
#'
#' @return A character string with the path to the Stan code
#'
#' @family stantools
#'
#' @export
pcd_stan_path <- function() {
  system.file("stan", package = "primarycensoreddist")
}

#' Extract function names or content from Stan code
#'
#' @param content Character vector containing Stan code
#'
#' @param extract_names Logical, if TRUE extract function names, otherwise
#' extract function content
#'
#' @param func_name Optional, function name to extract content for
#'
#' @return Character vector of function names or content
#' @keywords internal
.extract_stan_functions <- function(
    content, extract_names = TRUE, func_name = NULL) {
  func_pattern <- paste0(
    "^(real|vector|matrix|void|int)\\s+",
    "(\\w+)\\s*\\("
  )
  if (extract_names) {
    func_lines <- grep(func_pattern, content, value = TRUE)
    return(gsub(func_pattern, "\\2", func_lines))
  } else {
    start_line <- grep(paste0(func_pattern, ".*", func_name), content)
    if (length(start_line) > 0) {
      end_line <- which(
        cumsum(grepl("^\\s*\\{", content[start_line:length(content)])) ==
          cumsum(grepl("^\\s*\\}", content[start_line:length(content)]))
      )[1] + start_line - 1
      return(content[start_line:end_line])
    }
    return(character(0))
  }
}

#' Get Stan function names from Stan files
#'
#' This function reads all Stan files in the specified directory and extracts
#' the names of all functions defined in those files.
#'
#' @param stan_path Character string specifying the path to the directory
#' containing Stan files. Defaults to the Stan path of the primarycensoreddist
#' package.
#'
#' @return A character vector containing unique names of all functions found in
#' the Stan files.
#'
#' @export
#'
#' @family stantools
#'
#' @examples
#' \dontrun{
#' stan_functions <- pcd_stan_functions()
#' print(stan_functions)
#' }
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
    functions <- c(functions, .extract_stan_functions(content))
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
#'
#' @family stantools
#'
#' @export
pcd_load_stan_functions <- function(
    functions = NULL, stan_path = primarycensoreddist::pcd_stan_path(),
    wrap_in_block = FALSE, write_to_file = FALSE,
    output_file = "pcd_stan_functions.stan") {
  stan_files <- list.files(
    stan_path,
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
        func_content <- .extract_stan_functions(
          content,
          extract_names = FALSE, func_name = func
        )
        all_content <- c(all_content, func_content)
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
    message("Stan functions written to: ", output_file, "\n")
  }

  return(result)
}
