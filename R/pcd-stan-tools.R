#' Get the path to the Stan code
#'
#' @return A character string with the path to the Stan code
#'
#' @family stantools
#'
#' @export
pcd_stan_path <- function() {
  system.file("stan", "functions", package = "primarycensored")
}

#' Count the number of unmatched braces in a line
#' @noRd
.unmatched_braces <- function(line) {
  ifelse(
    grepl("{", line, fixed = TRUE),
    length(gregexpr("{", line, fixed = TRUE)), 0
  ) -
    ifelse(
      grepl("}", line, fixed = TRUE),
      length(gregexpr("}", line, fixed = TRUE)), 0
    )
}

#' Extract function names or content from Stan code
#'
#' @param content Character vector containing Stan code
#'
#' @param names_only Logical, if TRUE extract function names, otherwise
#' extract function content.
#'
#' @param functions Optional, character vector of function names to extract
#' content for.
#' @return Character vector of function names or content
#' @keywords internal
.extract_stan_functions <- function(
    content, names_only = FALSE, functions = NULL) {
  def_pattern <- "^(real|vector|matrix|void|int|array\\s*<\\s*(real|vector|matrix|int)\\s*>|tuple\\s*<\\s*.*\\s*>)\\s+" # nolint
  func_pattern <- paste0(
    def_pattern, "(\\w+)\\s*\\("
  )
  func_lines <- grep(func_pattern, content, value = TRUE)
  # remove the func_pattern
  func_lines <- sub(def_pattern, "", func_lines)
  # get the next complete word after the pattern until the first (
  func_names <- sub("\\s*\\(.*$", "", func_lines)
  if (!is.null(functions)) {
    func_names <- intersect(func_names, functions)
  }
  if (names_only) {
    return(func_names)
  } else {
    func_content <- character(0)
    for (func_name in func_names) {
      start_line <- grep(paste0(def_pattern, func_name, "\\("), content)
      if (length(start_line) == 0) next
      end_line <- start_line
      brace_count <- 0
      # Ensure we find the first opening brace
      repeat {
        line <- content[end_line]
        brace_count <- brace_count + .unmatched_braces(line)
        end_line <- end_line + 1
        if (brace_count > 0) break
      }
      # Continue until all braces are closed
      repeat {
        line <- content[end_line]
        brace_count <- brace_count + .unmatched_braces(line)
        if (brace_count == 0) break
        end_line <- end_line + 1
      }
      func_content <- c(
        func_content, paste(content[start_line:end_line], collapse = "\n")
      )
    }
    return(func_content)
  }
}

#' Get Stan function names from Stan files
#'
#' This function reads all Stan files in the specified directory and extracts
#' the names of all functions defined in those files.
#'
#' @param stan_path Character string specifying the path to the directory
#' containing Stan files. Defaults to the Stan path of the primarycensored
#' package.
#'
#' @return A character vector containing unique names of all functions found in
#' the Stan files.
#'
#' @export
#'
#' @family stantools
pcd_stan_functions <- function(
    stan_path = primarycensored::pcd_stan_path()) {
  stan_files <- list.files(
    stan_path,
    pattern = "\\.stan$", full.names = TRUE,
    recursive = TRUE
  )
  functions <- character(0)
  for (file in stan_files) {
    content <- readLines(file)
    functions <- c(
      functions, .extract_stan_functions(content, names_only = TRUE)
    )
  }
  unique(functions)
}

#' Get Stan files containing specified functions
#'
#' This function retrieves Stan files from a specified directory, optionally
#' filtering for files that contain specific functions.
#'
#' @param functions Character vector of function names to search for. If NULL,
#'   all Stan files are returned.
#' @inheritParams pcd_stan_functions
#'
#' @return A character vector of file paths to Stan files.
#'
#' @export
#'
#' @family stantools
pcd_stan_files <- function(
    functions = NULL,
    stan_path = primarycensored::pcd_stan_path()) {
  # List all Stan files in the directory
  all_files <- list.files(
    stan_path,
    pattern = "\\.stan$",
    full.names = TRUE,
    recursive = TRUE
  )

  if (is.null(functions)) {
    return(all_files)
  } else {
    # Initialize an empty vector to store matching files
    matching_files <- character(0)

    for (file in all_files) {
      content <- readLines(file)
      extracted_functions <- .extract_stan_functions(content, names_only = TRUE)

      if (any(functions %in% extracted_functions)) {
        matching_files <- c(matching_files, file)
      }
    }

    # remove the path from the file names
    matching_files <- sub(
      paste0(stan_path, "/"), "", matching_files
    )
    return(matching_files)
  }
}

#' Load Stan functions as a string
#'
#' @param functions Character vector of function names to load. Defaults to all
#' functions.
#'
#' @param stan_path Character string, the path to the Stan code. Defaults to the
#' path to the Stan code in the primarycensored package.
#'
#' @param wrap_in_block Logical, whether to wrap the functions in a
#' `functions{}` block. Default is FALSE.
#'
#' @param write_to_file Logical, whether to write the output to a file. Default
#' is FALSE.
#'
#' @param output_file Character string, the path to write the output file if
#' write_to_file is TRUE. Defaults to "pcd_functions.stan".
#'
#' @return A character string containing the requested Stan functions
#'
#' @family stantools
#'
#' @export
pcd_load_stan_functions <- function(
    functions = NULL, stan_path = primarycensored::pcd_stan_path(),
    wrap_in_block = FALSE, write_to_file = FALSE,
    output_file = "pcd_functions.stan") {
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
          names_only = FALSE,
          functions = func
        )
        all_content <- c(all_content, func_content)
      }
    }
  }

  # Add version comment
  version_comment <- paste(
    "// Stan functions from primarycensored version",
    utils::packageVersion("primarycensored")
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

#' Get distribution stan ID by name
#'
#' @param name String. Distribution name or alias
#' @param type String. "delay" or "primary" corresponding to the type of
#'  distribution to use as the look up. If delay then [pcd_distributions()]
#'  is used, if primary then [pcd_primary_distributions()] is used.
#'
#' @return Numeric distribution ID
#' @export
#' @family stantools
#' @examples
#' pcd_stan_dist_id("lnorm")
#' pcd_stan_dist_id("lognormal")
#' pcd_stan_dist_id("gamma")
#' pcd_stan_dist_id("weibull")
#' pcd_stan_dist_id("exp")
#' pcd_stan_dist_id("unif", type = "primary")
#' pcd_stan_dist_id("expgrowth", type = "primary")
pcd_stan_dist_id <- function(name, type = c("delay", "primary")) {
  type <- match.arg(type)
  df <- switch(type,
    delay = primarycensored::pcd_distributions,
    primary = primarycensored::pcd_primary_distributions
  )

  match_idx <- which(df$name == name | df$aliases == name)

  if (length(match_idx) == 0) {
    stop(
      "No ", type, " distribution found matching: ", name, "\n",
      .suggest_dist_name(name, type)
    )
  }

  df$stan_id[match_idx]
}
