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

#' Get the path to the C++ header for the Stan functions
#'
#' Some Stan functions have a C++ version with a hand-written gradient,
#' which is faster than the Stan version. Load the Stan functions with
#' `pcd_load_stan_functions(cpp = TRUE)` to declare these functions without
#' a body, and compile with this header so that the C++ versions are used.
#' With cmdstanr, pass it as `user_header`. [pcd_cmdstan_model()] does this
#' by default. Use [pcd_copy_stan_header()] to vendor the header into
#' another package.
#'
#' @return A character string with the path to the header.
#'
#' @family stantools
#'
#' @export
#' @examples
#' pcd_stan_header()
pcd_stan_header <- function() {
  system.file(
    "include", "primarycensored.hpp",
    package = "primarycensored", mustWork = TRUE
  )
}

#' Copy the C++ header for the Stan functions
#'
#' Copies the header returned by [pcd_stan_header()], and the files it
#' includes, into a directory. Use this with
#' `pcd_load_stan_functions(cpp = TRUE)` to vendor the Stan functions and
#' their C++ versions into another package together.
#'
#' @param path Directory to copy the header into. It is created if it does
#'   not exist.
#'
#' @return The paths of the copied files, invisibly. Include
#'   `primarycensored.hpp` from `path` before the model code.
#'
#' @family stantools
#'
#' @export
#' @examples
#' pcd_copy_stan_header(tempdir())
pcd_copy_stan_header <- function(path) {
  include_dir <- dirname(pcd_stan_header())
  files <- list.files(include_dir, recursive = TRUE)
  for (sub_dir in unique(file.path(path, dirname(files)))) {
    dir.create(sub_dir, showWarnings = FALSE, recursive = TRUE)
  }
  copied <- file.copy(
    file.path(include_dir, files), file.path(path, files),
    overwrite = TRUE
  )
  if (!all(copied)) {
    stop("Could not copy: ", toString(files[!copied]), call. = FALSE)
  }
  invisible(file.path(path, files))
}

#' Names of the Stan functions that have a C++ version
#'
#' @return Character vector of function names.
#' @noRd
.pcd_cpp_functions <- function() {
  "primarycensored_uniform_lcdf_vectorized"
}

#' Replace the bodies of functions that have a C++ version with declarations
#'
#' @param content Character string of Stan code.
#'
#' @return `content` with each function in `.pcd_cpp_functions()` declared
#'   without a body.
#' @noRd
.pcd_declare_cpp_functions <- function(content) {
  stan_lines <- strsplit(content, "\n", fixed = TRUE)[[1]]
  for (func in .pcd_cpp_functions()) {
    for (def in .extract_stan_functions(stan_lines, functions = func)) {
      def_start <- sub("\\s*\\{.*$", "", def)
      content <- sub(def, paste0(def_start, ";"), content, fixed = TRUE)
    }
  }
  content
}

#' Find dependencies of a Stan function
#'
#' @param func_content Character string containing the function body
#' @param all_func_names Character vector of all known function names
#'
#' @return Character vector of function names that are called
#' @noRd
.pcd_stan_function_deps <- function(func_content, all_func_names) {
  deps <- character(0)
  for (fname in all_func_names) {
    # Match function_name followed by ( with optional whitespace
    # Use word boundary to avoid partial matches
    pattern <- paste0("\\b", fname, "\\s*\\(")
    if (grepl(pattern, func_content)) {
      deps <- c(deps, fname)
    }
  }
  unique(deps)
}

#' Build dependency graph for all Stan functions
#'
#' @param stan_path Path to Stan functions directory
#'
#' @return Named list where each element is a vector of dependencies
#' @noRd
.pcd_build_dep_graph <- function(stan_path) {
  stan_files <- list.files(
    stan_path,
    pattern = "\\.stan$",
    full.names = TRUE,
    recursive = TRUE
  )

  # First pass: get all function names
  all_func_names <- character(0)
  for (file in stan_files) {
    content <- readLines(file)
    all_func_names <- c(
      all_func_names,
      .extract_stan_functions(content, names_only = TRUE)
    )
  }
  all_func_names <- unique(all_func_names)

  # Second pass: build dependency graph
  dep_graph <- list()
  for (file in stan_files) {
    content <- readLines(file)
    func_names <- .extract_stan_functions(content, names_only = TRUE)

    for (fname in func_names) {
      func_content <- .extract_stan_functions(
        content,
        names_only = FALSE,
        functions = fname
      )
      if (length(func_content) > 0) {
        # Find dependencies, excluding self-references
        deps <- .pcd_stan_function_deps(func_content, all_func_names)
        deps <- setdiff(deps, fname)
        dep_graph[[fname]] <- deps
      }
    }
  }

  dep_graph
}

#' Resolve all dependencies for a function recursively
#'
#' @param func_name Name of the function
#' @param dep_graph Dependency graph from .pcd_build_dep_graph
#' @param resolved Already resolved functions (for recursion)
#' @param visiting Currently visiting functions (for cycle detection)
#'
#' @return Character vector of all dependencies in topological order
#' @noRd
.pcd_resolve_deps <- function(func_name, dep_graph,
                              resolved = character(0),
                              visiting = character(0)) {
  if (func_name %in% resolved) {
    return(resolved)
  }
  if (func_name %in% visiting) {
    # Circular dependency - skip to avoid infinite loop
    return(resolved)
  }

  visiting <- c(visiting, func_name)
  deps <- dep_graph[[func_name]]

  if (length(deps) > 0) {
    for (dep in deps) {
      resolved <- .pcd_resolve_deps(dep, dep_graph, resolved, visiting)
    }
  }

  c(resolved, func_name)
}

#' Get dependencies for a Stan function
#'
#' Returns all Stan functions that the specified function depends on,
#' in topological order (dependencies before the functions that use them).
#'
#' @param function_name Character string, the name of the Stan function.
#' @inheritParams pcd_stan_functions
#'
#' @return A character vector of function names that the specified function
#' depends on, ordered so that dependencies come before functions that use
#' them. The requested function itself is included as the last element.
#'
#' @family stantools
#'
#' @export
#'
#' @examples
#' # See what primarycensored_lpmf depends on
#' pcd_stan_function_deps("primarycensored_lpmf")
#'
#' # A function with no dependencies
#' pcd_stan_function_deps("expgrowth_pdf")
pcd_stan_function_deps <- function(
    function_name,
    stan_path = primarycensored::pcd_stan_path()) {
  dep_graph <- .pcd_build_dep_graph(stan_path)

  if (!function_name %in% names(dep_graph)) {
    stop(
      "Function not found: ", function_name,
      call. = FALSE
    )
  }

  .pcd_resolve_deps(function_name, dep_graph)
}

#' Count the number of unmatched braces in a line
#' @noRd
.unmatched_braces <- function(line) {
  ifelse(
    grepl("{", line, fixed = TRUE),
    length(gregexpr("{", line, fixed = TRUE)),
    0
  ) -
    ifelse(
      grepl("}", line, fixed = TRUE),
      length(gregexpr("}", line, fixed = TRUE)),
      0
    )
}

#' Extract function names or content from Stan code
#'
#' @param content Character vector containing Stan code
#'
#' @param names_only Logical, if TRUE extract function names, otherwise
#'  extract function content.
#'
#' @param functions Optional, character vector of function names to extract
#'   content for.
#'
#' @return Character vector of function names or content
#'
#' @keywords internal
.extract_stan_functions <- function(
    content,
    names_only = FALSE,
    functions = NULL) {
  def_pattern <- "^(real|vector|matrix|void|int|array\\s*<\\s*(real|vector|matrix|int)\\s*>|tuple\\s*<\\s*.*\\s*>)\\s+" # nolint
  func_pattern <- paste0(
    def_pattern,
    "(\\w+)\\s*\\("
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
      # Find first opening brace
      repeat {
        brace_count <- brace_count + .unmatched_braces(content[end_line])
        end_line <- end_line + 1
        if (brace_count > 0) break
      }

      # Continue until all braces are closed
      repeat {
        brace_count <- brace_count + .unmatched_braces(content[end_line])
        if (brace_count == 0) break
        end_line <- end_line + 1
      }

      func_content <- c(
        func_content,
        paste(content[start_line:end_line], collapse = "\n")
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
#'  containing Stan files. Defaults to the Stan path of the primarycensored
#'  package.
#'
#' @return A character vector containing unique names of all functions found in
#'  the Stan files.
#'
#' @export
#'
#' @family stantools
pcd_stan_functions <- function(
    stan_path = primarycensored::pcd_stan_path()) {
  stan_files <- list.files(
    stan_path,
    pattern = "\\.stan$",
    full.names = TRUE,
    recursive = TRUE
  )
  functions <- character(0)
  for (file in stan_files) {
    content <- readLines(file)
    functions <- c(
      functions,
      .extract_stan_functions(content, names_only = TRUE)
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
      paste0(stan_path, "/"),
      "",
      matching_files
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
#' @param dependencies Logical, whether to include all functions that the
#' requested functions depend on. When TRUE, recursively finds and includes
#' all dependencies in the correct order (dependencies before the functions
#' that use them). Default is FALSE.
#'
#' @param cpp Logical, whether to declare the functions that have a C++
#' version without a body, so that the C++ versions are used. These are
#' faster. The model must then be compiled with the header from
#' [pcd_stan_header()]. With `dependencies = TRUE`, the Stan functions that
#' only these functions call are left out. Default is FALSE, which gives
#' Stan code that compiles on its own.
#'
#' @return A character string containing the requested Stan functions
#'
#' @family stantools
#'
#' @export
pcd_load_stan_functions <- function(
    functions = NULL,
    stan_path = primarycensored::pcd_stan_path(),
    wrap_in_block = FALSE,
    write_to_file = FALSE,
    output_file = "pcd_functions.stan",
    dependencies = FALSE,
    cpp = FALSE) {
  stan_files <- list.files(
    stan_path,
    pattern = "\\.stan$",
    full.names = TRUE,
    recursive = TRUE
  )

  # Resolve dependencies if requested
  if (dependencies && !is.null(functions)) {
    dep_graph <- .pcd_build_dep_graph(stan_path)
    if (cpp) {
      for (func in intersect(.pcd_cpp_functions(), names(dep_graph))) {
        dep_graph[[func]] <- character(0)
      }
    }

    # Validate that all requested functions exist
    available_funcs <- names(dep_graph)
    not_found <- setdiff(functions, available_funcs)
    if (length(not_found) > 0) {
      stop(
        "Function(s) not found: ", toString(not_found),
        call. = FALSE
      )
    }

    # Resolve all dependencies for each requested function
    all_funcs <- character(0)
    for (func in functions) {
      all_funcs <- .pcd_resolve_deps(func, dep_graph, all_funcs)
    }
    functions <- all_funcs
  }

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
  if (cpp) {
    result <- .pcd_declare_cpp_functions(result)
  }

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
pcd_stan_dist_id <- function(name, type = c("delay", "primary")) {
  type <- match.arg(type)
  lookup <- switch(type,
    delay = primarycensored::pcd_distributions,
    primary = primarycensored::pcd_primary_distributions
  )

  match_idx <- which(lookup$name == name | lookup$aliases == name)

  if (length(match_idx) == 0) {
    stop(
      "No ",
      type,
      " distribution found matching: ",
      name,
      "\n",
      .suggest_dist_name(name, type),
      call. = FALSE
    )
  }

  lookup$stan_id[match_idx]
}
