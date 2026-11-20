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
    # Match the whole name, whether called or passed to another function
    # (as with ODE right-hand sides)
    pattern <- paste0("\\b", fname, "\\b")
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
    func_names <- unique(.extract_stan_functions(content, names_only = TRUE))

    for (fname in func_names) {
      func_content <- .extract_stan_functions(
        content,
        names_only = FALSE,
        functions = fname
      )
      if (length(func_content) > 0) {
        # Find dependencies, excluding self-references
        deps <- .pcd_stan_function_deps(
          paste(func_content, collapse = "\n"), all_func_names
        )
        deps <- setdiff(deps, fname)
        dep_graph[[fname]] <- union(dep_graph[[fname]], deps)
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
  count <- function(brace) {
    lengths(regmatches(line, gregexpr(brace, line, fixed = TRUE)))
  }
  count("{") - count("}")
}

#' Locate the Stan function definitions in Stan code
#'
#' @param content Character vector containing Stan code
#'
#' @return A data frame with one row per definition, in order, giving the
#'   function `name`, the `start` and `end` lines of the definition and its
#'   `signature` (the header up to the opening brace).
#'
#' @keywords internal
.stan_function_defs <- function(content) {
  def_pattern <- "^(real|vector|matrix|void|int|array\\s*<\\s*(real|vector|matrix|int)\\s*>|tuple\\s*<\\s*.*\\s*>)\\s+" # nolint
  first <- grep(paste0(def_pattern, "\\w+\\s*\\("), content)
  name <- sub("\\s*\\(.*$", "", sub(def_pattern, "", content[first]))
  braces <- cumsum(vapply(content, .unmatched_braces, numeric(1),
                           USE.NAMES = FALSE))
  last <- integer(length(first))
  header <- character(length(first))
  for (i in seq_along(first)) {
    before <- if (first[i] > 1) braces[first[i] - 1] else 0
    rest <- first[i]:length(content)
    opened <- rest[grepl("{", content[rest], fixed = TRUE)][1]
    last[i] <- rest[rest >= opened & braces[rest] == before][1]
    header[i] <- paste(content[first[i]:opened], collapse = " ")
  }
  data.frame(
    name = name, start = first, end = last,
    signature = trimws(gsub("\\s+", " ", sub("\\{.*$", "", header))),
    stringsAsFactors = FALSE
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
#' @return Character vector of function names or content, with one element
#'   per definition, so overloaded names appear more than once.
#'
#' @keywords internal
.extract_stan_functions <- function(
    content,
    names_only = FALSE,
    functions = NULL) {
  defs <- .stan_function_defs(content)
  if (!is.null(functions)) {
    defs <- defs[defs$name %in% functions, , drop = FALSE]
  }
  if (names_only) {
    return(defs$name)
  }
  .stan_def_text(content, defs)
}

#' Get the code of Stan function definitions
#'
#' @inheritParams .stan_function_defs
#' @param defs A data frame of definitions from `.stan_function_defs()`.
#'
#' @return Character vector with the code of each definition.
#'
#' @keywords internal
.stan_def_text <- function(content, defs) {
  vapply(
    seq_len(nrow(defs)),
    function(i) paste(content[defs$start[i]:defs$end[i]], collapse = "\n"),
    character(1)
  )
}

#' Choose which definitions of overloaded Stan functions to keep
#'
#' Stan allows several definitions of a function with the same name and
#' different argument types, but cmdstanr cannot expose such functions to R.
#' This keeps one definition of each overloaded name, the first unless
#' `overloads` says otherwise.
#'
#' @inheritParams .stan_def_text
#' @inheritParams pcd_load_stan_functions
#'
#' @return A logical vector, `TRUE` for the definitions to keep.
#'
#' @keywords internal
.select_stan_overloads <- function(defs, overloads = list()) {
  unknown <- setdiff(names(overloads), defs$name)
  if (length(unknown) > 0) {
    stop("Function(s) not found: ", toString(unknown), call. = FALSE)
  }
  keep <- rep(TRUE, NROW(defs))
  for (name in unique(defs$name[duplicated(defs$name)])) {
    idx <- which(defs$name == name)
    choice <- overloads[[name]]
    if (is.null(choice)) {
      choice <- 1L
    }
    chosen <- if (is.character(choice)) {
      idx[grepl(choice, defs$signature[idx], fixed = TRUE)]
    } else {
      idx[choice]
    }
    if (length(chosen) == 0 || anyNA(chosen)) {
      stop(
        "No definition of ", name, " matches the overload choice ",
        toString(choice), ". Its definitions are:\n",
        paste(defs$signature[idx], collapse = "\n"),
        call. = FALSE
      )
    }
    keep[setdiff(idx, chosen)] <- FALSE
  }
  keep
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
#' @param overloads A named list choosing which definitions to keep of
#' functions that have several definitions with different argument types
#' (overloads). Each element is named after a function and is either the
#' position of the definitions to keep, in the order they appear in the Stan
#' files, or a string matched against their signatures. Only one definition
#' of each overloaded function is kept by default, the first, as cmdstanr
#' cannot expose overloaded functions to R. For example,
#' `list(primarycensored_lpmf = 2)` or
#' `list(primarycensored_lpmf = "array[] int n")` keeps the grouped
#' `primarycensored_lpmf()` in place of the per-delay one, and
#' `list(primarycensored_lpmf = 1:2)` keeps both.
#' The code does not compile if another loaded function calls a definition
#' that is left out, so to keep only a later definition load that function
#' with `dependencies = TRUE` rather than loading all functions.
#'
#' @return A character string containing the requested Stan functions
#'
#' @family stantools
#'
#' @export
#'
#' @examples
#' # The grouped overload of primarycensored_lpmf
#' pcd_load_stan_functions(
#'   "primarycensored_lpmf",
#'   overloads = list(primarycensored_lpmf = "array[] int n")
#' )
pcd_load_stan_functions <- function(
    functions = NULL,
    stan_path = primarycensored::pcd_stan_path(),
    wrap_in_block = FALSE,
    write_to_file = FALSE,
    output_file = "pcd_functions.stan",
    dependencies = FALSE,
    overloads = list()) {
  stan_files <- list.files(
    stan_path,
    pattern = "\\.stan$",
    full.names = TRUE,
    recursive = TRUE
  )

  # Resolve dependencies if requested
  if (dependencies && !is.null(functions)) {
    dep_graph <- .pcd_build_dep_graph(stan_path)

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

  contents <- lapply(stan_files, readLines)
  defs <- lapply(contents, .stan_function_defs)
  file_id <- rep(seq_along(defs), vapply(defs, nrow, integer(1)))
  keep <- split(
    .select_stan_overloads(do.call(rbind, defs), overloads),
    factor(file_id, levels = seq_along(defs))
  )

  all_content <- character(0)

  for (i in seq_along(contents)) {
    content <- contents[[i]]
    file_defs <- defs[[i]][keep[[i]], , drop = FALSE]
    if (is.null(functions)) {
      dropped <- defs[[i]][!keep[[i]], , drop = FALSE]
      drop_lines <- unlist(Map(seq, dropped$start, dropped$end))
      all_content <- c(all_content, content[setdiff(seq_along(content),
                                                     drop_lines)])
    } else {
      for (func in functions) {
        all_content <- c(all_content, .stan_def_text(
          content, file_defs[file_defs$name == func, , drop = FALSE]
        ))
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
