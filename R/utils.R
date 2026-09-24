#' Extract base function name
#'
#' This helper function extracts the base name of a function, removing
#' namespace prefixes.
#' Base R distribution functions are identified from the C routine they call.
#' Functions exported by another package (for example
#' `flexsurv::pgengamma.orig()`) are identified by their exported name.
#'
#' @inheritParams add_name_attribute
#'
#' @return Character string representing the base name of the function, or
#'  `"unknown"` if it cannot be determined.
#'
#' @keywords internal
.extract_function_name <- function(func) {
  bd <- grep(".Call", deparse(body(func)), value = TRUE, fixed = TRUE)
  if (length(bd) == 1) {
    return(sub("^.*\\.Call\\(C_(\\w+),.+$", "\\1", x = bd))
  }
  env <- environment(func)
  if (!is.null(env) && isNamespace(env)) {
    exports <- getNamespaceExports(env)
    matched <- vapply(
      mget(exports, envir = env, inherits = TRUE, ifnotfound = list(NULL)),
      identical, logical(1),
      y = func
    )
    if (sum(matched) == 1) {
      return(exports[matched])
    }
  }
  return("unknown")
}

#' Get the distribution name of a function
#'
#' Returns the `"name"` attribute of `func` if set. A `stats` function that
#' is identical to one named in [pcd_distributions] or
#' [pcd_primary_distributions] gets that name. Otherwise the name is found
#' with [.extract_function_name()], which deparses the function body and is
#' slower.
#'
#' @inheritParams add_name_attribute
#'
#' @return Character string with the name of the function, or `"unknown"`.
#'
#' @keywords internal
.dist_name <- function(func) {
  name <- attr(func, "name")
  if (!is.null(name)) {
    return(name)
  }
  name <- .registry_name(func)
  if (!is.null(name)) {
    return(name)
  }
  .extract_function_name(func)
}

#' Name a stats function found in the distribution registries
#'
#' Takes the C routine called at the end of the body of a `stats` function
#' (for example `C_pgamma` for [stats::pgamma()]) as the candidate name. The
#' name is returned if it is in [pcd_distributions] or
#' [pcd_primary_distributions] and `func` is identical to the `stats`
#' function of that name.
#'
#' @inheritParams add_name_attribute
#'
#' @return The registry name of `func`, or `NULL` if it is not found.
#'
#' @keywords internal
.registry_name <- function(func) {
  stats_ns <- asNamespace("stats")
  if (!identical(environment(func), stats_ns)) {
    return(NULL)
  }
  expr <- body(func)
  if (is.call(expr) && identical(expr[[1L]], as.name("{"))) {
    expr <- expr[[length(expr)]]
  }
  if (!is.call(expr) || !identical(expr[[1L]], as.name(".Call")) ||
    !is.name(expr[[2L]])) {
    return(NULL)
  }
  # Avoids regular expressions, which are slow relative to the rest
  name <- substring(as.character(expr[[2L]]), 3L)
  # Registered delays are listed by CDF, so match densities by their CDF
  cdf_name <- name
  if (startsWith(name, "d")) {
    cdf_name <- paste0("p", substring(name, 2L))
  }
  primaries <- primarycensored::pcd_primary_distributions
  known <- cdf_name %in% primarycensored::pcd_distributions$pdist ||
    name %in% primaries$dprimary || name %in% primaries$pprimary
  if (!known ||
    !identical(get0(name, envir = stats_ns, inherits = FALSE), func)) {
    return(NULL)
  }
  name
}

#' Helper method for custom distributions
#'
#' [pprimarycensored()] and related functions can identify which distributions
#' are provided via the `pdist` and `dprimary` arguments when those are base R
#' functions (e.g. `punif`, `dexp`) via the `name` attribute.
#'
#' If you need to use a non-base R implementation, but know the distribution
#' name, you can use this helper function to set it in a way that will be
#' detected by [pprimarycensored()] and related functions.
#'
#' This is useful as it enables the automatic use of analytical solutions for
#' distributions where they exist. You can check which analytical solutions are
#' available using `methods(pcens_cdf)` and check distribution names using
#' [pcd_dist_name()].
#'
#' @param func Function, for example the `p`- or `d`- form of a distribution
#' function.
#'
#' @param name Character string, starting with "p" or "d" indicating the
#' underlying distribution.
#'
#' @return Function, with a "name" attribute added
#' @family utils
#' @export
#' @examples
#' dist <- add_name_attribute(pnorm, "hello")
#' attr(dist, "name")
add_name_attribute <- function(func, name) {
  attr(func, "name") <- name
  func
}

#' Extract and Combine Distribution Names
#'
#' This helper function attempts to determine distribution names and uses those
#' to establish a class hierarchy for potential analytical solutions.
#'
#' @inheritParams pprimarycensored
#'
#' @param pdist_name Name of `pdist`, as given by [.dist_name()].
#'
#' @param dprim_name Name of `dprimary`, as given by [.dist_name()].
#'
#' @return A character vector of class names: specific (delay + primary),
#'   delay-only, and base class.
#'
#' @keywords internal
.format_class <- function(pdist, dprimary,
                          pdist_name = .dist_name(pdist),
                          dprim_name = .dist_name(dprimary)) {
  c(
    sprintf("pcens_%s_%s", pdist_name, dprim_name),
    sprintf("pcens_%s", pdist_name),
    "pcens"
  )
}

#' Look up the primary event CDF from the registry
#'
#' Given a primary event density function \code{dprimary}, looks up the
#' corresponding CDF function from \code{pcd_primary_distributions} using
#' the \code{"name"} attribute. Returns \code{NULL} silently when no match
#' is found so that callers can fall back to numerical integration.
#'
#' @param dprimary Function. The primary event density function.
#'
#' @param dprim_name Name of `dprimary`, as given by [.dist_name()].
#'
#' @return A function (the primary CDF) or \code{NULL}.
#'
#' @keywords internal
.lookup_pprimary <- function(dprimary, dprim_name = .dist_name(dprimary)) {
  if (is.null(dprim_name) || dprim_name == "unknown") {
    return(NULL)
  }
  registry <- primarycensored::pcd_primary_distributions
  idx <- which(
    registry$name == dprim_name |
      registry$aliases == dprim_name |
      registry$dprimary == dprim_name
  )
  if (length(idx) == 0L) {
    return(NULL)
  }
  pprimary_name <- registry$pprimary[idx[[1L]]]
  # nocov start: every shipped primary has a `pprimary` set
  if (is.na(pprimary_name)) {
    return(NULL)
  }
  # nocov end
  get0(pprimary_name, envir = asNamespace("primarycensored"))
}

#' Resolve a delay distribution function from a name or function
#'
#' Accepts either a function (returned as-is, with its existing
#' \code{"name"} attribute preserved) or a character string that is looked
#' up against \code{\link{pcd_distributions}}. When a string is supplied,
#' the corresponding base R \code{p<name>} function is returned with the
#' \code{"name"} attribute attached so analytical solutions can dispatch.
#'
#' @param pdist Either a function or a character string.
#' @param type Character string. \code{"p"} for CDF lookup, \code{"d"} for
#'   density. Defaults to \code{"p"}.
#'
#' @return A function with a \code{"name"} attribute.
#'
#' @keywords internal
.resolve_pdist <- function(pdist, type = c("p", "d")) {
  if (is.function(pdist)) {
    return(pdist)
  }
  type <- match.arg(type)
  if (!is.character(pdist) || length(pdist) != 1L) {
    stop(
      "pdist must be a function or a single character string.",
      call. = FALSE
    )
  }
  # Resolve names against `stats` first (covers all base R `p<name>` /
  # `d<name>` lookups in the registry), then `primarycensored` for package
  # extras (`pexpgrowth`, `pdiscretestep`, ...), and finally the global
  # search path. This avoids accidentally picking up an unrelated user
  # object that happens to share the function name.
  .lookup_dist_fn <- function(fn_name) {
    for (env in list(asNamespace("stats"), asNamespace("primarycensored"))) {
      if (exists(fn_name, envir = env, inherits = FALSE, mode = "function")) {
        return(get(fn_name, envir = env, inherits = FALSE, mode = "function"))
      }
    }
    tryCatch(get(fn_name, mode = "function"), error = function(e) NULL)
  }
  registry <- primarycensored::pcd_distributions
  idx <- which(registry$name == pdist | registry$aliases == pdist)
  if (length(idx) == 0L) {
    fn_name <- paste0(type, pdist)
    fn <- .lookup_dist_fn(fn_name)
    if (is.null(fn)) {
      stop(
        "No distribution found matching '", pdist, "'.",
        call. = FALSE
      )
    }
    return(add_name_attribute(fn, fn_name))
  }
  base <- registry$pdist[idx[[1L]]]
  if (is.na(base)) {
    stop(
      "Distribution '", pdist, "' has no base R implementation; ",
      "supply a function instead.",
      call. = FALSE
    )
  }
  fn_name <- if (type == "p") base else sub("^p", "d", base)
  fn <- .lookup_dist_fn(fn_name)
  # nocov start: defensive, the registry's `pdist` column is checked to
  # resolve at package build time so this branch is unreachable in
  # normal operation.
  if (is.null(fn)) {
    stop(
      "Could not find function '", fn_name, "' for distribution '",
      pdist, "'.",
      call. = FALSE
    )
  }
  # nocov end
  add_name_attribute(fn, fn_name)
}

#' Resolve \code{primary_args} / \code{dprimary_args} with a deprecation
#'
#' The resolver distinguishes "user supplied nothing" from "user supplied
#' an empty list" so the deprecated \code{dprimary_args} path can be
#' detected. The returned list is never \code{NULL}.
#'
#' The deprecation is soft (\code{lifecycle::deprecate_soft()}): a warning
#' is shown when the exported function is called from the global
#' environment or from the package under test, and calls from other
#' packages stay silent.
#'
#' @param primary_args The new argument value (or \code{NULL}).
#' @param dprimary_args The old argument value (or \code{NULL}).
#' @param fn Character string identifying the calling function (used in
#'   the deprecation message).
#' @param env Environment of the exported function that owns the
#'   deprecated argument. Defaults to the caller of this helper.
#' @param user_env Environment the exported function was called from.
#'   Defaults to the caller of the caller of this helper.
#'
#' @return A list (possibly empty) of primary distribution arguments.
#'
#' @keywords internal
.resolve_primary_args <- function(primary_args, dprimary_args, fn,
                                  env = parent.frame(),
                                  user_env = parent.frame(2)) {
  has_new <- !is.null(primary_args)
  has_old <- !is.null(dprimary_args)
  if (has_new && has_old) {
    stop(
      "Supply only one of `primary_args` or `dprimary_args`; ",
      "`dprimary_args` is deprecated.",
      call. = FALSE
    )
  }
  if (has_old) {
    lifecycle::deprecate_soft(
      when = "1.6.0",
      what = paste0(fn, "(dprimary_args)"),
      with = paste0(fn, "(primary_args)"),
      env = env,
      user_env = user_env
    )
    return(dprimary_args)
  }
  if (has_new) {
    return(primary_args)
  }
  list()
}

#' Resolve the primary CDF, validating against \code{dprimary} if both supplied
#'
#' Returns the primary CDF to use. If the user supplies \code{pprimary}
#' explicitly (either a function or a string name), it is returned (after
#' resolving the string via \code{\link{pcd_dist_name}}). When both
#' \code{dprimary} and \code{pprimary} carry a \code{"name"} attribute,
#' the names must agree on everything except the leading \code{d}/\code{p};
#' otherwise we error to catch typos like \code{dunif} + \code{pexpgrowth}.
#' If \code{pprimary} is not supplied, falls back to a registry lookup
#' against \code{dprimary} via \code{\link{.lookup_pprimary}}, which may
#' return \code{NULL}.
#'
#' @param dprimary The primary density function.
#' @param pprimary Optional user-supplied primary CDF (function or string).
#' @param d_name Name of `dprimary`, as given by [.dist_name()].
#'
#' @return A primary CDF function, or \code{NULL} if no match was found.
#'
#' @keywords internal
.resolve_pprimary <- function(dprimary, pprimary = NULL,
                              d_name = .dist_name(dprimary)) {
  if (is.null(pprimary)) {
    return(.lookup_pprimary(dprimary, d_name))
  }
  if (is.character(pprimary)) {
    if (length(pprimary) != 1L) {
      stop(
        "pprimary must be a function or a single character string.",
        call. = FALSE
      )
    }
    fn_name <- pcd_dist_name(pprimary, type = "primary")
    fn_name <- sub("^d", "p", fn_name)
    fn <- tryCatch(
      get(fn_name, envir = asNamespace("primarycensored")),
      error = function(e) tryCatch(get(fn_name), error = function(e2) NULL)
    )
    # nocov start: `pcd_dist_name` already errors on unknown names, so a
    # successful lookup followed by a missing p-counterpart is unreachable
    # for any registry entry shipped with the package.
    if (is.null(fn)) {
      stop(
        "Could not find primary CDF function '", fn_name, "'.",
        call. = FALSE
      )
    }
    # nocov end
    pprimary <- add_name_attribute(fn, fn_name)
  }
  if (!is.function(pprimary)) {
    stop(
      "pprimary must be a function or a single character string.",
      call. = FALSE
    )
  }
  .check_primary_names(d_name, pprimary)
  pprimary
}

#' Check that the primary density and CDF refer to the same distribution
#'
#' @param d_name Name of the primary density function, as given by
#'   [.dist_name()].
#' @param pprimary The primary CDF function.
#'
#' @return \code{NULL} invisibly. Called for its error.
#'
#' @keywords internal
.check_primary_names <- function(d_name, pprimary) {
  p_name <- .dist_name(pprimary)
  if (!is.null(d_name) && !is.null(p_name) &&
    d_name != "unknown" && p_name != "unknown" &&
    .strip_prefix(d_name, "d") != .strip_prefix(p_name, "p") &&
    !.same_primary(d_name, p_name)) {
    stop(
      "dprimary and pprimary refer to different distributions: '",
      d_name, "' vs '", p_name, "'.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Check whether two names refer to the same registry primary distribution
#'
#' @param d_name,p_name Names of a primary density and CDF. Each may be a
#'   name, alias, density or CDF name from [pcd_primary_distributions].
#'
#' @return `TRUE` if both names match the same row of
#'   [pcd_primary_distributions], otherwise `FALSE`.
#'
#' @keywords internal
.same_primary <- function(d_name, p_name) {
  registry <- primarycensored::pcd_primary_distributions
  row_of <- function(name) {
    which(
      registry$name == name | registry$aliases == name |
        registry$dprimary == name | registry$pprimary == name
    )[1L]
  }
  d_row <- row_of(d_name)
  !is.na(d_row) && identical(d_row, row_of(p_name))
}

#' Remove a one letter prefix from a distribution name
#'
#' Same as `sub(paste0("^", prefix), "", name)`, without a regular
#' expression for the usual single string.
#'
#' @param name Distribution name.
#' @param prefix Single character prefix, `"d"` or `"p"`.
#'
#' @return `name` without a leading `prefix`.
#'
#' @keywords internal
.strip_prefix <- function(name, prefix) {
  if (!is.character(name) || length(name) != 1L || is.na(name)) {
    return(sub(paste0("^", prefix), "", name))
  }
  if (startsWith(name, prefix)) substring(name, 2L) else name
}

#' Get distribution function cdf or pdf name
#'
#' @param name String. Distribution name or alias
#' @param type String. "delay" or "primary" corresponding to the type of
#'  distribution to use as the look up. If delay then [pcd_distributions()]
#'  is used, if primary then [pcd_primary_distributions()] is used.
#'
#' @return String distribution function name or NA if no base R implementation
#' @export
#' @family utils
#' @examples
#' pcd_dist_name("lnorm")
#' pcd_dist_name("lognormal")
#' pcd_dist_name("gamma")
#' pcd_dist_name("weibull")
#' pcd_dist_name("exp")
#' pcd_dist_name("unif", type = "primary")
#' pcd_dist_name("expgrowth", type = "primary")
pcd_dist_name <- function(name, type = c("delay", "primary")) {
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

  if (type == "delay") {
    lookup$pdist[match_idx]
  } else {
    lookup$dprimary[match_idx]
  }
}

#' @keywords internal
.suggest_dist_name <- function(input, type = "delay") {
  dist_names <- switch(type,
    delay = primarycensored::pcd_distributions$name,
    primary = primarycensored::pcd_primary_distributions$name
  )

  distances <- utils::adist(input, dist_names)
  min_dist <- min(distances)
  candidates <- dist_names[which(distances == min_dist)]

  if (min_dist <= 2 && length(candidates) > 0) {
    suggestions <- paste0(
      "Did you mean: ",
      toString(unique(candidates)),
      "?"
    )
  } else {
    suggestions <- paste0(
      "Available distributions:",
      toString(unique(dist_names))
    )
  }

  return(suggestions)
}
