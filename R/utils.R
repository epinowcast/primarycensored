#' Extract base function name
#'
#' This helper function extracts the base name of a function, removing
#' namespace prefixes.
#'
#' @inheritParams add_name_attribute
#'
#' @return Character string representing the base name of the function.
#'
#' @keywords internal
.extract_function_name <- function(func) {
  bd <- grep(".Call", deparse(body(func)), value = TRUE, fixed = TRUE)
  if (length(bd) == 1) {
    return(sub("^.*\\.Call\\(C_(\\w+),.+$", "\\1", x = bd))
  } else {
    return("unknown")
  }
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
#' @return A character vector of class names: specific (delay + primary),
#'   delay-only, and base class.
#'
#' @keywords internal
.format_class <- function(pdist, dprimary) {
  pdist_name <- attr(pdist, "name")
  if (is.null(pdist_name)) {
    pdist_name <- .extract_function_name(pdist)
  }
  dprim_name <- attr(dprimary, "name")
  if (is.null(dprim_name)) {
    dprim_name <- .extract_function_name(dprimary)
  }
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
#' @return A function (the primary CDF) or \code{NULL}.
#'
#' @keywords internal
.lookup_pprimary <- function(dprimary) {
  dprim_name <- attr(dprimary, "name")
  if (is.null(dprim_name)) {
    dprim_name <- .extract_function_name(dprimary)
  }
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
  if (is.na(pprimary_name)) {
    return(NULL)
  }
  fn <- tryCatch(get(pprimary_name, envir = asNamespace("primarycensored")),
    error = function(e) NULL
  )
  fn
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
  type <- match.arg(type)
  if (is.function(pdist)) {
    return(pdist)
  }
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
  if (is.null(fn)) {
    stop(
      "Could not find function '", fn_name, "' for distribution '",
      pdist, "'.",
      call. = FALSE
    )
  }
  add_name_attribute(fn, fn_name)
}

#' Resolve \code{primary_args} / \code{dprimary_args} with a deprecation
#'
#' Maps the deprecated \code{dprimary_args} formal to the new
#' \code{primary_args}. Emits a deprecation warning when the old name is
#' used, errors if both are supplied, and returns the resolved list.
#'
#' The new \code{primary_args} formal defaults to \code{NULL} (rather than
#' \code{list()}) so this resolver can distinguish "user supplied nothing"
#' from "user supplied an empty list" alongside the deprecated
#' \code{dprimary_args}. The returned list is never \code{NULL} (an empty
#' list is returned when neither argument is supplied), so downstream code
#' can treat the result as a list unconditionally.
#'
#' @param primary_args The new argument value (or \code{NULL}).
#' @param dprimary_args The old argument value (or \code{NULL}).
#' @param fn Character string identifying the calling function (used in
#'   the deprecation message).
#'
#' @return A list (possibly empty) of primary distribution arguments.
#'
#' @keywords internal
.resolve_primary_args <- function(primary_args, dprimary_args, fn) {
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
    if (requireNamespace("lifecycle", quietly = TRUE)) {
      lifecycle::deprecate_warn(
        when = "1.5.0",
        what = paste0(fn, "(dprimary_args)"),
        with = paste0(fn, "(primary_args)")
      )
    } else {
      warning(
        "Argument `dprimary_args` is deprecated; use `primary_args`.",
        call. = FALSE
      )
    }
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
#'
#' @return A primary CDF function, or \code{NULL} if no match was found.
#'
#' @keywords internal
.resolve_pprimary <- function(dprimary, pprimary = NULL) {
  if (is.null(pprimary)) {
    return(.lookup_pprimary(dprimary))
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
    if (is.null(fn)) {
      stop(
        "Could not find primary CDF function '", fn_name, "'.",
        call. = FALSE
      )
    }
    pprimary <- add_name_attribute(fn, fn_name)
  }
  if (!is.function(pprimary)) {
    stop(
      "pprimary must be a function or a single character string.",
      call. = FALSE
    )
  }
  d_name <- attr(dprimary, "name")
  if (is.null(d_name)) d_name <- .extract_function_name(dprimary)
  p_name <- attr(pprimary, "name")
  if (is.null(p_name)) p_name <- .extract_function_name(pprimary)
  if (!is.null(d_name) && !is.null(p_name) &&
    d_name != "unknown" && p_name != "unknown") {
    if (sub("^d", "", d_name) != sub("^p", "", p_name)) {
      stop(
        "dprimary and pprimary refer to different distributions: '",
        d_name, "' vs '", p_name, "'.",
        call. = FALSE
      )
    }
  }
  pprimary
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
