#' S3 class for primary event censored distribution computation
#'
#' @inheritParams pprimarycensored
#'
#' @param pprimary CDF of the primary event distribution. May be a function
#'   or a character string naming a primary distribution in
#'   \code{pcd_primary_distributions}. When \code{NULL} (the default), it is
#'   looked up automatically from the registry using the \code{"name"}
#'   attribute of \code{dprimary}. When both \code{dprimary} and
#'   \code{pprimary} carry a name, the two must agree on everything other
#'   than the leading \code{d}/\code{p} prefix; mismatches such as
#'   \code{dunif} + \code{pexpgrowth} raise an error.
#'
#' @param primary_args List of additional arguments to be passed to
#'   \code{dprimary} (and the looked-up \code{pprimary}). Replaces the
#'   deprecated \code{dprimary_args}.
#'
#' @param dprimary_args \[Deprecated\] Use \code{primary_args} instead.
#'
#' @return An object with class hierarchy
#'   \code{c("pcens_{pdist_name}_{dprimary_name}", "pcens_{pdist_name}",
#'   "pcens")}. It is a list with the fields:
#'   \describe{
#'     \item{\code{pdist}}{The delay distribution CDF.}
#'     \item{\code{dprimary}}{The primary event distribution density.}
#'     \item{\code{primary_args}}{A list of arguments passed to
#'       \code{dprimary} and \code{pprimary}.}
#'     \item{\code{dprimary_args}}{A copy of \code{primary_args}, kept for
#'       backward compatibility.}
#'     \item{\code{pprimary}}{The primary event CDF, or \code{NULL} if none
#'       is available.}
#'     \item{\code{args}}{A named list of the delay distribution parameters
#'       passed through \code{...}.}
#'   }
#'   It can be used with [pcens_cdf()] to compute the primary event
#'   censored CDF and with [pcens_pmf()] to compute the PMF. Use
#'   [update()][update.pcens()] to change its parameters.
#'
#' @family pcens
#'
#' @export
#' @examples
#' new_pcens(
#'   pdist = pgamma, dprimary = dunif,
#'   primary_args = list(min = 0, max = 1),
#'   shape = 1, scale = 1
#' )
new_pcens <- function(
    pdist,
    dprimary,
    primary_args = NULL,
    pprimary = NULL,
    dprimary_args = NULL,
    ...) {
  primary_args <- .resolve_primary_args(
    primary_args, dprimary_args, "new_pcens"
  )
  pprimary <- .resolve_pprimary(
    dprimary, pprimary
  )
  obj <- list(
    pdist = pdist,
    dprimary = dprimary,
    primary_args = primary_args,
    # Keep dprimary_args alias for backward compatibility with downstream
    # consumers that read object$dprimary_args directly.
    dprimary_args = primary_args,
    pprimary = pprimary,
    args = list(...)
  )
  class(obj) <- .format_class(pdist, dprimary)
  obj
}

#' Update the parameters of a pcens object
#'
#' Replaces the delay distribution parameters, and optionally the primary
#' event distribution arguments, of an existing `pcens` object.
#' The delay and primary distribution functions, the primary event CDF and
#' the class are kept as they are.
#' The updated object therefore dispatches to the same [pcens_cdf()] method
#' without looking up names or rebuilding the class.
#' This is cheaper than calling [new_pcens()] again when one distribution is
#' evaluated for many parameter sets, for example posterior draws.
#'
#' @param object A `pcens` object as created by [new_pcens()].
#'
#' @param ... Named delay distribution parameters. Each one replaces the
#'   entry of the same name in `object$args`, or is added if not present.
#'   Parameters that are not given keep their current values.
#'
#' @param primary_args Optional named list of primary event distribution
#'   arguments. These are merged into `object$primary_args` in the same way
#'   as `...` is merged into `object$args`. Defaults to `NULL`, which leaves
#'   the primary event distribution arguments unchanged.
#'
#' @details
#' Parameters are merged rather than replaced as a whole, so
#' `update(object, scale = 3)` changes `scale` and keeps all other
#' parameters. Parameters cannot be removed; use [new_pcens()] for that.
#'
#' A name in `...` that is not already in `object$args` must be an argument
#' of `object$pdist`, unless `pdist` takes `...`. Otherwise an error is
#' raised. Names in `primary_args` are not checked against `dprimary`.
#'
#' @return A `pcens` object with the same class as `object` and updated
#'   `args`, `primary_args` and `dprimary_args` fields. See [new_pcens()] for
#'   the fields of a `pcens` object.
#'
#' @family pcens
#'
#' @importFrom stats update
#' @export
#' @examples
#' obj <- new_pcens(
#'   pdist = pgamma, dprimary = dunif,
#'   primary_args = list(min = 0, max = 1),
#'   shape = 1, scale = 1
#' )
#' obj <- update(obj, shape = 2, scale = 3)
#' pcens_cdf(obj, q = c(1, 5, 10), pwindow = 1)
#'
#' # Update the primary event distribution arguments
#' obj <- new_pcens(
#'   pdist = pgamma, dprimary = dexpgrowth,
#'   primary_args = list(r = 0.2),
#'   shape = 2, scale = 3
#' )
#' obj <- update(obj, primary_args = list(r = 0.5))
#' pcens_pmf(obj, x = 0:5, pwindow = 1)
update.pcens <- function(object, ..., primary_args = NULL) {
  new_args <- list(...)
  if (length(new_args) > 0L) {
    .check_named_list(new_args, "Delay parameters passed to update()")
    unknown <- setdiff(names(new_args), names(object$args))
    if (length(unknown) > 0L) {
      pdist_args <- names(formals(object$pdist))
      if (!is.null(pdist_args) && !"..." %in% pdist_args) {
        unknown <- setdiff(unknown, pdist_args)
        if (length(unknown) > 0L) {
          stop(
            "Unknown delay parameter(s) for pdist: ", toString(unknown), ".",
            call. = FALSE
          )
        }
      }
    }
    object$args[names(new_args)] <- new_args
  }
  if (!is.null(primary_args)) {
    if (!is.list(primary_args)) {
      stop("primary_args must be a list.", call. = FALSE)
    }
    if (length(primary_args) > 0L) {
      .check_named_list(primary_args, "primary_args")
      object$primary_args[names(primary_args)] <- primary_args
      object$dprimary_args <- object$primary_args
    }
  }
  object
}

#' Check that every element of a list is named
#'
#' @param x A list.
#'
#' @param what Character string describing `x` for the error message.
#'
#' @return `NULL` invisibly. Called for its error.
#'
#' @keywords internal
.check_named_list <- function(x, what) {
  nms <- names(x)
  if (is.null(nms) || anyNA(nms) || !all(nzchar(nms))) {
    stop(what, " must be named.", call. = FALSE)
  }
  invisible(NULL)
}
