#' Compute primary event censored CDF
#'
#' This function dispatches to either analytical solutions (if available) or
#' numerical integration via the default method. To see which combinations have
#' analytical solutions implemented, use `methods(pcens_cdf)`. For example,
#' `pcens_cdf.gamma_unif` indicates an analytical solution exists for gamma
#' delay with uniform primary event distributions.
#'
#' @inheritParams pprimarycensored
#'
#' @param object A `primarycensored` object as created by
#' [new_pcens()].
#'
#' @param use_numeric Logical, if TRUE forces use of numeric integration
#' even for distributions with analytical solutions. This is primarily
#' useful for testing purposes or for settings where the analytical solution
#' breaks down.
#'
#' @details
#' When `pwindow = 0` the primary event time is known exactly and the
#' primary event censored CDF is the delay CDF. This case is handled before
#' dispatch, so every method returns `pdist(q)` without integrating over the
#' primary event window.
#'
#' @return Vector of computed primary event censored CDFs
#'
#' @family pcens
#'
#' @export
pcens_cdf <- function(
  object,
  q,
  pwindow,
  use_numeric = FALSE
) {
  if (.is_exact_window(pwindow)) {
    return(.delay_cdf(object, q))
  }
  UseMethod("pcens_cdf")
}

#' Test for a zero-width censoring window
#'
#' @param window A censoring window width.
#'
#' @return `TRUE` if `window` is a single value equal to zero.
#'
#' @keywords internal
.is_exact_window <- function(window) {
  length(window) == 1L && !is.na(window) && window == 0
}

#' Delay CDF of a pcens object
#'
#' Evaluates the delay distribution CDF of a `pcens` object with its stored
#' parameters. This is the primary event censored CDF when `pwindow = 0`.
#'
#' @param object A `pcens` object as created by [new_pcens()].
#'
#' @param q Vector of quantiles.
#'
#' @return Vector of delay CDF values.
#'
#' @keywords internal
.delay_cdf <- function(object, q) {
  do.call(object$pdist, c(list(q), object$args))
}

#' Default method for computing primary event censored CDF
#'
#' This method serves as a fallback for combinations of delay and primary
#' event distributions that don't have specific implementations. It uses
#' a numeric integration method.
#'
#' @inheritParams pcens_cdf
#' @inheritParams pprimarycensored
#'
#' @details
#' This method implements the numerical integration approach for computing
#' the primary event censored CDF. It uses the same mathematical formulation
#' as described in the details section of [pprimarycensored()], but
#' applies numerical integration instead of analytical solutions.
#'
#' @seealso [pprimarycensored()] for the mathematical details of the
#'  primary event censored CDF computation.
#'
#' @family pcens
#'
#' @inherit pcens_cdf return
#'
#' @export
#' @examples
#' # Create a primarycensored object with gamma delay and uniform primary
#' pcens_obj <- new_pcens(
#'   pdist = pgamma,
#'   dprimary = dunif,
#'   primary_args = list(min = 0, max = 1),
#'   shape = 3,
#'   scale = 2
#' )
#'
#' # Compute CDF for a single value
#' pcens_cdf(pcens_obj, q = 9, pwindow = 1)
#'
#' # Compute CDF for multiple values
#' pcens_cdf(pcens_obj, q = c(4, 6, 8), pwindow = 1)
pcens_cdf.default <- function(
  object,
  q,
  pwindow,
  use_numeric = FALSE
) {
  result <- vapply(
    q,
    function(d) {
      integrand <- function(p) {
        d_adj <- d - p
        do.call(object$pdist, c(list(q = d_adj), object$args)) *
          do.call(
            object$dprimary,
            c(list(x = p, min = 0, max = pwindow), object$dprimary_args)
          )
      }
      return(stats::integrate(integrand, lower = 0, upper = pwindow)$value)
    },
    numeric(1)
  )

  # Ensure the result is in [0, 1] (accounts for numerical errors)
  result <- pmin(1, pmax(0, result))

  return(result)
}

#' Method for step CDF delay with general primary event distribution
#'
#' Computes the analytic primary event censored CDF for a piecewise-constant
#' (step) delay distribution and an arbitrary primary event distribution
#' whose CDF \eqn{F_{primary}} is available via \code{object$pprimary}.
#'
#' The observation CDF is
#' \deqn{F_{obs}(q) = \int_0^{pwindow} F_{step}(q-p)\,dF_{primary}(p)}
#' Because \eqn{F_{step}} is piecewise constant, the integral reduces to
#' \deqn{F_{obs}(q) = \sum_k c_k \,[F_{primary}(p^{end}_k) -
#'   F_{primary}(p^{start}_k)]}
#' where \eqn{c_k} is the constant value of \eqn{F_{step}} on the
#' \eqn{k}-th sub-interval of the primary event window induced by the
#' step-function knots.
#' The partition is exact for any bin widths, so bins may be wider or
#' narrower than \code{pwindow}, and for boundaries that start below zero.
#'
#' Falls back to \code{pcens_cdf.default} when \code{use_numeric = TRUE}
#' or when no primary CDF is available on the object.
#'
#' @inheritParams pcens_cdf
#'
#' @family pcens
#'
#' @inherit pcens_cdf return
#'
#' @export
pcens_cdf.pcens_pdiscretestep <- function(
    object,
    q,
    pwindow,
    use_numeric = FALSE) {
  if (isTRUE(use_numeric)) {
    return(pcens_cdf.default(object, q, pwindow, use_numeric))
  }

  pprimary <- object$pprimary
  if (is.null(pprimary)) {
    return(pcens_cdf.default(object, q, pwindow, use_numeric))
  }

  boundaries <- object$args$boundaries
  pmf <- object$args$pmf

  if (is.null(boundaries) || is.null(pmf)) {
    stop(
      "boundaries and pmf are required for the step distribution.",
      call. = FALSE
    )
  }

  K <- length(pmf)
  cum_pmf <- cumsum(pmf)
  right_edges <- boundaries[-1L]
  primary_args <- object$primary_args

  # Evaluate F_primary(p) on [0, pwindow] via the stored pprimary function.
  # primary_args may contain r, shape, etc.; translate min/max to the
  # primary window.
  .F_primary <- function(p) {
    do.call(
      pprimary,
      c(list(q = p, min = 0, max = pwindow), primary_args)
    )
  }

  # Helper: evaluate F_step at a single point (right-continuous).
  .fstep <- function(x) {
    idx <- findInterval(x, right_edges, left.open = FALSE)
    if (idx == 0L) 0 else if (idx >= K) 1 else cum_pmf[idx]
  }

  result <- vapply(
    q,
    function(qi) {
      # Breakpoints in primary-time p where F_step(qi - p) changes value:
      # these are p = qi - right_edge[k], restricted to (0, pwindow).
      # Knots outside the window do not matter because F_step(qi - p) is
      # then constant on the whole sub-interval, whatever the bin width.
      p_knots <- qi - right_edges
      inside <- p_knots[p_knots > 0 & p_knots < pwindow]
      breaks <- sort(unique(c(0, inside, pwindow)))
      total <- 0
      for (j in seq_len(length(breaks) - 1L)) {
        a <- breaks[j]
        b <- breaks[j + 1L]
        mid <- 0.5 * (a + b)
        c_k <- .fstep(qi - mid)
        total <- total + c_k * (.F_primary(b) - .F_primary(a))
      }
      total
    },
    numeric(1)
  )

  pmin(1, pmax(0, result))
}

#' Method for hazard CDF delay with general primary event distribution
#'
#' Computes the analytic primary event censored CDF for a hazard-parameterised
#' piecewise-constant delay distribution. Converts hazards to a PMF via
#' [hazards_to_pmf()] then dispatches back through [pcens_cdf()] using a
#' freshly constructed step-distribution object. The same primary event
#' distribution and arguments are preserved.
#'
#' @inheritParams pcens_cdf
#'
#' @family pcens
#'
#' @inherit pcens_cdf return
#'
#' @export
pcens_cdf.pcens_pdiscretehazard <- function(
    object,
    q,
    pwindow,
    use_numeric = FALSE) {
  if (isTRUE(use_numeric)) {
    return(pcens_cdf.default(object, q, pwindow, use_numeric))
  }
  hazards <- object$args$hazards
  if (is.null(hazards)) {
    stop(
      "hazards are required for the hazard distribution.",
      call. = FALSE
    )
  }
  # `hazards_to_pmf()` is called exactly once per dispatch (cached into
  # `new_args$pmf`); the recursive `pcens_cdf()` call below then sweeps the
  # full q vector through the step method without re-running the
  # hazards->pmf conversion per element.
  new_args <- object$args
  new_args$hazards <- NULL
  new_args$pmf <- hazards_to_pmf(hazards)
  step_obj <- do.call(
    new_pcens,
    c(
      list(
        pdist = pdiscretestep,
        dprimary = object$dprimary,
        primary_args = object$dprimary_args,
        pprimary = object$pprimary
      ),
      new_args
    )
  )
  pcens_cdf(step_obj, q, pwindow, use_numeric)
}

#' Method for Gamma delay with uniform primary
#'
#' @inheritParams pcens_cdf
#'
#' @family pcens
#'
#' @inherit pcens_cdf return
#'
#' @export
pcens_cdf.pcens_pgamma_dunif <- function(
  object,
  q,
  pwindow,
  use_numeric = FALSE
) {
  if (isTRUE(use_numeric)) {
    return(pcens_cdf.default(object, q, pwindow, use_numeric))
  }
  .pcens_cdf_shared(.uniform_terms_gamma(object$args), q, pwindow)
}

#' Method for Log-Normal delay with uniform primary
#'
#' @inheritParams pcens_cdf
#'
#' @family pcens
#'
#' @inherit pcens_cdf return
#'
#' @export
pcens_cdf.pcens_plnorm_dunif <- function(
  object,
  q,
  pwindow,
  use_numeric = FALSE
) {
  if (isTRUE(use_numeric)) {
    return(pcens_cdf.default(object, q, pwindow, use_numeric))
  }
  .pcens_cdf_shared(.uniform_terms_lnorm(object$args), q, pwindow)
}

#' Method for Weibull delay with uniform primary
#'
#' @inheritParams pcens_cdf
#'
#' @family pcens
#'
#' @inherit pcens_cdf return
#'
#' @importFrom stats pgamma
#'
#' @export
pcens_cdf.pcens_pweibull_dunif <- function(
  object,
  q,
  pwindow,
  use_numeric = FALSE
) {
  if (isTRUE(use_numeric)) {
    return(pcens_cdf.default(object, q, pwindow, use_numeric))
  }
  .pcens_cdf_shared(.uniform_terms_weibull(object$args), q, pwindow)
}

#' Method for generalised gamma delay with uniform primary
#'
#' Analytical solution for the generalised gamma distribution in the Stacy
#' parameterisation used by `flexsurv::pgengamma.orig()`, with parameters
#' `shape`, `scale` and `k`.
#' The delay CDF is \eqn{F_T(t) = P(k, (t / \theta)^a)} with \eqn{P} the
#' regularised lower incomplete gamma function, \eqn{a} the `shape` and
#' \eqn{\theta} the `scale`.
#' The mean is \eqn{E[T] = \theta \Gamma(k + 1/a) / \Gamma(k)} and the partial
#' expectation distribution is \eqn{\tilde F_T(t) = P(k + 1/a, (t / \theta)^a)},
#' so the solution generalises the gamma (`shape = 1`) and Weibull (`k = 1`)
#' cases.
#' See `vignette("analytic-solutions")` for the derivation.
#'
#' @inheritParams pcens_cdf
#'
#' @family pcens
#'
#' @inherit pcens_cdf return
#'
#' @export
#' @examplesIf requireNamespace("flexsurv", quietly = TRUE)
#' pcens_obj <- new_pcens(
#'   pdist = flexsurv::pgengamma.orig,
#'   dprimary = dunif,
#'   dprimary_args = list(min = 0, max = 1),
#'   shape = 1.5,
#'   scale = 2,
#'   k = 0.8
#' )
#' pcens_cdf(pcens_obj, q = c(1, 4, 8), pwindow = 1)
pcens_cdf.pcens_pgengamma.orig_dunif <- function(
  object,
  q,
  pwindow,
  use_numeric = FALSE
) {
  if (isTRUE(use_numeric)) {
    return(pcens_cdf.default(object, q, pwindow, use_numeric))
  }
  .pcens_cdf_shared(
    .uniform_terms_gengamma(
      object$args$shape, object$args$scale, object$args$k
    ),
    q, pwindow
  )
}

#' Method for generalised gamma (Prentice parameterisation) delay with
#' uniform primary
#'
#' Analytical solution for the generalised gamma distribution in the Prentice
#' parameterisation used by `flexsurv::pgengamma()`, with parameters `mu`,
#' `sigma` and `Q`.
#' For `Q > 0` this is mapped to the Stacy parameterisation of
#' [pcens_cdf.pcens_pgengamma.orig_dunif()] via `shape = Q / sigma`,
#' `scale = exp(mu) * Q^(2 * sigma / Q)` and `k = 1 / Q^2`.
#' For `Q <= 0` (the lognormal and reflected cases) the numerical
#' [pcens_cdf.default()] method is used.
#'
#' @inheritParams pcens_cdf
#'
#' @family pcens
#'
#' @inherit pcens_cdf return
#'
#' @export
pcens_cdf.pcens_pgengamma_dunif <- function(
  object,
  q,
  pwindow,
  use_numeric = FALSE
) {
  spec <- if (!isTRUE(use_numeric)) .uniform_terms_prentice(object$args)
  if (is.null(spec)) {
    return(pcens_cdf.default(object, q, pwindow, use_numeric))
  }
  .pcens_cdf_shared(spec, q, pwindow)
}
