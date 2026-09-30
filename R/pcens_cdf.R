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

#' Uniform primary event censored CDF from terms at the window ends
#'
#' For a delay with mean \eqn{E} and a uniform primary event window of width
#' \eqn{w}, the CDF at \eqn{d} is \eqn{(G(d) - G(q)) / w} with
#' \eqn{q = \max(d - w, 0)} and \eqn{G(t) = t F(t) - E \tilde F(t)}.
#' \eqn{F} is the delay CDF and \eqn{\tilde F} is the CDF of the partial
#' expectation distribution.
#'
#' Delays at or below zero give 0, delays of `Inf` give 1 and missing delays
#' are an error.
#' Zero-width windows give the delay CDF.
#' `q` and `pwindow` are recycled against each other.
#'
#' @inheritParams pcens_cdf
#'
#' @param terms_fn Function of delays greater than or equal to zero returning
#'  \eqn{G}.
#'
#' @param upper_fn Function of delays greater than or equal to zero returning
#'  \eqn{H(t) = t S(t) - E \tilde S(t)}, where \eqn{S} and \eqn{\tilde S} are
#'  the upper tails of \eqn{F} and \eqn{\tilde F}.
#'
#' @param mean Mean of the delay distribution, \eqn{E}.
#'
#' @param delay_cdf Function of delays returning the delay CDF.
#'
#' @return Numeric vector of CDF values in \[0, 1\].
#'
#' @keywords internal
.pcens_cdf_uniform <- function(q, pwindow, terms_fn, upper_fn, mean,
                               delay_cdf) {
  .check_pwindow(q, pwindow)
  if (length(q) == 0L) {
    return(numeric(0))
  }
  n <- max(length(q), length(pwindow))
  if (n %% length(q) != 0L || n %% length(pwindow) != 0L) {
    warning(
      "longer object length is not a multiple of shorter object length",
      call. = FALSE
    )
  }
  q <- rep_len(q, n)
  pwindow <- rep_len(pwindow, n)
  active <- q > 0 & q < Inf
  exact <- active & pwindow == 0
  if (!any(exact) && all(active)) {
    result <- .uniform_window_cdf(
      q, pwindow, terms_fn, upper_fn, mean, delay_cdf
    )
    return(pmin.int(1, pmax.int(0, result)))
  }

  result <- numeric(n)
  result[q == Inf] <- 1
  if (any(exact)) {
    result[exact] <- delay_cdf(q[exact])
    active <- active & !exact
  }
  if (any(active)) {
    result[active] <- .uniform_window_cdf(
      q[active], pwindow[active], terms_fn, upper_fn, mean, delay_cdf
    )
  }
  pmin.int(1, pmax.int(0, result))
}

# Errors if q has missing values or pwindow is empty, missing or negative
.check_pwindow <- function(q, pwindow) {
  if (anyNA(q)) {
    stop("q must not contain missing values.", call. = FALSE)
  }
  if (length(pwindow) == 0L || anyNA(pwindow) || min(pwindow) < 0) {
    stop(
      "pwindow must be non-negative with no missing values.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

# Delays more than this many windows wide may use the survival form
.upper_tail_ratio <- 1e3
# Delays more than this many windows wide use the mean of the delay CDF
.narrow_window_ratio <- 1e6

# Uniform primary CDF for finite delays and windows greater than 0, with one
# window per delay. Values are not clamped to [0, 1].
.uniform_window_cdf <- function(d, w, terms_fn, upper_fn, mean, delay_cdf) {
  # Both forms lose about 1e-16 * d / w for narrow windows (1e-10 at the
  # cutoff), so the mean of the delay CDF over the window is used instead
  narrow <- d > .narrow_window_ratio * w
  if (any(narrow)) {
    if (all(narrow)) {
      return(.narrow_window_cdf(d, w, delay_cdf))
    }
    result <- numeric(length(d))
    result[narrow] <- .narrow_window_cdf(d[narrow], w[narrow], delay_cdf)
    result[!narrow] <- .uniform_window_cdf(
      d[!narrow], w[!narrow], terms_fn, upper_fn, mean, delay_cdf
    )
    return(result)
  }
  lo <- pmax.int(d - w, 0)
  # In the far upper tail G(d) - G(lo) cancels, unlike 1 - (H(d) - H(lo)) / w
  # with G(t) = t - mean - H(t). This only pays off past the delay mean.
  up <- lo > mean & d > .upper_tail_ratio * w
  if (!any(up)) {
    return((terms_fn(d) - terms_fn(lo)) / w)
  }
  if (all(up)) {
    return(1 - (upper_fn(d) - upper_fn(lo)) / w)
  }
  result <- numeric(length(d))
  result[up] <- 1 - (upper_fn(d[up]) - upper_fn(lo[up])) / w[up]
  result[!up] <- (terms_fn(d[!up]) - terms_fn(lo[!up])) / w[!up]
  result
}

# Mean of the delay CDF over [d - w, d] from a 5 point Gauss-Legendre rule
.narrow_window_cdf <- function(d, w, delay_cdf) {
  # Weights are halved so they sum to 1
  nodes <- c(-0.906179845938664, -0.5384693101056831, 0,
             0.5384693101056831, 0.906179845938664)
  gl_weights <- c(0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
                  0.4786286704993665, 0.2369268850561891) / 2
  n <- length(d)
  half <- w / 2
  x <- rep(d - half, 5) + rep(nodes, each = n) * rep(half, 5)
  drop(matrix(delay_cdf(x), nrow = n) %*% gl_weights)
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
    return(
      pcens_cdf.default(object, q, pwindow, use_numeric)
    )
  }
  # Extract Gamma distribution parameters
  shape <- object$args$shape
  scale <- object$args$scale
  rate <- object$args$rate
  # if we don't have scale get fromm rate
  if (is.null(scale) && !is.null(rate)) {
    scale <- 1 / rate
  }
  if (is.null(shape)) {
    stop("shape parameter is required for Gamma distribution", call. = FALSE)
  }
  if (is.null(scale)) {
    stop(
      "scale or rate parameter is required for Gamma distribution",
      call. = FALSE
    )
  }

  # The partial expectation distribution is a gamma with shape + 1
  E_T <- shape * scale
  terms_fn <- function(t) {
    t * pgamma(t, shape, scale = scale) -
      E_T * pgamma(t, shape + 1, scale = scale)
  }
  upper_fn <- function(t) {
    t * pgamma(t, shape, scale = scale, lower.tail = FALSE) -
      E_T * pgamma(t, shape + 1, scale = scale, lower.tail = FALSE)
  }
  delay_cdf <- function(t) pgamma(t, shape, scale = scale)

  .pcens_cdf_uniform(q, pwindow, terms_fn, upper_fn, E_T, delay_cdf)
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
    return(
      pcens_cdf.default(object, q, pwindow, use_numeric)
    )
  }

  # Extract Log-Normal distribution parameters
  mu <- object$args$meanlog
  sigma <- object$args$sdlog
  if (is.null(mu)) {
    stop(
      "meanlog parameter is required for Log-Normal distribution",
      call. = FALSE
    )
  }
  if (is.null(sigma)) {
    stop(
      "sdlog parameter is required for Log-Normal distribution",
      call. = FALSE
    )
  }

  # The partial expectation distribution has meanlog + sdlog^2, which gives
  # z - sdlog for the standardised log delay z
  E_T <- exp(mu + 0.5 * sigma^2)
  terms_fn <- function(t) {
    z <- (log(t) - mu) / sigma
    t * stats::pnorm(z) - E_T * stats::pnorm(z - sigma)
  }
  upper_fn <- function(t) {
    z <- (log(t) - mu) / sigma
    t * stats::pnorm(z, lower.tail = FALSE) -
      E_T * stats::pnorm(z - sigma, lower.tail = FALSE)
  }
  delay_cdf <- function(t) stats::plnorm(t, mu, sigma)

  .pcens_cdf_uniform(q, pwindow, terms_fn, upper_fn, E_T, delay_cdf)
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
    return(
      pcens_cdf.default(object, q, pwindow, use_numeric)
    )
  }

  # Extract Weibull distribution parameters
  shape <- object$args$shape
  scale <- object$args$scale
  if (is.null(shape)) {
    stop("shape parameter is required for Weibull distribution", call. = FALSE)
  }
  if (is.null(scale)) {
    stop("scale parameter is required for Weibull distribution", call. = FALSE)
  }

  # Precompute constants
  inv_scale <- 1 / scale
  a <- 1 + 1 / shape
  lgamma_a <- lgamma(a)

  # E[T] F~_T(t) = scale * gamma(a, x) with x = (t / scale)^shape, formed on
  # the log scale from the regularised pgamma to avoid overflow
  terms_fn <- function(t) {
    x <- (t * inv_scale)^shape
    t * -expm1(-x) -
      scale * exp(pgamma(x, a, log.p = TRUE) + lgamma_a)
  }
  upper_fn <- function(t) {
    x <- (t * inv_scale)^shape
    t * exp(-x) -
      scale * exp(
        pgamma(x, a, lower.tail = FALSE, log.p = TRUE) + lgamma_a
      )
  }
  delay_cdf <- function(t) -expm1(-(t * inv_scale)^shape)

  .pcens_cdf_uniform(
    q, pwindow, terms_fn, upper_fn, scale * exp(lgamma_a), delay_cdf
  )
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
    return(
      pcens_cdf.default(object, q, pwindow, use_numeric)
    )
  }

  # Extract generalised gamma (Stacy) distribution parameters
  shape <- object$args$shape
  scale <- object$args$scale
  k <- object$args$k
  if (is.null(shape)) {
    stop(
      "shape parameter is required for generalised gamma distribution",
      call. = FALSE
    )
  }
  if (is.null(scale)) {
    stop(
      "scale parameter is required for generalised gamma distribution",
      call. = FALSE
    )
  }
  if (is.null(k)) {
    stop(
      "k parameter is required for generalised gamma distribution",
      call. = FALSE
    )
  }

  return(.pcens_cdf_gengamma_unif(q, pwindow, shape, scale, k))
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
  # flexsurv defaults for mu and sigma
  mu <- object$args$mu
  sigma <- object$args$sigma
  Q <- object$args$Q
  if (is.null(mu)) {
    mu <- 0
  }
  if (is.null(sigma)) {
    sigma <- 1
  }
  if (isTRUE(use_numeric) || is.null(Q) || Q <= 0) {
    return(
      pcens_cdf.default(object, q, pwindow, use_numeric)
    )
  }

  return(.pcens_cdf_gengamma_unif(
    q, pwindow,
    shape = Q / sigma, scale = exp(mu) * Q^(2 * sigma / Q), k = Q^-2
  ))
}

#' Analytical primary event censored CDF for the generalised gamma
#'
#' Shared implementation for the Stacy parameterisation used by both
#' generalised gamma [pcens_cdf()] methods.
#'
#' @inheritParams pcens_cdf
#'
#' @param shape,scale,k Generalised gamma parameters in the Stacy
#'  parameterisation of `flexsurv::pgengamma.orig()`.
#'
#' @inherit pcens_cdf return
#'
#' @keywords internal
.pcens_cdf_gengamma_unif <- function(q, pwindow, shape, scale, k) {
  # The partial expectation distribution has shape k + 1 / shape
  k_shift <- k + 1 / shape
  E_T <- scale * exp(lgamma(k_shift) - lgamma(k))

  terms_fn <- function(t) {
    x <- (t / scale)^shape
    t * pgamma(x, k) - E_T * pgamma(x, k_shift)
  }
  upper_fn <- function(t) {
    x <- (t / scale)^shape
    t * pgamma(x, k, lower.tail = FALSE) -
      E_T * pgamma(x, k_shift, lower.tail = FALSE)
  }
  delay_cdf <- function(t) pgamma((t / scale)^shape, k)

  .pcens_cdf_uniform(q, pwindow, terms_fn, upper_fn, E_T, delay_cdf)
}
