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
#' \eqn{w}, the primary event censored CDF at \eqn{d} is
#' \eqn{(G(d) - G(q)) / w} with \eqn{q = \max(d - w, 0)} and
#' \eqn{G(t) = t F(t) - E \tilde F(t)}, where \eqn{F} is the delay CDF and
#' \eqn{\tilde F} is the CDF of the partial expectation distribution.
#' A delay distribution with an analytical solution of this form is added by
#' supplying `terms_fn`, which is called once for the delays and once for the
#' start of each window.
#'
#' Far into the upper tail \eqn{G(d)} and \eqn{G(q)} are both close to
#' \eqn{t - E}, so their difference loses all precision and can return 0 where
#' the answer is 1.
#' With \eqn{H(t) = t S(t) - E \tilde S(t)}, where \eqn{S = 1 - F} and
#' \eqn{\tilde S = 1 - \tilde F} are the upper tails, \eqn{G(t) = t - E - H(t)}
#' and for \eqn{d \ge w} the CDF is \eqn{1 - (H(d) - H(q)) / w}.
#' This removes the far tail cancellation.
#' The lower tail form has an absolute error of about
#' \eqn{10^{-16} d / w}, which is negligible unless \eqn{d} is many windows
#' wide.
#' The upper tail form also cancels when \eqn{H(q)} is not small, as for a
#' heavy tailed delay near the switch, so both forms lose about
#' \eqn{10^{-16} \max(|G|, |H|) / w} when the window is narrow relative to
#' \eqn{d}.
#' When `upper_fn`, which returns \eqn{H}, is supplied it is used for the
#' delays whose window starts above `switch_at`, taken as the delay mean, and
#' that are more than 1000 windows wide.
#' `terms_fn` is used for the rest.
#' The two forms agree where they are both accurate.
#'
#' Delays more than \eqn{10^6} windows wide are not given by either form.
#' If `delay_cdf` is supplied, the CDF is instead the mean of \eqn{F} over
#' the window, \eqn{(1 / w) \int_{d - w}^{d} F(u) du}, from a 5 point
#' Gauss-Legendre rule.
#' \eqn{F} is smooth over so narrow a window, so the rule is accurate to
#' rounding, and the cancellation error is at most about \eqn{10^{-10}}
#' for the forms used below the cutoff.
#' Without `delay_cdf` the two forms are always used.
#'
#' Delays at or below zero give 0, delays of `Inf` give 1 and missing delays
#' are an error.
#' Elements of `pwindow` equal to 0 are exact primary events and give the
#' delay CDF, from `delay_cdf`.
#' `q` and `pwindow` are recycled against each other, with a warning if the
#' longer length is not a multiple of the shorter.
#'
#' @inheritParams pcens_cdf
#'
#' @param terms_fn Function of a numeric vector of delays, all greater than or
#'  equal to zero, returning \eqn{G} at each.
#'  Sharing work between \eqn{F} and \eqn{\tilde F}, such as the standardised
#'  delay, is done inside `terms_fn`.
#'
#' @param upper_fn Optional function of a numeric vector of delays, all
#'  greater than or equal to zero, returning \eqn{H} at each, from the upper
#'  tails of the delay and partial expectation distributions.
#'
#' @param switch_at Delay above which the window start may use `upper_fn`, see
#'  Details.
#'  Only used if `upper_fn` is supplied.
#'
#' @param delay_cdf Optional function of a numeric vector of delays returning
#'  the delay CDF.
#'  Required if `pwindow` contains zeros, and used for delays more than
#'  \eqn{10^6} windows wide, see Details.
#'
#' @return Numeric vector of CDF values in \[0, 1\].
#'
#' @keywords internal
.pcens_cdf_uniform <- function(q, pwindow, terms_fn, upper_fn = NULL,
                               switch_at = Inf, delay_cdf = NULL) {
  if (anyNA(q)) {
    stop("q must not contain missing values.", call. = FALSE)
  }
  n_window <- length(pwindow)
  if (n_window == 0L || anyNA(pwindow) || min(pwindow) < 0) {
    stop(
      "pwindow must be non-negative with no missing values.",
      call. = FALSE
    )
  }
  vector_window <- n_window > 1L
  if (vector_window && length(q) != n_window) {
    if (length(q) == 0L) {
      return(numeric(0))
    }
    recycled <- .recycle_window(q, pwindow)
    q <- recycled$q
    pwindow <- recycled$pwindow
  }
  has_exact <- if (vector_window) any(pwindow == 0) else pwindow == 0
  active <- q > 0 & q < Inf
  if (!has_exact && all(active)) {
    # The common case, with every delay positive and finite
    result <- .uniform_window_cdf(
      q, pwindow, terms_fn, upper_fn, switch_at, delay_cdf
    )
    return(pmin.int(1, pmax.int(0, result)))
  }

  result <- numeric(length(q))
  result[q == Inf] <- 1
  if (has_exact) {
    if (is.null(delay_cdf)) {
      stop("delay_cdf is required when pwindow contains zeros.", call. = FALSE)
    }
    # Exact primary events are the delay CDF. The scalar case is dispatched
    # before this, see pcens_cdf().
    exact <- if (vector_window) active & pwindow == 0 else active
    if (any(exact)) {
      result[exact] <- delay_cdf(q[exact])
    }
    active <- active & !exact
  }
  if (any(active)) {
    w <- if (vector_window) pwindow[active] else pwindow
    result[active] <- .uniform_window_cdf(
      q[active], w, terms_fn, upper_fn, switch_at, delay_cdf
    )
  }
  # Ensure the result is in [0, 1] (accounts for numerical errors)
  pmin.int(1, pmax.int(0, result))
}

#' Recycle delays and a vector of windows against each other
#'
#' @inheritParams .pcens_cdf_uniform
#'
#' @return A list with `q` and `pwindow` of the same length.
#' A warning is given if the longer length is not a multiple of the shorter.
#'
#' @keywords internal
.recycle_window <- function(q, pwindow) {
  n <- max(length(q), length(pwindow))
  if (n %% length(q) != 0L || n %% length(pwindow) != 0L) {
    warning(
      "longer object length is not a multiple of shorter object length",
      call. = FALSE
    )
  }
  list(q = rep_len(q, n), pwindow = rep_len(pwindow, n))
}

#' Uniform primary CDF for positive, finite delays and positive windows
#'
#' Chooses between the lower tail, upper tail and narrow window forms of
#' `.pcens_cdf_uniform()` for each delay.
#'
#' @param d Numeric vector of finite delays greater than 0.
#'
#' @param w Window width, a single value or one per delay, greater than 0.
#'
#' @inheritParams .pcens_cdf_uniform
#'
#' @return Numeric vector of unclamped CDF values.
#'
#' @keywords internal
.uniform_window_cdf <- function(d, w, terms_fn, upper_fn, switch_at,
                                delay_cdf = NULL) {
  # Neither tail form is accurate for windows narrower than 1e-6 of the
  # delay, where they lose about 1e-16 * d / w, so use the mean of the delay
  # CDF over the window. This needs d > w, so the window starts above 0.
  narrow <- if (is.null(delay_cdf)) FALSE else d > 1e6 * w
  if (any(narrow)) {
    if (all(narrow)) {
      return(.narrow_window_cdf(d, w, delay_cdf))
    }
    vector_window <- length(w) > 1L
    result <- numeric(length(d))
    result[narrow] <- .narrow_window_cdf(
      d[narrow], if (vector_window) w[narrow] else w, delay_cdf
    )
    result[!narrow] <- .uniform_window_cdf(
      d[!narrow], if (vector_window) w[!narrow] else w,
      terms_fn, upper_fn, switch_at
    )
    return(result)
  }
  lo <- pmax.int(d - w, 0)
  # G(d) - G(lo) loses about d / w digits, so the upper tail form is only
  # needed for delays many windows wide, which keeps it off the usual path.
  up <- if (is.null(upper_fn)) {
    FALSE
  } else {
    lo > switch_at & d > 1e3 * w
  }
  if (!any(up)) {
    return((terms_fn(d) - terms_fn(lo)) / w)
  }
  if (all(up)) {
    return(1 - (upper_fn(d) - upper_fn(lo)) / w)
  }
  vector_window <- length(w) > 1L
  result <- numeric(length(d))
  result[up] <- 1 -
    (upper_fn(d[up]) - upper_fn(lo[up])) / (if (vector_window) w[up] else w)
  result[!up] <- (terms_fn(d[!up]) - terms_fn(lo[!up])) /
    (if (vector_window) w[!up] else w)
  result
}

#' Mean of the delay CDF over a narrow uniform primary window
#'
#' The uniform primary CDF is \eqn{(1 / w) \int_{d - w}^{d} F(u) du}.
#' This applies a 5 point Gauss-Legendre rule, which is exact for polynomials
#' of degree 9 and so accurate to rounding when the window is a small fraction
#' of the delay, see `.pcens_cdf_uniform()`.
#'
#' @param d Numeric vector of delays with `d > w`.
#'
#' @param w Window width, a single value or one per delay.
#'
#' @inheritParams .pcens_cdf_uniform
#'
#' @return Numeric vector of CDF values.
#'
#' @keywords internal
.narrow_window_cdf <- function(d, w, delay_cdf) {
  # Nodes and weights of the 5 point rule on [-1, 1], with weights scaled to
  # sum to 1 so the result is the mean of F over the window
  nodes <- c(-0.906179845938664, -0.5384693101056831, 0,
             0.5384693101056831, 0.906179845938664)
  gl_weights <- c(0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
               0.4786286704993665, 0.2369268850561891) / 2
  n <- length(d)
  half <- rep_len(w, n) / 2
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

  # G(t) = t F_T(t) - E_T F~_T(t), where F~_T is the CDF of the partial
  # expectation distribution, a gamma with shape + 1, and E_T = shape * scale.
  E_T <- shape * scale
  terms_fn <- function(t) {
    t * pgamma(t, shape, scale = scale) -
      E_T * pgamma(t, shape + 1, scale = scale)
  }
  # H(t) = t S(t) - E_T S~(t) from the upper tails, for the far upper tail
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

  # G(t) = t F_T(t) - E_T F~_T(t), where F~_T is the lognormal CDF shifted to
  # meanlog + sdlog^2 and E_T = exp(meanlog + sdlog^2 / 2). Both CDFs share
  # the standardised log delay z, as (log(t) - meanlog - sdlog^2) / sdlog
  # is z - sdlog.
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

  # G(t) = t F_T(t) - scale g(t), with E[T] F~_T(t) = scale g(t), where
  # g(t) = gamma(a, x) is the lower incomplete gamma function at
  # x = (t / scale)^shape and F_T(t) = 1 - exp(-x) shares x. g is formed on
  # the log scale from the regularised stats::pgamma, which is stable for
  # large x where an unregularised series expansion would overflow.
  terms_fn <- function(t) {
    x <- (t * inv_scale)^shape
    t * -expm1(-x) -
      scale * exp(pgamma(x, a, log.p = TRUE) + lgamma_a)
  }
  # H(t) = t S(t) - scale * Gamma(a, x) from the upper tails, for the far
  # upper tail
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
  # F_T(t; k) = P(k, x) with x = (t / scale)^shape, and the partial
  # expectation distribution is F_T(t; k + 1 / shape). Both share x.
  k_shift <- k + 1 / shape
  # E[T] = scale * Gamma(k + 1 / shape) / Gamma(k).
  E_T <- scale * exp(lgamma(k_shift) - lgamma(k))

  # G(t) = t F_T(t) - E_T F~_T(t)
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
