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
#' `q`, `pwindow` and the elements of `params` are recycled against each
#' other.
#'
#' @inheritParams pcens_cdf
#'
#' @param terms_fn Function of delays greater than or equal to zero and
#'  `params` returning \eqn{G}.
#'
#' @param upper_fn Function of delays greater than or equal to zero and
#'  `params` returning \eqn{H(t) = t S(t) - E \tilde S(t)}, where \eqn{S} and
#'  \eqn{\tilde S} are the upper tails of \eqn{F} and \eqn{\tilde F}.
#'
#' @param delay_cdf Function of delays and `params` returning the delay CDF.
#'
#' @param params Named list of delay parameters, including the delay mean
#'  \eqn{E} as `mean`.
#'
#' @return Numeric vector of CDF values in \[0, 1\].
#'
#' @keywords internal
.pcens_cdf_uniform <- function(q, pwindow, terms_fn, upper_fn, delay_cdf,
                               params) {
  .check_pwindow(q, pwindow)
  if (length(q) == 0L) {
    return(numeric(0))
  }
  param_lengths <- lengths(params)
  n <- max(length(q), length(pwindow), param_lengths)
  if (any(n %% c(length(q), length(pwindow), param_lengths) != 0L)) {
    warning(
      "longer object length is not a multiple of shorter object length",
      call. = FALSE
    )
  }
  scalar_params <- all(param_lengths == 1L)
  share <- scalar_params && length(pwindow) == 1L
  q <- rep_len(q, n)
  pwindow <- rep_len(pwindow, n)
  if (!scalar_params) {
    params <- lapply(params, rep_len, n)
  }
  active <- q > 0 & q < Inf
  exact <- active & pwindow == 0
  if (!any(exact) && all(active)) {
    result <- .uniform_window_cdf(
      q, pwindow, terms_fn, upper_fn, delay_cdf, params, share
    )
    return(pmin.int(1, pmax.int(0, result)))
  }

  result <- numeric(n)
  result[q == Inf] <- 1
  if (any(exact)) {
    result[exact] <- delay_cdf(q[exact], .index_params(params, exact))
    active <- active & !exact
  }
  if (any(active)) {
    result[active] <- .uniform_window_cdf(
      q[active], pwindow[active], terms_fn, upper_fn, delay_cdf,
      .index_params(params, active), share
    )
  }
  pmin.int(1, pmax.int(0, result))
}

#' Check the delays and windows of a uniform primary CDF
#'
#' @inheritParams pcens_cdf
#'
#' @return Invisible NULL if valid, otherwise stops with an error message.
#'
#' @keywords internal
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

#' Index vector delay parameters
#'
#' @inheritParams .pcens_cdf_uniform
#'
#' @param i Logical or integer index.
#'
#' @return `params` with each element longer than one indexed by `i`.
#'
#' @keywords internal
.index_params <- function(params, i) {
  lapply(params, function(x) if (length(x) == 1L) x else x[i])
}

# Delays more than this many windows wide may use the survival form
.upper_tail_ratio <- 1e3
# Delays more than this many windows wide use the mean of the delay CDF
.narrow_window_ratio <- 1e6

#' Uniform primary CDF for positive finite delays and windows
#'
#' @inheritParams .pcens_cdf_uniform
#'
#' @param d Positive finite delays.
#'
#' @param w Primary event windows greater than zero, one per delay.
#'
#' @param share Whether terms may be shared between window ends, as in
#'  [.uniform_endpoint_diff()]. Only valid for a scalar window and scalar
#'  `params`.
#'
#' @return Numeric vector of CDF values, not clamped to \[0, 1\].
#'
#' @keywords internal
.uniform_window_cdf <- function(d, w, terms_fn, upper_fn, delay_cdf, params,
                                share) {
  # Both forms lose about 1e-16 * d / w for narrow windows (1e-10 at the
  # cutoff), so the mean of the delay CDF over the window is used instead
  narrow <- d > .narrow_window_ratio * w
  if (any(narrow)) {
    if (all(narrow)) {
      return(.narrow_window_cdf(d, w, delay_cdf, params))
    }
    result <- numeric(length(d))
    result[narrow] <- .narrow_window_cdf(
      d[narrow], w[narrow], delay_cdf, .index_params(params, narrow)
    )
    result[!narrow] <- .uniform_window_cdf(
      d[!narrow], w[!narrow], terms_fn, upper_fn, delay_cdf,
      .index_params(params, !narrow), share
    )
    return(result)
  }
  lo <- pmax.int(d - w, 0)
  # In the far upper tail G(d) - G(lo) cancels, unlike 1 - (H(d) - H(lo)) / w
  # with G(t) = t - mean - H(t). This only pays off past the delay mean.
  up <- lo > params$mean & d > .upper_tail_ratio * w
  if (!any(up)) {
    return(.uniform_endpoint_diff(terms_fn, d, lo, params, share) / w)
  }
  if (all(up)) {
    return(1 - .uniform_endpoint_diff(upper_fn, d, lo, params, share) / w)
  }
  result <- numeric(length(d))
  result[up] <- 1 - .uniform_endpoint_diff(
    upper_fn, d[up], lo[up], .index_params(params, up), share
  ) / w[up]
  result[!up] <- .uniform_endpoint_diff(
    terms_fn, d[!up], lo[!up], .index_params(params, !up), share
  ) / w[!up]
  result
}

#' Difference of terms between the ends of uniform primary windows
#'
#' Returns `fn(d, params) - fn(lo, params)`.
#' With `share = TRUE` and at least 64 delays, lower ends that equal a delay
#' reuse the terms at that delay, so `fn` is evaluated once per distinct
#' window end.
#' This is tried only when at least 3 of 6 sampled lower ends match a delay,
#' which skips the full match when few ends coincide, and used only when at
#' most half as many unmatched distinct lower ends as delays remain.
#'
#' @inheritParams .uniform_window_cdf
#'
#' @param fn Function of delays and `params`, such as `terms_fn` or
#'  `upper_fn`.
#'
#' @param lo Lower window ends, one per delay.
#'
#' @return Numeric vector of differences, one per delay.
#'
#' @keywords internal
.uniform_endpoint_diff <- function(fn, d, lo, params, share) {
  n <- length(d)
  # Matching costs more than it saves for short inputs
  if (share && n >= 64L) {
    probe <- round(seq.int(1, n, length.out = 6L))
    if (sum(lo[probe] %in% d) >= 3L) {
      idx <- match(lo, d)
      unmatched <- is.na(idx)
      extra <- unique(lo[unmatched])
      if (2L * length(extra) <= n) {
        idx[unmatched] <- n + match(lo[unmatched], extra)
        at_d <- fn(d, params)
        return(at_d - c(at_d, fn(extra, params))[idx])
      }
    }
  }
  fn(d, params) - fn(lo, params)
}

#' Mean of the delay CDF over narrow uniform primary windows
#'
#' Uses a 5 point Gauss-Legendre rule over \[d - w, d\].
#'
#' @inheritParams .uniform_window_cdf
#'
#' @return Numeric vector of CDF values, one per delay.
#'
#' @keywords internal
.narrow_window_cdf <- function(d, w, delay_cdf, params) {
  # Weights are halved so they sum to 1
  nodes <- c(-0.906179845938664, -0.5384693101056831, 0,
             0.5384693101056831, 0.906179845938664)
  gl_weights <- c(0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
                  0.4786286704993665, 0.2369268850561891) / 2
  n <- length(d)
  half <- w / 2
  x <- rep(d - half, 5) + rep(nodes, each = n) * rep(half, 5)
  cdf <- delay_cdf(x, .index_params(params, rep(seq_len(n), 5)))
  drop(matrix(cdf, nrow = n) %*% gl_weights)
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
  terms_fn <- function(t, p) {
    t * pgamma(t, p$shape, scale = p$scale) -
      p$mean * pgamma(t, p$shape + 1, scale = p$scale)
  }
  upper_fn <- function(t, p) {
    t * pgamma(t, p$shape, scale = p$scale, lower.tail = FALSE) -
      p$mean * pgamma(t, p$shape + 1, scale = p$scale, lower.tail = FALSE)
  }
  delay_cdf <- function(t, p) pgamma(t, p$shape, scale = p$scale)

  .pcens_cdf_uniform(
    q, pwindow, terms_fn, upper_fn, delay_cdf,
    list(shape = shape, scale = scale, mean = shape * scale)
  )
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
  terms_fn <- function(t, p) {
    z <- (log(t) - p$mu) / p$sigma
    t * stats::pnorm(z) - p$mean * stats::pnorm(z - p$sigma)
  }
  upper_fn <- function(t, p) {
    z <- (log(t) - p$mu) / p$sigma
    t * stats::pnorm(z, lower.tail = FALSE) -
      p$mean * stats::pnorm(z - p$sigma, lower.tail = FALSE)
  }
  delay_cdf <- function(t, p) stats::plnorm(t, p$mu, p$sigma)

  .pcens_cdf_uniform(
    q, pwindow, terms_fn, upper_fn, delay_cdf,
    list(mu = mu, sigma = sigma, mean = exp(mu + 0.5 * sigma^2))
  )
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

  # E[T] F~_T(t) = scale * gamma(a, x) with a = 1 + 1 / shape and
  # x = (t / scale)^shape, formed on the log scale from the regularised
  # pgamma to avoid overflow
  terms_fn <- function(t, p) {
    x <- (t / p$scale)^p$shape
    t * -expm1(-x) -
      p$scale * exp(pgamma(x, p$a, log.p = TRUE) + p$lgamma_a)
  }
  upper_fn <- function(t, p) {
    x <- (t / p$scale)^p$shape
    t * exp(-x) -
      p$scale * exp(
        pgamma(x, p$a, lower.tail = FALSE, log.p = TRUE) + p$lgamma_a
      )
  }
  delay_cdf <- function(t, p) -expm1(-(t / p$scale)^p$shape)

  a <- 1 + 1 / shape
  lgamma_a <- lgamma(a)
  .pcens_cdf_uniform(
    q, pwindow, terms_fn, upper_fn, delay_cdf,
    list(
      shape = shape, scale = scale, a = a, lgamma_a = lgamma_a,
      mean = scale * exp(lgamma_a)
    )
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
  terms_fn <- function(t, p) {
    x <- (t / p$scale)^p$shape
    t * pgamma(x, p$k) - p$mean * pgamma(x, p$k_shift)
  }
  upper_fn <- function(t, p) {
    x <- (t / p$scale)^p$shape
    t * pgamma(x, p$k, lower.tail = FALSE) -
      p$mean * pgamma(x, p$k_shift, lower.tail = FALSE)
  }
  delay_cdf <- function(t, p) pgamma((t / p$scale)^p$shape, p$k)

  k_shift <- k + 1 / shape
  .pcens_cdf_uniform(
    q, pwindow, terms_fn, upper_fn, delay_cdf,
    list(
      shape = shape, scale = scale, k = k, k_shift = k_shift,
      mean = scale * exp(lgamma(k_shift) - lgamma(k))
    )
  )
}

#' Method for Exponential delay with uniform primary
#'
#' Analytical solution for the exponential distribution, which is the gamma
#' solution with `shape = 1` without the incomplete gamma function.
#' Delay arguments other than `rate`, such as `lower.tail`, use the numerical
#' [pcens_cdf.default()] method.
#' See `vignette("analytic-solutions")` for the derivation.
#'
#' @inheritParams pcens_cdf
#'
#' @family pcens
#'
#' @inherit pcens_cdf return
#'
#' @export
#' @examples
#' pcens_obj <- new_pcens(
#'   pdist = pexp,
#'   dprimary = dunif,
#'   primary_args = list(min = 0, max = 1),
#'   rate = 0.5
#' )
#' pcens_cdf(pcens_obj, q = c(1, 4, 8), pwindow = 1)
pcens_cdf.pcens_pexp_dunif <- function(
  object,
  q,
  pwindow,
  use_numeric = FALSE
) {
  delay_args <- .delay_args(object, "rate", use_numeric)
  if (is.null(delay_args)) {
    return(pcens_cdf.default(object, q, pwindow, use_numeric))
  }
  rate <- if (is.null(delay_args$rate)) 1 else delay_args$rate

  # G(t) = 0 for t <= 0 as F_T(t) = 0 there, so clamp t at 0
  .pcens_cdf_antiderivative(q, pwindow, function(t) {
    .expon_shortfall(rate * pmax(t, 0)) / rate
  })
}

#' Method for Normal delay with uniform primary
#'
#' Analytical solution for the normal distribution, which has support on
#' the reals.
#' The primary event window \eqn{[d - w_P, d]} is not clipped at zero.
#' Below \eqn{z = -10} an asymptotic series avoids cancellation.
#' Delay arguments other than `mean` and `sd`, such as `lower.tail`, use the
#' numerical [pcens_cdf.default()] method.
#' See `vignette("analytic-solutions")` for the derivation.
#'
#' @inheritParams pcens_cdf
#'
#' @family pcens
#'
#' @inherit pcens_cdf return
#'
#' @export
#' @examples
#' pcens_obj <- new_pcens(
#'   pdist = pnorm,
#'   dprimary = dunif,
#'   primary_args = list(min = 0, max = 1),
#'   mean = 5,
#'   sd = 2
#' )
#' pcens_cdf(pcens_obj, q = c(-1, 2, 5, 8), pwindow = 1)
pcens_cdf.pcens_pnorm_dunif <- function(
  object,
  q,
  pwindow,
  use_numeric = FALSE
) {
  delay_args <- .delay_args(object, c("mean", "sd"), use_numeric)
  if (is.null(delay_args)) {
    return(pcens_cdf.default(object, q, pwindow, use_numeric))
  }
  mu <- if (is.null(delay_args$mean)) 0 else delay_args$mean
  sigma <- if (is.null(delay_args$sd)) 1 else delay_args$sd

  .pcens_cdf_antiderivative(q, pwindow, function(t) {
    sigma * .norm_shortfall((t - mu) / sigma)
  })
}

#' Method for Chi-square delay with uniform primary
#'
#' The chi-square distribution with `df` degrees of freedom is the gamma
#' distribution with `shape = df / 2` and `scale = 2`, so this uses
#' [pcens_cdf.pcens_pgamma_dunif()].
#' A non-central chi-square (`ncp` not zero) has no such form and uses the
#' numerical [pcens_cdf.default()] method.
#' So do delay arguments other than `df` and `ncp`, such as `lower.tail`.
#'
#' @inheritParams pcens_cdf
#'
#' @family pcens
#'
#' @inherit pcens_cdf return
#'
#' @export
#' @examples
#' pcens_obj <- new_pcens(
#'   pdist = pchisq,
#'   dprimary = dunif,
#'   primary_args = list(min = 0, max = 1),
#'   df = 4
#' )
#' pcens_cdf(pcens_obj, q = c(1, 4, 8), pwindow = 1)
pcens_cdf.pcens_pchisq_dunif <- function(
  object,
  q,
  pwindow,
  use_numeric = FALSE
) {
  delay_args <- .delay_args(object, c("df", "ncp"), use_numeric)
  if (is.null(delay_args)) {
    return(pcens_cdf.default(object, q, pwindow, use_numeric))
  }
  if (is.null(delay_args$df)) {
    stop("df parameter is required for Chi-square distribution", call. = FALSE)
  }

  object$pdist <- pgamma
  object$args <- list(shape = delay_args$df / 2, scale = 2)
  pcens_cdf.pcens_pgamma_dunif(object, q, pwindow)
}

#' Method for Beta delay with uniform primary
#'
#' Analytical solution for the beta distribution, which has support on
#' \eqn{[0, 1]}.
#' A non-central beta (`ncp` not zero) uses the numerical
#' [pcens_cdf.default()] method.
#' So do delay arguments other than `shape1`, `shape2` and `ncp`, such as
#' `lower.tail`.
#'
#' @inheritParams pcens_cdf
#'
#' @family pcens
#'
#' @inherit pcens_cdf return
#'
#' @export
#' @examples
#' pcens_obj <- new_pcens(
#'   pdist = pbeta,
#'   dprimary = dunif,
#'   primary_args = list(min = 0, max = 1),
#'   shape1 = 2,
#'   shape2 = 3
#' )
#' pcens_cdf(pcens_obj, q = c(0.2, 0.6, 1.5), pwindow = 0.5)
pcens_cdf.pcens_pbeta_dunif <- function(
  object,
  q,
  pwindow,
  use_numeric = FALSE
) {
  delay_args <- .delay_args(
    object, c("shape1", "shape2", "ncp"), use_numeric
  )
  if (is.null(delay_args)) {
    return(pcens_cdf.default(object, q, pwindow, use_numeric))
  }
  a <- delay_args$shape1
  b <- delay_args$shape2
  if (is.null(a) || is.null(b)) {
    stop(
      "shape1 and shape2 parameters are required for Beta distribution",
      call. = FALSE
    )
  }
  E_T <- a / (a + b)

  # Clamp t to [0, 1] and add the linear part above the support
  .pcens_cdf_antiderivative(q, pwindow, function(t) {
    t_in <- pmin(pmax(t, 0), 1)
    t_in * stats::pbeta(t_in, a, b) -
      E_T * stats::pbeta(t_in, a + 1, b) + pmax(t - 1, 0)
  })
}

#' Delay distribution arguments an analytical solution handles
#'
#' Arguments the solution does not use, such as `lower.tail`, would be
#' ignored, and a non-zero `ncp` changes the distribution, so the numerical
#' [pcens_cdf.default()] method is used instead.
#'
#' @inheritParams pcens_cdf
#'
#' @param allowed Character vector of the argument names the solution uses.
#'
#' @return The named list of delay arguments, or `NULL` if `use_numeric` is
#'  `TRUE`, any argument is unnamed or not in `allowed`, or `ncp` is not
#'  zero.
#'
#' @keywords internal
.delay_args <- function(object, allowed, use_numeric = FALSE) {
  if (isTRUE(use_numeric)) {
    return(NULL)
  }
  delay_args <- object$args
  if (length(delay_args) == 0L) {
    return(delay_args)
  }
  if (
    is.null(names(delay_args)) || !all(names(delay_args) %in% allowed) ||
      !isTRUE(all(delay_args[["ncp"]] == 0))
  ) {
    return(NULL)
  }
  delay_args
}

#' Primary event censored CDF from an antiderivative of the delay CDF
#'
#' The CDF is \eqn{(G(d) - G(d - w_P)) / w_P} for an antiderivative
#' \eqn{G} of the delay CDF.
#'
#' @inheritParams pcens_cdf
#'
#' @param G Function giving an antiderivative of the delay CDF at a vector
#'  of times.
#'
#' @inherit pcens_cdf return
#'
#' @keywords internal
.pcens_cdf_antiderivative <- function(q, pwindow, G) {
  result <- (G(q) - G(q - pwindow)) / pwindow
  # Both antiderivatives are infinite at q = Inf, and the CDF is 1
  result[which(q == Inf)] <- 1

  # Ensure the result is in [0, 1] (accounts for numerical errors)
  pmin(1, pmax(0, result))
}

#' Evaluate x - 1 + exp(-x) for x >= 0
#'
#' A series is used for `x < 0.1` where the direct form cancels.
#'
#' @param x Non-negative numeric vector.
#'
#' @return Vector of `x - 1 + exp(-x)`.
#'
#' @keywords internal
.expon_shortfall <- function(x) {
  out <- x + expm1(-x)
  small <- which(x < 0.1)
  if (length(small) > 0L) {
    # x^2 / 2 * sum_k 2 (-x)^k / (k + 2)!
    xs <- x[small]
    term <- 1
    series <- 1
    for (k in 1:10) {
      term <- term * -xs / (k + 2)
      series <- series + term
    }
    out[small] <- xs^2 / 2 * series
  }
  out
}

#' Expected shortfall of a standard normal
#'
#' Evaluates \eqn{g(z) = z \Phi(z) + \phi(z)}.
#' Below `z = -10` the asymptotic series
#' \eqn{\phi(z) / z^2 \sum_n (-1)^n (2n + 1)!! / z^{2n}} is used with 20
#' terms, as the direct form cancels.
#'
#' @param z Numeric vector of standardised values.
#'
#' @return Vector of `g(z)`.
#'
#' @keywords internal
.norm_shortfall <- function(z) {
  out <- stats::dnorm(z) + z * stats::pnorm(z)
  far <- which(z < -10)
  if (length(far) > 0L) {
    zf <- z[far]
    y <- 1 / zf^2
    term <- 1
    series <- 1
    for (k in 1:20) {
      term <- term * -(2 * k + 1) * y
      series <- series + term
    }
    out[far] <- stats::dnorm(zf) * y * series
  }
  out
}
