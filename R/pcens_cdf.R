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

  partial_pgamma <- function(q) {
    pgamma(q, shape = shape, scale = scale)
  }
  partial_pgamm_k_1 <- function(q) {
    pgamma(q, shape = shape + 1, scale = scale)
  }
  # Adjust q so that we have [q-pwindow, q]
  q <- q - pwindow
  # Handle cases where q + pwindow <= 0
  zero_cases <- q + pwindow <= 0
  result <- ifelse(zero_cases, 0, NA)

  # Process non-zero cases only if there are any
  if (!all(zero_cases)) {
    non_zero_q <- q[!zero_cases]
    d <- non_zero_q + pwindow

    # Compute delay CDF at the interval endpoints and at the shifted (k+1)
    # distribution for the mean-shift term E[T] = shape * scale.
    F_T_q <- partial_pgamma(non_zero_q)
    F_T_d <- partial_pgamma(d)
    F_T_q_kp1 <- partial_pgamm_k_1(non_zero_q)
    F_T_d_kp1 <- partial_pgamm_k_1(d)

    E_T <- shape * scale

    # Direct CDF form:
    #   F_{S+}(d) = ( d F_T(d) - q F_T(q) - E_T (F~_T(d) - F~_T(q)) ) / w_P
    non_zero_result <-
      (d * F_T_d - non_zero_q * F_T_q -
        E_T * (F_T_d_kp1 - F_T_q_kp1)) / pwindow

    # Assign non-zero results back to the main result vector
    result[!zero_cases] <- non_zero_result
  }

  # Ensure the result is in [0, 1] (accounts for numerical errors)
  result <- pmin(1, pmax(0, result))

  return(result)
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

  partial_plnorm <- function(q) {
    stats::plnorm(q, meanlog = mu, sdlog = sigma)
  }
  partial_plnorm_sigma2 <- function(q) {
    stats::plnorm(q, meanlog = mu + sigma^2, sdlog = sigma)
  }
  # Adjust q so that we have [q-pwindow, q]
  q <- q - pwindow

  # Handle cases where q + pwindow <= 0
  zero_cases <- q + pwindow <= 0
  result <- ifelse(zero_cases, 0, NA)

  # Process non-zero cases only if there are any
  if (!all(zero_cases)) {
    non_zero_q <- q[!zero_cases]
    d <- non_zero_q + pwindow

    # Compute delay CDF at the interval endpoints and at the shifted
    # (meanlog + sigma^2) distribution for the mean-shift term
    # E[T] = exp(mu + sigma^2 / 2).
    F_T_q <- partial_plnorm(non_zero_q)
    F_T_d <- partial_plnorm(d)
    F_T_q_shift <- partial_plnorm_sigma2(non_zero_q)
    F_T_d_shift <- partial_plnorm_sigma2(d)

    E_T <- exp(mu + 0.5 * sigma^2)

    # Direct CDF form:
    #   F_{S+}(d) = ( d F_T(d) - q F_T(q) - E_T (F~_T(d) - F~_T(q)) ) / w_P
    non_zero_result <-
      (d * F_T_d - non_zero_q * F_T_q -
        E_T * (F_T_d_shift - F_T_q_shift)) / pwindow

    # Assign non-zero results back to the main result vector
    result[!zero_cases] <- non_zero_result
  }

  # Ensure the result is in [0, 1] (accounts for numerical errors)
  result <- pmin(1, pmax(0, result))

  return(result)
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

  partial_pweibull <- function(q) {
    stats::pweibull(q, shape = shape, scale = scale)
  }

  # Precompute constants
  inv_shape <- 1 / shape
  inv_scale <- 1 / scale
  a <- 1 + inv_shape
  lgamma_a <- lgamma(a)

  # Lower incomplete gamma gamma(a, x) via the regularised form from
  # stats::pgamma, which is numerically stable for large x where an
  # unregularised series expansion would overflow.
  g <- function(t) {
    x <- (t * inv_scale)^shape
    exp(pgamma(x, shape = a, scale = 1, log.p = TRUE) + lgamma_a)
  }

  # Adjust q so that we have [q-pwindow, q]
  q <- q - pwindow

  # Handle cases where q + pwindow <= 0
  zero_cases <- q + pwindow <= 0
  result <- ifelse(zero_cases, 0, NA)

  # Process non-zero cases only if there are any
  if (!all(zero_cases)) {
    non_zero_q <- q[!zero_cases]
    d <- non_zero_q + pwindow
    # Clamp to zero for evaluating F_T and g (both undefined / zero on R_-).
    # The products q * F_T(q) and scale * g(q) are then zero when q < 0,
    # matching F_T(q) = 0 and g(q) = 0 for q <= 0.
    q_pos <- pmax(non_zero_q, 0)
    d_pos <- pmax(d, 0)

    # Compute delay CDF and helper g at the interval endpoints.
    F_T_q <- partial_pweibull(q_pos)
    F_T_d <- partial_pweibull(d_pos)
    g_q <- g(q_pos)
    g_d <- g(d_pos)

    # Direct CDF form (with E[T] = scale and the shifted-CDF role played
    # by g / scale, so that E[T] * (F~_T(d) - F~_T(q)) = scale * (g(d) - g(q))):
    #   F_{S+}(d) = ( d F_T(d) - q F_T(q) - scale (g(d) - g(q)) ) / w_P
    non_zero_result <-
      (d_pos * F_T_d - q_pos * F_T_q - scale * (g_d - g_q)) / pwindow

    # Assign non-zero results back to the main result vector
    result[!zero_cases] <- non_zero_result
  }

  # Ensure the result is in [0, 1] (accounts for numerical errors)
  result <- pmin(1, pmax(0, result))

  return(result)
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
  # F_T(t; k) = P(k, (t / scale)^shape) and the partial expectation
  # distribution is F_T(t; k + 1 / shape), both on the transformed scale.
  partial_pgengamma <- function(t, a) {
    pgamma((t / scale)^shape, shape = a)
  }
  k_shift <- k + 1 / shape

  # Adjust q so that we have [q-pwindow, q]
  q <- q - pwindow

  # Handle cases where q + pwindow <= 0
  zero_cases <- q + pwindow <= 0
  result <- ifelse(zero_cases, 0, NA)

  # Process non-zero cases only if there are any
  if (!all(zero_cases)) {
    non_zero_q <- q[!zero_cases]
    d <- non_zero_q + pwindow
    # Clamp to zero as F_T(t) = 0 for t <= 0 and (t / scale)^shape is
    # undefined for t < 0. The product q * F_T(q) is then zero when q < 0.
    q_pos <- pmax(non_zero_q, 0)
    d_pos <- pmax(d, 0)

    # Compute delay CDF at the interval endpoints and at the shifted
    # (k + 1 / shape) distribution for the mean-shift term
    # E[T] = scale * Gamma(k + 1 / shape) / Gamma(k).
    F_T_q <- partial_pgengamma(q_pos, k)
    F_T_d <- partial_pgengamma(d_pos, k)
    F_T_q_shift <- partial_pgengamma(q_pos, k_shift)
    F_T_d_shift <- partial_pgengamma(d_pos, k_shift)

    E_T <- scale * exp(lgamma(k_shift) - lgamma(k))

    # Direct CDF form:
    #   F_{S+}(d) = ( d F_T(d) - q F_T(q) - E_T (F~_T(d) - F~_T(q)) ) / w_P
    non_zero_result <-
      (d_pos * F_T_d - q_pos * F_T_q -
        E_T * (F_T_d_shift - F_T_q_shift)) / pwindow

    # Assign non-zero results back to the main result vector
    result[!zero_cases] <- non_zero_result
  }

  # Ensure the result is in [0, 1] (accounts for numerical errors)
  result <- pmin(1, pmax(0, result))

  return(result)
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
  delay_args <- .delay_args(object, "rate")
  if (isTRUE(use_numeric) || is.null(delay_args)) {
    return(
      pcens_cdf.default(object, q, pwindow, use_numeric)
    )
  }

  # pexp defaults to a rate of 1
  rate <- delay_args$rate
  if (is.null(rate)) {
    rate <- 1
  }

  # G(t) = 0 for t <= 0 as F_T(t) = 0 there, so clamp t at 0
  G <- function(t) {
    .expon_shortfall(rate * pmax(t, 0)) / rate
  }

  .pcens_cdf_antiderivative(q, pwindow, G)
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
  delay_args <- .delay_args(object, c("mean", "sd"))
  if (isTRUE(use_numeric) || is.null(delay_args)) {
    return(
      pcens_cdf.default(object, q, pwindow, use_numeric)
    )
  }

  # pnorm defaults to a standard normal
  mu <- delay_args$mean
  sigma <- delay_args$sd
  if (is.null(mu)) {
    mu <- 0
  }
  if (is.null(sigma)) {
    sigma <- 1
  }

  G <- function(t) {
    sigma * .norm_shortfall((t - mu) / sigma)
  }

  .pcens_cdf_antiderivative(q, pwindow, G)
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
  delay_args <- .delay_args(object, c("df", "ncp"))
  degrees <- delay_args$df
  if (
    isTRUE(use_numeric) || is.null(delay_args) ||
      .is_noncentral(delay_args$ncp)
  ) {
    return(
      pcens_cdf.default(object, q, pwindow, use_numeric)
    )
  }
  if (is.null(degrees)) {
    stop("df parameter is required for Chi-square distribution", call. = FALSE)
  }

  gamma_obj <- object
  gamma_obj$pdist <- pgamma
  gamma_obj$args <- list(shape = degrees / 2, scale = 2)
  pcens_cdf.pcens_pgamma_dunif(gamma_obj, q, pwindow)
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
  delay_args <- .delay_args(object, c("shape1", "shape2", "ncp"))
  a <- delay_args$shape1
  b <- delay_args$shape2
  if (
    isTRUE(use_numeric) || is.null(delay_args) ||
      .is_noncentral(delay_args$ncp)
  ) {
    return(
      pcens_cdf.default(object, q, pwindow, use_numeric)
    )
  }
  if (is.null(a)) {
    stop("shape1 parameter is required for Beta distribution", call. = FALSE)
  }
  if (is.null(b)) {
    stop("shape2 parameter is required for Beta distribution", call. = FALSE)
  }

  E_T <- a / (a + b)

  # Clamp t to [0, 1] and add the linear part above the support
  G <- function(t) {
    t_in <- pmin(pmax(t, 0), 1)
    t_in * stats::pbeta(t_in, a, b) -
      E_T * stats::pbeta(t_in, a + 1, b) + pmax(t - 1, 0)
  }

  .pcens_cdf_antiderivative(q, pwindow, G)
}

#' Delay distribution arguments an analytical solution handles
#'
#' Arguments the solution does not use, such as `lower.tail`, would be
#' ignored, so the numerical [pcens_cdf.default()] method is used instead.
#'
#' @param object A `pcens` object as created by [new_pcens()].
#'
#' @param allowed Character vector of the argument names the solution uses.
#'
#' @return The named list of delay arguments, or `NULL` if any argument is
#'  unnamed or not in `allowed`.
#'
#' @keywords internal
.delay_args <- function(object, allowed) {
  delay_args <- object$args
  arg_names <- names(delay_args)
  if (length(delay_args) == 0L) {
    return(delay_args)
  }
  if (is.null(arg_names) || !all(arg_names %in% allowed)) {
    return(NULL)
  }
  delay_args
}

#' Test for a non-central delay distribution
#'
#' @param ncp The `ncp` argument of a delay CDF, or `NULL`.
#'
#' @return `TRUE` if `ncp` is supplied and not zero.
#'
#' @keywords internal
.is_noncentral <- function(ncp) {
  !is.null(ncp) && !isTRUE(all(ncp == 0))
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
