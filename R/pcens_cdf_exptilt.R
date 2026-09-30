#' Methods for delays with an exponentially tilted primary
#'
#' Analytical primary event censored CDFs for exponential, gamma, normal and
#' lognormal delay distributions with an exponentially tilted primary event
#' window, the [dexpgrowth()] primary distribution with `r` equal to the tilt
#' \eqn{\rho}.
#' They honour `use_numeric`, and use the numerical method of
#' [pcens_cdf.default()] when no closed form applies.
#'
#' @inheritParams pcens_cdf
#'
#' @details
#' With the window width \eqn{w}, delay CDF \eqn{F} and the transform
#' \eqn{J(x) = T_f(-\rho; x)} of [tilt_transform], the primary event censored
#' CDF at \eqn{q} is
#' \deqn{
#' F_\rho(q) = F(q - w) + \frac{e^{\rho q} \{J(q) - J(q - w)\} -
#'   \{F(q) - F(q - w)\}}{e^{\rho w} - 1},
#' }
#' which is the same as the expression of
#' \eqn{\{e^{\rho w} F(q - w) - F(q) + e^{\rho q} (J(q) - J(q - w))\} /
#' (e^{\rho w} - 1)}.
#' For delays on the non-negative reals \eqn{F(x) = J(x) = 0} for
#' \eqn{x \le 0}. The normal delay has full support and its transforms start
#' at \eqn{-\infty}.
#' All terms depend on one endpoint, \eqn{q} or \eqn{q - w}, so each
#' endpoint is evaluated once and reused when several `q` share an endpoint,
#' as for the integer delays of [pcens_pmf()].
#'
#' Closed forms of the transform, for \eqn{\xi = -\rho}:
#' * Exponential with rate \eqn{\lambda}:
#'   \eqn{T_f(\xi; \tau) = \lambda / (\lambda - \xi)
#'   \{1 - e^{-(\lambda - \xi) \tau}\}}.
#' * Gamma with shape \eqn{\alpha} and rate \eqn{\lambda}:
#'   \eqn{T_f(\xi; \tau) = (\lambda / (\lambda - \xi))^\alpha
#'   P(\alpha, (\lambda - \xi) \tau)} with \eqn{P} the regularised lower
#'   incomplete gamma function.
#' * Normal with mean \eqn{\mu} and standard deviation \eqn{\sigma}:
#'   \eqn{T_f(\xi; \tau) = e^{\xi \mu + \xi^2 \sigma^2 / 2}
#'   \Phi((\tau - \mu - \xi \sigma^2) / \sigma)}.
#'
#' * Lognormal with `meanlog` \eqn{\mu} and `sdlog` \eqn{\sigma}: no closed
#'   form. The transform is evaluated by quadrature for \eqn{\rho > 0} and by
#'   a series for \eqn{\rho < 0}, see [tilt_transform_lognormal].
#'
#' The terms are evaluated on the log scale, and differences are taken
#' between the lower or the upper tail representations, whichever loses less
#' precision, see `.log_diff_exp()`. Results agree with numerical integration
#' to a relative difference of about 1e-9 or better away from the two regimes
#' listed under **Precision**.
#'
#' **Admissibility.** The exponential and gamma forms need
#' \eqn{\lambda - \xi = \lambda + \rho > 0}, so that the tilted delay
#' distribution exists. Otherwise the method falls back to
#' [pcens_cdf.default()] for every `q`. The normal form has no restriction.
#' In Stan the fallback is the ODE, whose accuracy is lower in the lower tail
#' of a gamma with a shape below 1, where the density is singular at 0. For
#' shape 0.3 with \eqn{\rho = -1} and rate 1 the relative error of the ODE is
#' about 4e-2 at \eqn{q = 10^{-6}} and \eqn{10^{-3}}, against 2e-5 for the
#' numerical method in R. Stable forms for these tilts, such as
#' the expm1 closed form for the exponential and a Kummer series for the gamma,
#' are not implemented, see #388.
#' The lognormal has no tilted delay for \eqn{\rho < 0} and needs none, as
#' the transform is truncated at \eqn{t}. It falls back to the numerical
#' method for every `q` where \eqn{\rho \sigma^2 e^\mu} overflows, above
#' about \eqn{e^{690}}. For \eqn{\rho < 0} the series needs about
#' \eqn{|\rho| q + 9 \sqrt{|\rho| q} + 30} terms per quantile, so it is slower
#' than the numerical method beyond \eqn{|\rho| q} of about 200. The numerical
#' method is used for the `q` above that, unless \eqn{|\rho| w} is above 2,
#' where the numerical method loses accuracy in the lower tail (a relative
#' error of up to 5e-3 at \eqn{|\rho| w} of 600). The series is then kept up
#' to 20000 terms, which is for \eqn{|\rho| q} up to about 18700, and the
#' numerical method is used beyond that.
#'
#' **Tilts close to zero.** The expression above cancels as \eqn{\rho \to 0}.
#' Two forms replace it where the cancellation would lose precision. With
#' \eqn{G_k(t) = \int (t - u)^k f(u) du} over the support up to \eqn{t}:
#' * If \eqn{|\rho| w < 10^{-4}} and \eqn{|\rho| (|q| + w) < 0.1}, the
#'   uniform window limit with its first order correction in \eqn{\rho} is
#'   used,
#'   \deqn{F_\rho(q) = \frac{G_1(q) - G_1(q - w)}{w} +
#'     \rho \frac{G_2(q) - w G_1(q) - G_2(q - w) - w G_1(q - w)}{2 w} +
#'     O((\rho w)^2).}
#'   At \eqn{\rho = 0} this is the uniform window solution.
#' * For delays on the non-negative reals with \eqn{q < w}, where
#'   \eqn{F(q - w) = J(q - w) = 0}, and \eqn{|\rho| q < 10^{-4}}, the form
#'   \eqn{F_\rho(q) = \rho \{G_1(q) + \rho G_2(q) / 2\} /
#'   (e^{\rho w} - 1) + O((\rho q)^2)} is used. It keeps precision for `q`
#'   close to zero.
#'
#' Away from these regions the error of the direct form is below 1e-9. The
#' truncation error of the two forms is below 1e-9 at their thresholds.
#' The first form subtracts terms of the size of \eqn{|q| + w}, so its
#' rounding error grows as \eqn{10^{-14} (|q| + w) / w (1 + |\rho| (|q| +
#' w))}. The second factor makes it worse than the direct form, whose error
#' is about \eqn{10^{-15} / (|\rho| w)} whatever the distance from the
#' origin, once \eqn{|\rho| (|q| + w)} is above about 0.1. The direct form
#' is used there.
#'
#' **Precision.** The R CDF agrees with numerical integration (`integrate()`
#' at a relative tolerance of 1e-13) to a relative difference of about 1e-9 or
#' better for delays that are not deep in the tails. For tilts close to zero
#' the upper tail, \eqn{1 - F_\rho(q)}, has an absolute difference of up to
#' about \eqn{10^{-15} G_1(q) / w}, which is about 1e-14 for `pwindow = 1`
#' and up to about 1e-11 for `pwindow = 1e-3`.
#' The Stan CDF agrees to a relative difference of about 1e-9 or better for
#' windows of 0.5 or more, and of up to about 5e-9 for windows of 1e-3 with
#' delays up to 25.
#' The Stan gradients of the gamma solution are exact for any shape. The
#' incomplete gamma function is taken from a series below the shape and from
#' a continued fraction above it, as the shape derivatives of Stan's
#' `gamma_lcdf()` and `gamma_lccdf()` are inaccurate in parts of the bulk
#' and NaN for shapes of about 200 or more.
#' The normal solution takes the lower tail beyond 37 standard deviations from
#' an asymptotic series for the same reason.
#' Three regimes are less accurate.
#' * The small window form for a normal delay in the deep lower tail, about 30
#'   to 35 standard deviations below the mean, has a relative error of up to
#'   about 1e-7 in both R and Stan (3e-8 at 35 standard deviations for
#'   `pwindow = 1e-3`).
#' * A window much smaller than the delay loses precision to the difference of
#'   the terms at its two ends, up to about 1e-13 / `pwindow` in relative
#'   terms, as the uniform window solutions do. The first order correction of
#'   the small window form adds about \eqn{|\rho| q} times that. For a gamma
#'   delay with shape 2.5 and rate 0.4 and `pwindow = 1e-8` the relative error
#'   is 6e-7 for \eqn{\rho = 0.5} and 6e-5 for \eqn{\rho = 50}. Use
#'   `pwindow = 0` for a primary event at a known time.
#' * The small tilt forms have a Stan gradient in \eqn{\rho} that is less
#'   accurate than their value, which has a truncation error below 1e-9. The
#'   second order term of the expansion is not included. The relative error
#'   of the derivative in \eqn{\rho} is about \eqn{|\rho| w / 6} for the
#'   small window form, and of at most the same order in \eqn{|\rho| q} for
#'   the small delay form. It is up to about 2e-5 at the thresholds of
#'   \eqn{10^{-4}}. The gradients in the delay parameters have the error of
#'   the value.
#'
#' **Speed.** For many values of `q` the analytical method is faster than
#' `use_numeric = TRUE`, as each endpoint is evaluated once. It is about 2 to
#' 2.5 times faster for 10 values and about 10 to 17 times faster for 100
#' values. For a single `q` it is about 1.5 to 3 times slower, as the call
#' evaluates four transforms at two endpoints. It is kept for single values
#' as it is accurate where the numerical method has errors of up to 1e-2 in
#' the tails of a normal delay.
#'
#' **Lognormal precision and speed.** The quadrature is accurate to an
#' absolute difference of about 1e-11 in the log transform for `sdlog` up to
#' 1.8, and 1e-13 for `sdlog` of 1 or below. Each range of the quadrature is
#' split into `ceiling(sdlog / 1.8)` panels, which keeps the CDF accurate to
#' a relative difference of 1e-9 or better for `sdlog` up to 15, as tested.
#' The CDF agrees with numerical integration to a relative difference of 1e-7
#' or better over the tested grid of tilts, windows and quantiles. The
#' largest differences are deep in the lower tail, where the CDF is below
#' 1e-100 and the direct form cancels
#' by about \eqn{1 / (\rho q)} times the gap between `q` and the mean of the
#' delays below it.
#' The quadrature and the series have a fixed cost that the numerical method
#' beats for fewer than 10 quantiles, so `pcens_cdf()` uses the numerical
#' method for fewer than 10 `q` and the transform for 10 or more, which is
#' 1.2 times faster at 12 and about 3 times faster at 40 in a benchmark.
#' The numerical method is used for fewer than 10 `q` only where
#' \eqn{|\rho| w} is at most 1. Its relative error is about 1e-6 there, 1e-4
#' at 50, and it fails at 1000, where the transform is accurate to 1e-9, so
#' the transform is used for any number of `q` above that.
#' The Stan solution takes about as long as the ODE for one delay, and is 2 to
#' 4 times faster with the shared terms of the vectorised PMF, see the NEWS.
#'
#' **Extending.** A new delay distribution is supported by defining
#' `.pcens_tilt_lower()`, `.pcens_tilt_available()`, `.pcens_tilt_transform()`
#' and, for the small tilt forms, `.pcens_tilt_moments()` for its class, and
#' a `pcens_cdf` method for the class that calls `.pcens_cdf_exptilt()`. See
#' [tilt_transform]. Other primary event windows reuse the same transforms
#' with their own tilts and coefficients.
#'
#' @family pcens
#'
#' @inherit pcens_cdf return
#'
#' @name pcens_cdf_exptilt
#'
#' @examples
#' # Exponential delay with a growing primary event process
#' pexp_obj <- new_pcens(
#'   pdist = pexp, dprimary = dexpgrowth,
#'   primary_args = list(r = 0.3), rate = 0.5
#' )
#' pcens_cdf(pexp_obj, q = c(1, 4, 8), pwindow = 2)
#'
#' # Gamma delay with a declining process
#' pgamma_obj <- new_pcens(
#'   pdist = pgamma, dprimary = dexpgrowth,
#'   primary_args = list(r = -0.3), shape = 2, rate = 1
#' )
#' pcens_cdf(pgamma_obj, q = c(1, 4, 8), pwindow = 2)
#'
#' # Normal delay, for example a difference of event times
#' pnorm_obj <- new_pcens(
#'   pdist = pnorm, dprimary = dexpgrowth,
#'   primary_args = list(r = 0.1), mean = 3, sd = 2
#' )
#' pcens_cdf(pnorm_obj, q = c(-1, 3, 8), pwindow = 2)
NULL

#' @rdname pcens_cdf_exptilt
#' @export
pcens_cdf.pcens_pexp_dexpgrowth <- function(
  object,
  q,
  pwindow,
  use_numeric = FALSE
) {
  .pcens_cdf_exptilt(object, q, pwindow, use_numeric)
}

#' @rdname pcens_cdf_exptilt
#' @export
pcens_cdf.pcens_pgamma_dexpgrowth <- function(
  object,
  q,
  pwindow,
  use_numeric = FALSE
) {
  .pcens_cdf_exptilt(object, q, pwindow, use_numeric)
}

#' @rdname pcens_cdf_exptilt
#' @export
pcens_cdf.pcens_pnorm_dexpgrowth <- function(
  object,
  q,
  pwindow,
  use_numeric = FALSE
) {
  .pcens_cdf_exptilt(object, q, pwindow, use_numeric)
}

#' @rdname pcens_cdf_exptilt
#' @export
pcens_cdf.pcens_plnorm_dexpgrowth <- function(
  object,
  q,
  pwindow,
  use_numeric = FALSE
) {
  .pcens_cdf_exptilt(
    object, q, pwindow, use_numeric,
    min_q = .lnorm_exptilt_min_q, min_xw = .lnorm_exptilt_min_xw
  )
}

# Tilts with |rho| times the window (or q) below this use the small tilt
# forms. The direct form loses about 1e-14 / (|rho| w) of relative precision
# and the small tilt forms have a truncation error of about (|rho| w)^2 / 12,
# which balance at about 1e-4 with both below 1e-9.
.exptilt_small <- 1e-4

# The small window form cancels by about 1e-14 * (|q| + w) / w * (1 + |rho|
# (|q| + w)), and the direct form by about 1e-15 / (|rho| w). The small
# window form has the smaller error while |rho| (|q| + w) is below about 0.1,
# whatever the window. Beyond that the direct form is used.
.exptilt_small_reach <- 0.1

#' Tilt of the exponentially tilted primary of a pcens object
#'
#' @inheritParams pcens_cdf
#'
#' @return The tilt `r`, a single finite number.
#'
#' @keywords internal
.exptilt_rho <- function(object) {
  rho <- object$primary_args$r
  if (is.null(rho)) {
    stop(
      "r parameter is required for the exponential growth primary ",
      "distribution",
      call. = FALSE
    )
  }
  if (!is.numeric(rho) || length(rho) != 1L || !is.finite(rho)) {
    stop(
      "r must be a single finite number for the exponential growth primary ",
      "distribution",
      call. = FALSE
    )
  }
  rho
}

#' Primary event censored CDF for an exponentially tilted primary
#'
#' Shared implementation of the [pcens_cdf_exptilt] methods. It dispatches
#' on the delay class of `object` through the generics of [tilt_transform].
#'
#' @inheritParams pcens_cdf
#'
#' @param min_q Number of quantiles below which the numerical method is
#'   used, for delays whose transform has a fixed cost. The default 0 always
#'   uses the closed forms.
#'
#' @param min_xw Largest \eqn{|\rho| w} at which `min_q` applies. The
#'   numerical method loses accuracy for a larger tilt times the window, so
#'   the closed forms are used there for any number of quantiles.
#'
#' @return Vector of computed primary event censored CDFs.
#'
#' @keywords internal
.pcens_cdf_exptilt <- function(
  object,
  q,
  pwindow,
  use_numeric = FALSE,
  min_q = 0L,
  min_xw = Inf
) {
  if (isTRUE(use_numeric)) {
    return(pcens_cdf.default(object, q, pwindow, use_numeric))
  }
  rho <- .exptilt_rho(object)
  # The closed forms are for a single window and need the tilted delay. A
  # transform that is evaluated by quadrature has a fixed cost that the
  # numerical method beats for a few quantiles, where it is accurate.
  if (length(pwindow) != 1L || !is.finite(pwindow) || pwindow <= 0 ||
    !.pcens_tilt_available(object, -rho)) {
    return(pcens_cdf.default(object, q, pwindow, use_numeric))
  }
  if (length(q) < min_q && abs(rho) * pwindow <= min_xw) {
    return(pcens_cdf.default(object, q, pwindow, use_numeric))
  }

  result <- rep(NA_real_, length(q))
  result[!is.na(q) & q == Inf] <- 1
  result[!is.na(q) & q == -Inf] <- 0
  finite <- which(is.finite(q))
  if (length(finite) > 0L) {
    # Transforms that cannot be evaluated at a large q use the numerical
    # method for that q alone
    fits <- .pcens_tilt_fits(object, -rho, q[finite], pwindow)
    result[finite[fits]] <- .exptilt_cdf_finite(
      object, q[finite[fits]], pwindow, rho
    )
    if (!all(fits)) {
      result[finite[!fits]] <- pcens_cdf.default(
        object, q[finite[!fits]], pwindow, use_numeric
      )
    }
  }
  result
}

#' Exponentially tilted CDF at finite points
#'
#' Chooses the form for each `q` and evaluates the transforms at the unique
#' endpoints the forms need.
#'
#' @inheritParams pcens_cdf
#'
#' @param rho The tilt, the `r` of the exponential growth primary.
#'
#' @return Vector of CDFs, clamped to \[0, 1\].
#'
#' @keywords internal
.exptilt_cdf_finite <- function(object, q, pwindow, rho) {
  if (length(q) == 0L) {
    return(numeric(0))
  }
  lower <- .pcens_tilt_lower(object)
  positive <- is.finite(lower)
  log_cdf <- rep(-Inf, length(q))

  # Below the support of the delay no mass has arrived
  active <- !positive | q > lower
  if (abs(rho) * pwindow < .exptilt_small) {
    # The small window form cancels far from the origin, see
    # `.exptilt_small_reach`
    small_window <- active &
      abs(rho) * (abs(q) + pwindow) < .exptilt_small_reach
    tiny_delay <- rep(FALSE, length(q))
  } else {
    small_window <- rep(FALSE, length(q))
    tiny_delay <- positive & active & q < pwindow &
      abs(rho) * q < .exptilt_small
  }
  direct <- active & !small_window & !tiny_delay

  if (any(small_window)) {
    log_cdf[small_window] <- .exptilt_lcdf_small_window(
      object, q[small_window], pwindow, rho, lower
    )
  }
  if (any(tiny_delay)) {
    log_cdf[tiny_delay] <- .exptilt_lcdf_tiny_delay(
      object, q[tiny_delay], pwindow, rho
    )
  }
  if (any(direct)) {
    log_cdf[direct] <- .exptilt_lcdf_direct(
      object, q[direct], pwindow, rho, lower
    )
  }
  pmin(1, exp(log_cdf))
}

#' Endpoints at which an exponentially tilted CDF needs the transforms
#'
#' The CDF at `q` needs the terms at `q` and at `q - pwindow`. The union is
#' taken so each endpoint is evaluated once, even where neighbouring `q`
#' share one. For delays with support bounded below, all endpoints at or
#' below `lower` have the same terms and share one entry.
#'
#' @param q Numeric vector of finite quantiles.
#'
#' @param pwindow Primary event window.
#'
#' @param lower Lower end of the support of the delay, 0 or `-Inf`.
#'
#' @return Sorted numeric vector of unique endpoints. The first is `lower`
#'   if it is finite.
#'
#' @keywords internal
.exptilt_endpoints <- function(q, pwindow, lower) {
  endpoints <- sort.int(unique(c(q, q - pwindow)))
  if (is.finite(lower)) {
    endpoints <- c(lower, endpoints[endpoints > lower])
  }
  endpoints
}

#' Position of endpoints in the output of `.exptilt_endpoints()`
#'
#' @param t Numeric vector of endpoints.
#'
#' @inheritParams .exptilt_endpoints
#'
#' @param endpoints Output of `.exptilt_endpoints()`.
#'
#' @return Integer vector of positions.
#'
#' @keywords internal
.exptilt_index <- function(t, endpoints, lower) {
  if (is.finite(lower)) {
    t <- pmax(t, lower)
  }
  match(t, endpoints)
}

#' Direct form of the exponentially tilted log CDF
#'
#' Evaluates
#' \eqn{F(q - w) + \{e^{\rho q} (J(q) - J(q - w)) - (F(q) - F(q - w))\} /
#' (e^{\rho w} - 1)} on the log scale from the transforms at the endpoints.
#' Each of the differences \eqn{F(q) - F(q - w)} and
#' \eqn{J(q) - J(q - w)} is taken between the lower tail terms or between the
#' upper tail terms, whichever has the smaller ratio of the two terms.
#' This avoids cancellation in the upper tail, where for positive tilts
#' \eqn{e^{\rho q} J(q)} is much larger than the result.
#'
#' @inheritParams .exptilt_cdf_finite
#'
#' @param lower Lower end of the support of the delay.
#'
#' @return Vector of log CDFs.
#'
#' @keywords internal
.exptilt_lcdf_direct <- function(object, q, pwindow, rho, lower) {
  endpoints <- .exptilt_endpoints(q, pwindow, lower)
  endpoint_terms <- cbind(
    .pcens_tilt_pair(object, endpoints, 0),
    .pcens_tilt_pair(object, endpoints, -rho)
  )
  at_q <- endpoint_terms[
    .exptilt_index(q, endpoints, lower), ,
    drop = FALSE
  ]
  at_y <- endpoint_terms[
    .exptilt_index(q - pwindow, endpoints, lower), ,
    drop = FALSE
  ]
  log_diff_f <- .exptilt_tail_diff(
    at_q[, 1], at_y[, 1], at_q[, 2], at_y[, 2]
  )
  log_diff_j <- .exptilt_tail_diff(
    at_q[, 3], at_y[, 3], at_q[, 4], at_y[, 4]
  )
  if (rho > 0) {
    log_num <- .log_diff_exp(rho * q + log_diff_j, log_diff_f)
    log_den <- .log_diff_exp(rho * pwindow, 0)
  } else {
    log_num <- .log_diff_exp(log_diff_f, rho * q + log_diff_j)
    log_den <- .log1m_exp(rho * pwindow)
  }
  .log_sum_exp(at_y[, 1], log_num - log_den)
}

#' Log of a difference between two points from either tail
#'
#' For a lower tail quantity `L` and upper tail quantity `U` with
#' `L(t) + U(t)` constant, the difference `L(q) - L(y)` equals
#' `U(y) - U(q)`. The relative precision of a difference is best when the
#' subtracted term is small, so this takes the representation with the smaller
#' ratio of the terms.
#'
#' @param lower_q,lower_y Log lower tail quantity at `q` and at `y`.
#'
#' @param upper_q,upper_y Log upper tail quantity at `q` and at `y`.
#'
#' @return Vector of the log differences.
#'
#' @keywords internal
.exptilt_tail_diff <- function(lower_q, lower_y, upper_q, upper_y) {
  use_upper <- lower_y - lower_q > upper_q - upper_y
  # Terms that underflow on both sides give NaN, which is a zero difference.
  # An upper tail that diverges (`Inf`) is never used.
  use_upper[is.na(use_upper) | is.infinite(upper_q)] <- FALSE
  out <- .log_diff_exp(lower_q, lower_y)
  if (any(use_upper)) {
    out[use_upper] <- .log_diff_exp(upper_y[use_upper], upper_q[use_upper])
  }
  out
}

#' Small window form of the exponentially tilted log CDF
#'
#' The uniform window limit with its first order correction in the tilt,
#' used for \eqn{|\rho| w < 10^{-4}} and \eqn{|\rho| (|q| + w) < 0.1}.
#' It needs only the moments
#' \eqn{G_1} and \eqn{G_2} at the endpoints, see [pcens_cdf_exptilt].
#'
#' @inheritParams .exptilt_lcdf_direct
#'
#' @return Vector of log CDFs.
#'
#' @keywords internal
.exptilt_lcdf_small_window <- function(object, q, pwindow, rho, lower) {
  endpoints <- .exptilt_endpoints(q, pwindow, lower)
  moments <- .pcens_tilt_moments(object, endpoints)
  at_q <- moments[.exptilt_index(q, endpoints, lower), , drop = FALSE]
  at_y <- moments[.exptilt_index(q - pwindow, endpoints, lower), ,
    drop = FALSE
  ]
  # Scale by G_1(q) so nothing underflows when the CDF is small
  scale <- at_q[, 1]
  g1_y <- exp(at_y[, 1] - scale)
  g2_q <- exp(at_q[, 2] - scale)
  g2_y <- exp(at_y[, 2] - scale)
  relative <- (1 - g1_y) +
    0.5 * rho * (g2_q - pwindow - g2_y - pwindow * g1_y)
  out <- scale + log(pmax(relative, 0)) - log(pwindow)
  out[!is.finite(scale)] <- -Inf
  out
}

#' Small delay form of the exponentially tilted log CDF
#'
#' For delays on the non-negative reals with `q < pwindow` and
#' \eqn{|\rho| q < 10^{-4}}, where the terms at `q - pwindow` are zero.
#' See [pcens_cdf_exptilt].
#'
#' @inheritParams .exptilt_lcdf_direct
#'
#' @return Vector of log CDFs.
#'
#' @keywords internal
.exptilt_lcdf_tiny_delay <- function(object, q, pwindow, rho) {
  moments <- .pcens_tilt_moments(object, q)
  if (rho > 0) {
    log_den <- .log_diff_exp(rho * pwindow, 0)
  } else {
    log_den <- .log1m_exp(rho * pwindow)
  }
  moments[, 1] + log(abs(rho)) - log_den +
    log1p(0.5 * rho * exp(moments[, 2] - moments[, 1]))
}
