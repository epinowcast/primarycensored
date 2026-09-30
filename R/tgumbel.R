#' Truncated Gumbel distribution functions
#'
#' Density, distribution function, and random generation for a Gumbel
#' distribution truncated to `[min, max]`. It is a primary event distribution
#' for exposure that builds up or decays unevenly within the window, for
#' example environmental sources in outbreak reconstruction.
#'
#' @inheritParams expgrowth
#'
#' @param mu Location of the Gumbel distribution before truncation. It is on
#'  the same scale as `x`, so for the default `min = 0` it is the time from
#'  the start of the window at which the exposure rate peaks.
#'
#' @param beta Scale of the Gumbel distribution before truncation, positive.
#'
#' @return `dtgumbel` gives the density, `ptgumbel` gives the distribution
#' function, and `rtgumbel` generates random deviates.
#'
#' The length of the result is determined by `n` for `rtgumbel`, and is the
#' maximum of the lengths of the numerical arguments for the other functions.
#'
#' @details
#' With \eqn{G(z) = \exp[-\exp\{-(z - \mu) / \beta\}]}, the CDF of the
#' Gumbel distribution, the density on \[min, max\] is
#' \deqn{f(x) = \frac{G'(x)}{G(max) - G(min)}}
#' and the cumulative distribution function is
#' \deqn{F(x) = \frac{G(x) - G(min)}{G(max) - G(min)}.}
#' The density rises to a peak at \eqn{x = \mu} and falls more slowly
#' above it, so the window is skewed. If \eqn{\mu} is far below `min` the
#' density decays exponentially with rate \eqn{1 / \beta}, as for
#' [dexpgrowth()] with the rate \eqn{r = -1 / \beta}.
#' If it is far above `max` the density rises very steeply towards `max`.
#'
#' Write \eqn{s(x) = \exp\{-(x - \mu) / \beta\}}, so that
#' \eqn{G(x) = \exp\{-s(x)\}}. The CDF and its complement are evaluated as
#' \deqn{\log F(x) = \log(1 - e^{-\{s(min) - s(x)\}}) -
#'   \{s(x) - s(max)\} - \log(1 - e^{-\{s(min) - s(max)\}})}
#' and \deqn{\log(1 - F(x)) = \log(1 - e^{-\{s(x) - s(max)\}}) -
#'   \log(1 - e^{-\{s(min) - s(max)\}}),}
#' with the differences of \eqn{s} written as, for example,
#' \eqn{s(x) (e^{(x - min) / \beta} - 1)}.
#' This avoids the cancellation of \eqn{G(max) - G(min)} when
#' \eqn{\mu} is far from the window and keeps the tails precise.
#'
#' Random numbers are drawn by inverting the upper tail form, which adds
#' \eqn{s(x) - s(max)} to \eqn{s(max)} and so has no cancellation.
#'
#' The Stan equivalents are `tgumbel_lpdf()`, `tgumbel_lcdf()` and
#' `tgumbel_rng()`, with the primary distribution identifier 4 and
#' `primary_params = [mu, beta]`.
#' Analytical primary event censored CDFs for this primary are described in
#' [pcens_cdf_gumbel].
#'
#' @family primaryeventdistributions
#'
#' @examples
#' x <- seq(0, 2, by = 0.25)
#' dens <- dtgumbel(x, min = 0, max = 2, mu = 0.5, beta = 0.3)
#' cumprobs <- ptgumbel(x, min = 0, max = 2, mu = 0.5, beta = 0.3)
#' samples <- rtgumbel(100, min = 0, max = 2, mu = 0.5, beta = 0.3)
#'
#' @name tgumbel
NULL

#' Check the parameters of the truncated Gumbel distribution
#'
#' @inheritParams tgumbel
#'
#' @return `NULL` invisibly, called for its error.
#'
#' @keywords internal
.check_tgumbel <- function(min, max, mu, beta) {
  finite <- function(x) is.numeric(x) && !anyNA(x) && all(is.finite(x))
  if (!finite(mu) || !finite(beta) || any(beta <= 0)) {
    stop(
      "mu must be finite and beta must be finite and positive for the ",
      "truncated Gumbel distribution",
      call. = FALSE
    )
  }
  if (anyNA(min) || anyNA(max) || any(min >= max)) {
    stop(
      "min must be less than max for the truncated Gumbel distribution",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Log of `exp(x) - 1` for non-negative `x`
#'
#' @param x Numeric vector, at least 0.
#'
#' @return `log(expm1(x))`, `-Inf` at 0, without overflow for large `x`.
#'
#' @keywords internal
.log_expm1 <- function(x) {
  x + .log1m_exp(-x)
}

#' Log of the width of the Gumbel window on the scale of `s`
#'
#' With \eqn{s(x) = \exp\{-(x - \mu) / \beta\}} the truncated Gumbel CDF is
#' a ratio of \eqn{1 - e^{-\delta}} terms, where \eqn{\delta} is a difference
#' of \eqn{s} at two points. The window difference
#' \eqn{s(min) - s(max)} is computed as \eqn{s(max) \{e^{(max - min) /
#' \beta} - 1\}}, on the log scale.
#'
#' @inheritParams tgumbel
#'
#' @return The log of `s(min) - s(max)`.
#'
#' @keywords internal
.tgumbel_log_delta_window <- function(min, max, mu, beta) {
  -(max - mu) / beta + .log_expm1((max - min) / beta)
}

#' Log of the Gumbel differences of `s` at points in the window
#'
#' @inheritParams tgumbel
#'
#' @param x Numeric vector of points in `[min, max]`.
#'
#' @return A list with the log of `lower`, `s(min) - s(x)`, and `upper`,
#'   `s(x) - s(max)`, computed as in `.tgumbel_log_delta_window()`.
#'
#' @keywords internal
.tgumbel_log_deltas <- function(x, min, max, mu, beta) {
  list(
    lower = -(x - mu) / beta + .log_expm1((x - min) / beta),
    upper = -(max - mu) / beta + .log_expm1((max - x) / beta)
  )
}

#' @rdname tgumbel
#' @export
dtgumbel <- function(x, min = 0, max = 1, mu, beta, log = FALSE) {
  .check_tgumbel(min, max, mu, beta)
  log_s <- -(x - mu) / beta
  log_delta <- .tgumbel_log_delta_window(min, max, mu, beta)
  # log(G(max) - G(min)), the normalisation
  log_norm <- -exp(-(max - mu) / beta) + .log1m_exp(-exp(log_delta))
  result <- -log(beta) + log_s - exp(log_s) - log_norm
  result[is.na(x)] <- NA_real_
  result[!is.na(x) & (x < min | x > max)] <- -Inf
  if (log) {
    result
  } else {
    exp(result)
  }
}

#' @rdname tgumbel
#' @export
ptgumbel <- function(
  q,
  min = 0,
  max = 1,
  mu,
  beta,
  lower.tail = TRUE,
  log.p = FALSE
) {
  .check_tgumbel(min, max, mu, beta)
  inside <- !is.na(q) & q >= min & q <= max
  # Evaluate inside the window only
  x <- pmin(pmax(q, min), max)
  deltas <- .tgumbel_log_deltas(x, min, max, mu, beta)
  log_norm <- .log1m_exp(
    -exp(.tgumbel_log_delta_window(min, max, mu, beta))
  )
  # The CDF is (e^delta_lower - 1) / (e^delta_window - 1) and the upper tail
  # is (1 - e^-delta_upper) / (1 - e^-delta_window), where the window
  # difference is the sum of the lower and upper differences
  log_cdf <- .log1m_exp(-exp(deltas$lower)) - exp(deltas$upper) - log_norm
  log_ccdf <- .log1m_exp(-exp(deltas$upper)) - log_norm
  result <- if (lower.tail) log_cdf else log_ccdf
  below <- !is.na(q) & q < min
  above <- !is.na(q) & q > max
  result[below] <- if (lower.tail) -Inf else 0
  result[above] <- if (lower.tail) 0 else -Inf
  result[is.na(q)] <- NA_real_
  # Guard rounding just outside [-Inf, 0]
  result[inside] <- pmin(result[inside], 0)
  if (log.p) {
    result
  } else {
    exp(result)
  }
}

#' @rdname tgumbel
#' @importFrom stats runif
#' @export
rtgumbel <- function(n, min = 0, max = 1, mu, beta) {
  .check_tgumbel(min, max, mu, beta)
  u <- runif(n)
  log_delta_window <- .tgumbel_log_delta_window(min, max, mu, beta)
  log_norm <- .log1m_exp(-exp(log_delta_window))
  # Invert the upper tail, 1 - F(x) = 1 - u, which is
  # (1 - e^-delta_upper) / (1 - e^-delta_window) with
  # delta_upper = s(x) - s(max). Adding delta_upper to s(max) has no
  # cancellation, and the lower tail is never subtracted from
  delta_upper <- -.log1m_exp(log1p(-u) + log_norm)
  s <- exp(-(max - mu) / beta) + delta_upper
  samples <- mu - beta * log(s)
  pmin(pmax(samples, min), max)
}

attr(dtgumbel, "name") <- "dtgumbel"
attr(ptgumbel, "name") <- "ptgumbel"
