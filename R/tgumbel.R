#' Truncated Gumbel distribution functions
#'
#' Density, distribution function, and random generation for a Gumbel
#' distribution truncated to `[min, max]`. It is a primary event distribution
#' for exposure that builds up or decays unevenly within the window.
#'
#' @inheritParams expgrowth
#'
#' @param mu Location of the Gumbel distribution before truncation, on the
#'  same scale as `x`.
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
#' The density peaks at \eqn{x = \mu} and falls more slowly above it.
#' If \eqn{\mu} is far below `min` it decays with rate \eqn{1 / \beta}, as for
#' [dexpgrowth()] with \eqn{r = -1 / \beta}.
#' If \eqn{\mu} is far above `max` it is a narrow spike at `max`.
#'
#' Differences of \eqn{s(x) = \exp\{-(x - \mu) / \beta\}} are evaluated on the
#' log scale, so the tails stay precise when \eqn{\mu} is far from the window.
#' Random numbers are drawn by inverting the upper tail.
#'
#' The Stan equivalents are `tgumbel_lpdf()`, `tgumbel_lcdf()` and
#' `tgumbel_rng()`, with the primary distribution identifier 4 and
#' `primary_params = [mu, beta]`.
#' Analytical primary event censored CDFs are described in
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

# Checks mu, beta and the bounds, called for its error
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

# log(1 - exp(-exp(x))), which is x below -37 where exp(x) would underflow
.log1m_exp_neg_exp <- function(x) {
  out <- x
  large <- which(x >= -37)
  out[large] <- .log1m_exp(-exp(x[large]))
  out
}

# Log of s(min) - s(max), with s(x) = exp(-(x - mu) / beta)
.tgumbel_log_delta_window <- function(min, max, mu, beta) {
  -(max - mu) / beta + .log_diff_exp((max - min) / beta, 0)
}

# Logs of s(min) - s(x) and s(x) - s(max)
.tgumbel_log_deltas <- function(x, min, max, mu, beta) {
  list(
    lower = -(x - mu) / beta + .log_diff_exp((x - min) / beta, 0),
    upper = -(max - mu) / beta + .log_diff_exp((max - x) / beta, 0)
  )
}

#' @rdname tgumbel
#' @export
dtgumbel <- function(x, min = 0, max = 1, mu, beta, log = FALSE) {
  .check_tgumbel(min, max, mu, beta)
  inside_x <- pmin(pmax(x, min), max)
  log_delta <- .tgumbel_log_delta_window(min, max, mu, beta)
  # The exponent s(x) - s(max) is a difference, so it does not cancel where
  # s(max) is large
  log_upper <- .tgumbel_log_deltas(inside_x, min, max, mu, beta)$upper
  result <- -log(beta) - (inside_x - mu) / beta - exp(log_upper) -
    .log1m_exp_neg_exp(log_delta)
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
  x <- pmin(pmax(q, min), max)
  deltas <- .tgumbel_log_deltas(x, min, max, mu, beta)
  log_norm <- .log1m_exp_neg_exp(
    .tgumbel_log_delta_window(min, max, mu, beta)
  )
  log_cdf <- .log1m_exp_neg_exp(deltas$lower) - exp(deltas$upper) - log_norm
  log_ccdf <- .log1m_exp_neg_exp(deltas$upper) - log_norm
  result <- if (lower.tail) log_cdf else log_ccdf
  below <- !is.na(q) & q < min
  above <- !is.na(q) & q > max
  result[below] <- if (lower.tail) -Inf else 0
  result[above] <- if (lower.tail) 0 else -Inf
  result[is.na(q)] <- NA_real_
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
  log_norm <- .log1m_exp_neg_exp(.tgumbel_log_delta_window(min, max, mu, beta))
  # Invert the upper tail, x = max - beta log(1 + (s(x) - s(max)) / s(max))
  log_delta_upper <- log(-.log1m_exp(log1p(-u) + log_norm))
  log_ratio <- log_delta_upper + (max - mu) / beta
  samples <- max - beta * (pmax(log_ratio, 0) + log1p(exp(-abs(log_ratio))))
  pmin(pmax(samples, min), max)
}

attr(dtgumbel, "name") <- "dtgumbel"
attr(ptgumbel, "name") <- "ptgumbel"
