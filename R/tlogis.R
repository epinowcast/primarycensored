#' Truncated logistic distribution functions
#'
#' Density, distribution function, and random generation for the logistic
#' distribution truncated to the interval \[min, max\]. As a primary event
#' distribution it describes a cumulative exposure probability that increases
#' smoothly and then saturates, for example an environmental source that
#' becomes established.
#'
#' @param x,q Vector of quantiles.
#'
#' @param n Number of observations. If `length(n) > 1`, the length is taken to
#'  be the number required.
#'
#' @param min Minimum value of the distribution range. Default is 0.
#'
#' @param max Maximum value of the distribution range. Default is 1.
#'
#' @param location Location of the logistic distribution before truncation,
#'  its midpoint, on the same scale as `min` and `max`. Default is 0.
#'
#' @param scale Scale of the logistic distribution before truncation. A
#'  larger scale gives a slower transition. It must be positive. Default is 1.
#'
#' @param log,log.p Logical; if TRUE, probabilities p are given as log(p).
#'
#' @param lower.tail Logical; if TRUE (default), probabilities are P\[X <= x\],
#'  otherwise, P\[X > x\].
#'
#' @return `dtlogis` gives the density, `ptlogis` gives the distribution
#' function, and `rtlogis` generates random deviates.
#'
#' The length of the result is determined by `n` for `rtlogis`, and is the
#' maximum of the lengths of the numerical arguments for the other functions.
#'
#' @details
#' With the logistic CDF
#' \eqn{L(z) = 1 / (1 + \exp\{-(z - \text{location}) / \text{scale}\})} and
#' the mass \eqn{D_L = L(max) - L(min)} of the interval, the probability
#' density function on \[min, max\] is
#'
#' \deqn{f(x) = \frac{L'(x)}{D_L}}
#'
#' and the cumulative distribution function is
#'
#' \deqn{F(x) = \frac{L(x) - L(min)}{D_L}.}
#'
#' As a primary event distribution the window is \[0, `pwindow`\], so `min`
#' is 0 and `max` is the window.
#'
#' The differences of the logistic CDF are computed on the log scale from
#' the lower tails where the interval is below `location`, and from the upper
#' tails where it is above. This keeps \eqn{D_L} and the CDF accurate when the
#' interval is far from `location` and \eqn{L(max) - L(min)} would otherwise
#' round to zero.
#'
#' For random number generation, we use the inverse transform method with the
#' same tails. A uniform draw \eqn{u} gives
#' \eqn{x = L^{-1}(L(min) + u D_L)}, evaluated on the log scale.
#' With a large `scale` relative to the interval the distribution approaches
#' the uniform distribution.
#'
#' @family primaryeventdistributions
#'
#' @examples
#' x <- seq(0, 1, by = 0.1)
#' dens <- dtlogis(x, location = 0.5, scale = 0.2)
#' cumprobs <- ptlogis(x, location = 0.5, scale = 0.2)
#' samples <- rtlogis(100, location = 0.5, scale = 0.2)
#'
#' @name tlogis
NULL

#' Check the arguments of the truncated logistic distribution functions
#'
#' @inheritParams tlogis
#'
#' @return `NULL` invisibly. Called for its error.
#'
#' @noRd
.tlogis_check <- function(min, max, location, scale) {
  if (anyNA(min) || anyNA(max) || any(min >= max)) {
    stop("min must be less than max.", call. = FALSE)
  }
  if (anyNA(location) || !all(is.finite(location))) {
    stop("location must be finite.", call. = FALSE)
  }
  if (anyNA(scale) || !all(is.finite(scale)) || any(scale <= 0)) {
    stop("scale must be positive and finite.", call. = FALSE)
  }
  invisible(NULL)
}

#' Log of a difference of logistic CDFs
#'
#' Evaluates \eqn{\log(L(hi) - L(lo))} for `lo <= hi` from the lower tails
#' where `lo` is below `location` and from the upper tails otherwise, so the
#' difference is not lost to rounding far from `location`.
#'
#' @param lo,hi Numeric vectors of interval ends with `lo <= hi`.
#'
#' @inheritParams tlogis
#'
#' @return Numeric vector of log differences, `-Inf` where they are zero.
#'
#' @noRd
.tlogis_log_diff <- function(lo, hi, location, scale) {
  n <- max(length(lo), length(hi), length(location), length(scale))
  lo <- rep_len(lo, n)
  hi <- rep_len(hi, n)
  location <- rep_len(location, n)
  scale <- rep_len(scale, n)
  above <- lo >= location
  out <- numeric(n)
  below <- which(!above)
  above <- which(above)
  if (length(below) > 0L) {
    out[below] <- .log_diff_exp(
      stats::plogis(
        hi[below], location[below], scale[below],
        log.p = TRUE
      ),
      stats::plogis(
        lo[below], location[below], scale[below],
        log.p = TRUE
      )
    )
  }
  if (length(above) > 0L) {
    out[above] <- .log_diff_exp(
      stats::plogis(
        lo[above], location[above], scale[above],
        lower.tail = FALSE, log.p = TRUE
      ),
      stats::plogis(
        hi[above], location[above], scale[above],
        lower.tail = FALSE, log.p = TRUE
      )
    )
  }
  out
}

#' @rdname tlogis
#' @export
dtlogis <- function(x, min = 0, max = 1, location = 0, scale = 1,
                    log = FALSE) {
  .tlogis_check(min, max, location, scale)
  result <- stats::dlogis(x, location, scale, log = TRUE) -
    .tlogis_log_diff(min, max, location, scale)
  result[x < min | x > max] <- -Inf
  if (log) {
    result
  } else {
    exp(result)
  }
}

#' @rdname tlogis
#' @export
ptlogis <- function(q, min = 0, max = 1, location = 0, scale = 1,
                    lower.tail = TRUE, log.p = FALSE) {
  .tlogis_check(min, max, location, scale)
  # The interval above q (or below it) is empty outside the window
  inside <- pmin(pmax(q, min), max)
  log_mass <- .tlogis_log_diff(min, max, location, scale)
  result <- if (lower.tail) {
    .tlogis_log_diff(min, inside, location, scale) - log_mass
  } else {
    .tlogis_log_diff(inside, max, location, scale) - log_mass
  }
  if (log.p) {
    result
  } else {
    exp(result)
  }
}

#' @rdname tlogis
#' @importFrom stats runif
#' @export
rtlogis <- function(n, min = 0, max = 1, location = 0, scale = 1) {
  .tlogis_check(min, max, location, scale)
  n <- if (length(n) > 1L) length(n) else n
  u <- runif(n)
  min <- rep_len(min, n)
  max <- rep_len(max, n)
  location <- rep_len(location, n)
  scale <- rep_len(scale, n)
  # log(1 - u) and log(u) mix the tails of the ends of the interval. The
  # quantile is taken from the lower tails where the interval is below the
  # location and from the upper tails otherwise.
  above <- min >= location
  samples <- numeric(n)
  for (upper in c(FALSE, TRUE)) {
    idx <- which(above == upper)
    if (length(idx) == 0L) {
      next
    }
    log_tail_min <- stats::plogis(
      min[idx], location[idx], scale[idx],
      lower.tail = !upper, log.p = TRUE
    )
    log_tail_max <- stats::plogis(
      max[idx], location[idx], scale[idx],
      lower.tail = !upper, log.p = TRUE
    )
    log_tail <- .log_sum_exp(
      log1p(-u[idx]) + log_tail_min, log(u[idx]) + log_tail_max
    )
    samples[idx] <- stats::qlogis(
      log_tail, location[idx], scale[idx],
      lower.tail = !upper, log.p = TRUE
    )
  }
  pmin(pmax(samples, min), max)
}

attr(dtlogis, "name") <- "dtlogis"
attr(ptlogis, "name") <- "ptlogis"

# Multiples of the scale from the centre of the density at which the numerical
# method breaks the integral. The mass beyond 30 scales is below 1e-13.
.tlogis_break_multiples <- c(-30, -12, -6, -3, -1.5, 0, 1.5, 3, 6, 12, 30)

#' Location and scale of the truncated logistic primary of a pcens object
#'
#' @param object A `pcens` object with a truncated logistic primary.
#'
#' @return A list with `location` and `scale`, defaulting as in [dtlogis()].
#'
#' @noRd
.tlogis_primary_args <- function(object) {
  m <- object$primary_args$location
  s <- object$primary_args$scale
  if (is.null(m)) {
    m <- 0
  }
  if (is.null(s)) {
    s <- 1
  }
  if (!is.numeric(m) || length(m) != 1L || !is.finite(m)) {
    stop(
      "location must be a single finite number for the truncated logistic ",
      "primary distribution",
      call. = FALSE
    )
  }
  if (!is.numeric(s) || length(s) != 1L || !is.finite(s) || s <= 0) {
    stop(
      "scale must be a single positive finite number for the truncated ",
      "logistic primary distribution",
      call. = FALSE
    )
  }
  list(location = m, scale = s)
}
