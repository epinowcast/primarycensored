# Tilt transforms of Weibull and generalised gamma delays, the Stacy family
# with CDF P(k, (x / scale)^shape) for the regularised lower incomplete gamma
# function P. The Weibull has k = 1. The transform is the series
# sum_n (xi scale)^n / n! Gamma(k + n / shape) / Gamma(k) *
#   P(k + n / shape, (t / scale)^shape).

# Terms evaluated together as a matrix, and the most terms summed
.stacy_block <- 32L
.stacy_max_terms <- 400L

#' Stacy family parameters of a pcens object
#'
#' Covers the Weibull and the generalised gamma of `flexsurv::pgengamma.orig()`
#' and, for `Q > 0`, of `flexsurv::pgengamma()`, as in
#' [pcens_cdf.pcens_pgengamma_dunif()].
#'
#' @param object A `pcens` object.
#'
#' @return A list with `shape`, `scale` and `k`, or `NULL` if the parameters
#'   are not positive and finite.
#'
#' @noRd
.stacy_params <- function(object) {
  delay_args <- object$args
  scale <- if (is.null(delay_args$scale)) 1 else delay_args$scale
  p <- if (inherits(object, "pcens_pweibull")) {
    list(shape = delay_args$shape, scale = scale, k = 1)
  } else if (inherits(object, "pcens_pgengamma.orig")) {
    list(shape = delay_args$shape, scale = scale, k = delay_args$k)
  } else {
    mu <- if (is.null(delay_args$mu)) 0 else delay_args$mu
    sigma <- if (is.null(delay_args$sigma)) 1 else delay_args$sigma
    Q <- delay_args$Q
    if (is.numeric(Q) && length(Q) == 1L && is.finite(Q) && Q > 0) {
      list(
        shape = Q / sigma, scale = exp(mu) * Q^(2 * sigma / Q), k = Q^-2
      )
    }
  }
  valid <- vapply(
    p,
    function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0,
    logical(1)
  )
  if (length(p) > 0L && all(valid)) p
}

#' Log of the row sums of the exponentials of a matrix
#'
#' @param m Numeric matrix on the log scale.
#'
#' @return Numeric vector, `-Inf` for a row of `-Inf` or a matrix with no
#'   columns.
#'
#' @noRd
.row_log_sum_exp <- function(m) {
  if (ncol(m) == 0L) {
    return(rep(-Inf, nrow(m)))
  }
  top <- m[cbind(seq_len(nrow(m)), max.col(m, ties.method = "first"))]
  out <- top + log(rowSums(exp(m - top)))
  out[top == -Inf] <- -Inf
  out
}

#' Incomplete gamma series of the tilt transform of a Stacy family delay
#'
#' The terms for negative `xi` alternate in sign and their absolute values sum
#' to the transform at `-xi`. Positive and negative terms are summed apart
#' until the remainder, bounded by a geometric series once `n + 1 > |xi| t`,
#' is below 1e-17 of that sum.
#'
#' @param t Numeric vector of positive points.
#'
#' @param xi Tilt, not zero.
#'
#' @param shape,scale,k Parameters of the Stacy family.
#'
#' @return Numeric vector of log transforms with the attribute `log_loss`, the
#'   log of the sum of the absolute terms over the result. Both are `NaN`
#'   where that factor is more than `.exptilt_max_loss`, `|xi| t` is 400 or
#'   more, or 400 terms do not converge.
#'
#' @noRd
.stacy_series <- function(t, xi, shape, scale, k) {
  x <- (t / scale)^shape
  alternating <- xi < 0
  log_z <- log(abs(xi) * scale)
  log_w <- log(abs(xi) * t)
  w <- exp(log_w)
  log_max_loss <- log(.exptilt_max_loss)
  log_k <- lgamma(k)
  m <- length(t)
  log_pos <- log_neg <- log_first <- rep(-Inf, m)
  out <- rep(NaN, m)
  log_loss <- rep(NaN, m)
  active <- which(w < .stacy_max_terms)
  n_first <- 0L
  while (length(active) > 0L && n_first <= .stacy_max_terms) {
    n <- n_first:min(n_first + .stacy_block - 1L, .stacy_max_terms)
    n_active <- length(active)
    s <- k + n / shape
    log_terms <- matrix(
      stats::pgamma(
        rep(x[active], times = length(n)), rep(s, each = n_active),
        log.p = TRUE
      ) + rep(n * log_z - lgamma(n + 1) + lgamma(s) - log_k, each = n_active),
      nrow = n_active
    )
    if (n_first == 0L) {
      log_first[active] <- log_terms[, 1]
    }
    negative <- alternating & n %% 2L == 1L
    log_pos[active] <- .log_sum_exp(
      log_pos[active], .row_log_sum_exp(log_terms[, !negative, drop = FALSE])
    )
    log_neg[active] <- .log_sum_exp(
      log_neg[active], .row_log_sum_exp(log_terms[, negative, drop = FALSE])
    )
    log_abs <- .log_sum_exp(log_pos[active], log_neg[active])
    n_last <- max(n)
    n_first <- n_last + 1L
    # A CDF that underflows gives no mass. The transform is at most the CDF
    # for negative xi.
    underflow <- log_first[active] == -Inf
    too_lossy <- alternating & log_abs - log_first[active] > log_max_loss
    converged <- !too_lossy & n_last + 1 > w[active] &
      log_terms[, length(n)] + log_w[active] -
        log(pmax(n_last + 1 - w[active], 0)) < log_abs + log(1e-17)
    done <- active[(converged | underflow) & !too_lossy]
    if (length(done) > 0L) {
      result <- .stacy_sum(log_pos[done], log_neg[done], log_max_loss)
      out[done] <- result
      log_loss[done] <- attr(result, "log_loss")
    }
    active <- active[!(converged | underflow) & !too_lossy]
  }
  structure(out, log_loss = log_loss)
}

#' Combine the positive and negative terms of the series
#'
#' @param log_pos,log_neg Log of the sum of the positive and of the negative
#'   terms, `-Inf` if there are none.
#'
#' @param log_max_loss Log of the largest accepted ratio of the sum of the
#'   absolute terms to the result.
#'
#' @return Numeric vector of log sums, `NaN` where the terms cancel too much,
#'   with the attribute `log_loss`.
#'
#' @noRd
.stacy_sum <- function(log_pos, log_neg, log_max_loss) {
  out <- .log_diff_exp(log_pos, log_neg)
  out[log_neg > -Inf & log_pos <= log_neg] <- NaN
  loss <- .log_sum_exp(log_pos, log_neg) - out
  # No mass is a zero sum without loss
  loss[log_pos == -Inf] <- 0
  out[is.na(loss) | loss > log_max_loss] <- NaN
  structure(out, log_loss = loss)
}

#' Log tilt transform of a Stacy family delay
#'
#' The upper transform is the survival function for `xi = 0` and `NaN`
#' otherwise, as the series has no tail form.
#'
#' @inheritParams .stacy_series
#'
#' @param upper If `TRUE` return the transform over `(t, Inf)`.
#'
#' @return Numeric vector of log transforms with the attribute `log_loss`.
#'
#' @noRd
.stacy_log_transform <- function(t, xi, shape, scale, k, upper = FALSE) {
  out <- rep(if (upper) 0 else -Inf, length(t))
  log_loss <- rep(0, length(t))
  if (xi != 0 && upper) {
    return(structure(rep(NaN, length(t)), log_loss = log_loss))
  }
  positive <- which(!is.na(t) & t > 0)
  if (length(positive) > 0L && xi == 0) {
    out[positive] <- stats::pgamma(
      (t[positive] / scale)^shape, k,
      lower.tail = !upper, log.p = TRUE
    )
  } else if (length(positive) > 0L) {
    series <- .stacy_series(t[positive], xi, shape, scale, k)
    out[positive] <- series
    log_loss[positive] <- attr(series, "log_loss")
  }
  structure(out, log_loss = log_loss)
}

#' Log moments of a Stacy family delay about a point
#'
#' As [.gamma_moments()], from the partial moments
#' `scale^j Gamma(k + j / shape) / Gamma(k)` times
#' `P(k + j / shape, (t / scale)^shape)`.
#'
#' @param t Numeric vector of finite points.
#'
#' @inheritParams .stacy_series
#'
#' @return A matrix with columns `G1`, `G2` and `G3`, `-Inf` for `t <= 0`.
#'
#' @noRd
.stacy_moments <- function(t, shape, scale, k) {
  positive <- t > 0
  tp <- pmax(t, 0)
  x <- (tp / scale)^shape
  log_t <- log(tp)
  log_m <- function(j) {
    j * log(scale) + lgamma(k + j / shape) - lgamma(k) +
      stats::pgamma(x, k + j / shape, log.p = TRUE)
  }
  log_g1 <- .log_diff_exp(log_t + log_m(0), log_m(1))
  log_h <- .log_diff_exp(log_t + log_m(1), log_m(2))
  log_g2 <- .log_diff_exp(log_t + log_g1, log_h)
  log_a <- .log_diff_exp(log_t + log_m(2), log_m(3))
  log_b <- .log_diff_exp(log_t + log_h, log_a)
  log_g3 <- .log_diff_exp(log_t + log_g2, log_b)
  cbind(
    G1 = ifelse(positive, log_g1, -Inf),
    G2 = ifelse(positive, log_g2, -Inf),
    G3 = ifelse(positive, log_g3, -Inf)
  )
}

.stacy_tilt_lower <- function(object) 0

.stacy_tilt_available <- function(object, xi) !is.null(.stacy_params(object))

.stacy_tilt_transform <- function(object, t, xi, upper = FALSE) {
  p <- .stacy_params(object)
  .stacy_log_transform(t, xi, p$shape, p$scale, p$k, upper)
}

.stacy_tilt_moments <- function(object, t) {
  p <- .stacy_params(object)
  .stacy_moments(t, p$shape, p$scale, p$k)
}

#' @exportS3Method
.pcens_tilt_lower.pcens_pweibull <- .stacy_tilt_lower

#' @exportS3Method
.pcens_tilt_lower.pcens_pgengamma.orig <- .stacy_tilt_lower

#' @exportS3Method
.pcens_tilt_lower.pcens_pgengamma <- .stacy_tilt_lower

#' @exportS3Method
.pcens_tilt_available.pcens_pweibull <- .stacy_tilt_available

#' @exportS3Method
.pcens_tilt_available.pcens_pgengamma.orig <- .stacy_tilt_available

#' @exportS3Method
.pcens_tilt_available.pcens_pgengamma <- .stacy_tilt_available

#' @exportS3Method
.pcens_tilt_transform.pcens_pweibull <- .stacy_tilt_transform

#' @exportS3Method
.pcens_tilt_transform.pcens_pgengamma.orig <- .stacy_tilt_transform

#' @exportS3Method
.pcens_tilt_transform.pcens_pgengamma <- .stacy_tilt_transform

#' @exportS3Method
.pcens_tilt_moments.pcens_pweibull <- .stacy_tilt_moments

#' @exportS3Method
.pcens_tilt_moments.pcens_pgengamma.orig <- .stacy_tilt_moments

#' @exportS3Method
.pcens_tilt_moments.pcens_pgengamma <- .stacy_tilt_moments
