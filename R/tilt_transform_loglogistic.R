# Largest |xi| t for the series, and smallest shape for the series. The
# series needs powers of 1 / shape that overflow when tiny.
.loglogistic_tilt_limit <- 10
.loglogistic_min_shape <- 0.2

#' Number of terms of the transform series
#'
#' @noRd
.loglogistic_terms_needed <- function(z) {
  n <- 1L
  weight <- z
  while (n < 400L && (n <= z || weight > 1e-19)) {
    n <- n + 1L
    weight <- weight * z / n
  }
  n
}

#' Sum of the transform series
#'
#' @noRd
.loglogistic_psi <- function(xi, t, shape, log_A) {
  z <- abs(xi) * t
  n <- seq_len(.loglogistic_terms_needed(max(z)))
  ratios <- .loglogistic_ratio(n / shape, log_A)
  weight <- exp(
    outer(log(z), n) - matrix(lgamma(n + 1), length(t), length(n), byrow = TRUE)
  )
  if (xi < 0) {
    weight <- weight * matrix((-1)^n, length(t), length(n), byrow = TRUE)
  }
  1 + rowSums(weight * ratios)
}

#' Log of the moments G_k(t) of the delay about t
#'
#' @noRd
.loglogistic_log_moments <- function(t, shape, scale, orders) {
  out <- matrix(-Inf, length(t), length(orders))
  colnames(out) <- paste0("G", orders)
  positive <- t > 0
  if (!any(positive)) {
    return(out)
  }
  log_t <- log(t[positive])
  log_ratio <- shape * (log_t - log(scale))
  log_cdf <- -.log1p_exp(-log_ratio)
  r <- .loglogistic_ratio(seq_len(max(orders)) / shape, log_ratio)
  out[positive, 1L] <- log_t + log_cdf + log1p(-r[, 1L])
  if (2L %in% orders) {
    out[positive, 2L] <- 2 * log_t + log_cdf +
      log(1 - 2 * r[, 1L] + r[, 2L])
  }
  if (3L %in% orders) {
    out[positive, 3L] <- 3 * log_t + log_cdf +
      log(1 - 3 * r[, 1L] + 3 * r[, 2L] - r[, 3L])
  }
  out
}

#' @exportS3Method
.pcens_tilt_lower.pcens_pllogis <- function(object) {
  0
}

#' @exportS3Method
.pcens_tilt_available.pcens_pllogis <- function(object, xi) {
  .loglogistic_shape_scale(object)$shape >= .loglogistic_floor_shape
}

#' @exportS3Method
.pcens_tilt_transform.pcens_pllogis <- function(object, t, xi, upper = FALSE) {
  # The upper transform is only available for xi = 0, and the lower transform
  # is NaN where the series is not accurate
  p <- .loglogistic_shape_scale(object)
  out <- rep(if (upper) if (xi == 0) 0 else NaN else -Inf, length(t))
  positive <- t > 0
  if (!any(positive) || (upper && xi != 0)) {
    return(out)
  }
  log_ratio <- p$shape * (log(t[positive]) - log(p$scale))
  if (xi == 0) {
    out[positive] <- -.log1p_exp(if (upper) log_ratio else -log_ratio)
    return(out)
  }
  inside <- positive & abs(xi) * t <= .loglogistic_tilt_limit &
    p$shape >= .loglogistic_min_shape
  out[positive] <- NaN
  if (any(inside)) {
    log_ratio <- p$shape * (log(t[inside]) - log(p$scale))
    out[inside] <- -.log1p_exp(-log_ratio) +
      log(.loglogistic_psi(xi, t[inside], p$shape, log_ratio))
  }
  out
}

#' @exportS3Method
.pcens_tilt_ill_conditioned.pcens_pllogis <- function(
  object, q, pwindow, rho, log_cdf, small_window
) {
  p <- .loglogistic_shape_scale(object)
  if (small_window) {
    bounds <- .loglogistic_window_error(object, q, pwindow)
    return(.log_error_exceeds(
      bounds$log_error, log_cdf, bounds$log_density
    ))
  }
  ill <- rep(TRUE, length(q))
  inside <- abs(rho) * q <= .loglogistic_tilt_limit &
    p$shape >= .loglogistic_min_shape
  if (any(inside)) {
    ill[inside] <- .loglogistic_direct_ill(
      object, q[inside], pwindow, rho, log_cdf[inside]
    )
  }
  ill
}

#' Check if the direct form is ill conditioned at each point
#'
#' An error in a transform, bounded by 4e-16 F(t) exp(z) / sqrt(1 + z) for
#' z = |rho| t, is multiplied by \eqn{e^{\rho q} / |e^{\rho w} - 1|} in the
#' combination.
#'
#' @noRd
.loglogistic_direct_ill <- function(object, q, pwindow, rho, log_cdf) {
  p <- .loglogistic_shape_scale(object)
  log_error <- function(t) {
    out <- rep(-Inf, length(t))
    positive <- t > 0
    z <- abs(rho) * t[positive]
    out[positive] <- log(.loglogistic_error_constant) + z - 0.5 * log1p(z) -
      .log1p_exp(-p$shape * (log(t[positive]) - log(p$scale)))
    out
  }
  log_den <- if (rho > 0) {
    .log_diff_exp(rho * pwindow, 0)
  } else {
    .log1m_exp(rho * pwindow)
  }
  .log_error_exceeds(
    rho * q + .log_sum_exp(log_error(q), log_error(q - pwindow)) - log_den,
    log_cdf
  )
}

#' @exportS3Method
.pcens_tilt_numeric.pcens_pllogis <- function(object, q, pwindow) {
  .loglogistic_numeric_cdf(object, q, pwindow)
}

#' @exportS3Method
.pcens_tilt_moments.pcens_pllogis <- function(object, t) {
  p <- .loglogistic_shape_scale(object)
  .loglogistic_log_moments(t, p$shape, p$scale, 1:3)
}

#' @rdname pcens_cdf_loglogistic
#' @export
pcens_cdf.pcens_pllogis_dexpgrowth <- function(
  object,
  q,
  pwindow,
  use_numeric = FALSE
) {
  .pcens_cdf_exptilt(object, q, pwindow, use_numeric)
}

#' @rdname pcens_cdf_loglogistic
#' @export
pcens_cdf.pcens_pllogis_dexpgrowth <- function(
  object,
  q,
  pwindow,
  use_numeric = FALSE
) {
  .pcens_cdf_exptilt(object, q, pwindow, use_numeric)
}
