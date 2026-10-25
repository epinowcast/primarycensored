#' Method for a normal delay with a truncated Gumbel primary
#'
#' Analytical primary event censored CDF for a normal delay distribution with
#' a truncated Gumbel primary event window, the [dtgumbel()] primary
#' distribution. It honours `use_numeric`, and uses a numerical method where
#' the series below is not accurate.
#'
#' @inheritParams pcens_cdf
#'
#' @details
#' Write \eqn{G(z) = \exp[-\exp\{-(z - \mu) / \beta\}]} and
#' \eqn{D = G(w) - G(0)} for a window of width \eqn{w}, and let \eqn{F} be
#' the delay CDF and \eqn{T_f(\xi; \tau)} the transform of [pcens_cdf_exptilt].
#' With \eqn{c(q) = \exp\{-(q - \mu) / \beta\}} the reversed window
#' \eqn{G(q - u) = \exp\{-c(q) e^{u / \beta}\}} is a power series in
#' \eqn{e^{u / \beta}}, which gives
#' \deqn{
#' F_G(q) = F(q - w) + \frac{1}{D} \Big[
#'   (1 - G(0)) \{F(q) - F(q - w)\} +
#'   \sum_{n \ge 1} \frac{\{-c(q)\}^n}{n!} \Delta_n(q) \Big],
#' }
#' with \eqn{\Delta_n(q) = T_f(n / \beta; q) - T_f(n / \beta; q - w)}.
#' The transforms at each endpoint are evaluated once at every tilt and reused
#' when several `q` share an endpoint.
#' Each \eqn{\Delta_n} is taken between the lower or the upper tail
#' transforms, whichever loses less precision, as for [pcens_cdf_exptilt].
#'
#' The terms alternate and are bounded by \eqn{s_0^n / n!} times the delay mass
#' in the window, where \eqn{s_0 = e^{\mu / \beta}}.
#' They are summed on the log scale with the odd and even terms separate.
#' The series is truncated where the bound is below 1e-18.
#'
#' **Fallback to the numerical method.** The relative error of the series is
#' estimated for each `q`.
#' Where it is above 1e-8 the numerical method is used for that `q`.
#' The numerical method is used for every `q` where \eqn{\mu / \beta} is
#' above `log(15)` or where the window is narrow relative to \eqn{\beta}.
#'
#' **Numerical method.** The window density is a spike of width about
#' \eqn{\beta} when \eqn{\mu / \beta} is large, which the integration of
#' [pcens_cdf.default()] can miss.
#' So [pcens_cdf.default()] integrates this primary in a variable in which
#' the integrand is smooth, over the range that holds the mass of the window.
#' This applies to any delay with this primary.
#'
#' @inherit pcens_cdf return
#'
#' @family pcens
#'
#' @name pcens_cdf_gumbel
#'
#' @examples
#' pnorm_obj <- new_pcens(
#'   pdist = pnorm, dprimary = dtgumbel,
#'   primary_args = list(mu = 0.5, beta = 0.5), mean = 3, sd = 2
#' )
#' pcens_cdf(pnorm_obj, q = c(-1, 3, 8), pwindow = 2)
NULL

#' @rdname pcens_cdf_gumbel
#' @export
pcens_cdf.pcens_pnorm_dtgumbel <- function(
  object,
  q,
  pwindow,
  use_numeric = FALSE
) {
  primary <- .gumbel_primary_args(object)
  if (length(pwindow) != 1L || !is.finite(pwindow) || pwindow <= 0) {
    return(pcens_cdf.default(object, q, pwindow, use_numeric))
  }
  .gumbel_cdf(
    object, q, pwindow, primary$mu, primary$beta,
    series = !isTRUE(use_numeric) &&
      .gumbel_available(object, primary$mu, primary$beta)
  )
}

# The series is not used above s0 = exp(mu / beta) of 15, as the terms reach
# exp(s0) and their number grows with s0
.gumbel_max_s0 <- 15

# Largest estimated relative error for which the series is used, about the
# accuracy of the numerical method
.gumbel_tol <- 1e-8

# Largest rounding estimate for which the series is tried
.gumbel_screen_tol <- 1e-7

# Relative tolerance of each piece of the numerical integral
.gumbel_num_tol <- 1e-10

# Largest u = s(z) - s(w) integrated, beyond which exp(-u) underflows
.gumbel_u_max <- 745

# Scales above the peak at which the density is cut, exp(-46) of its peak
.gumbel_z_cut <- 46

#' Number of terms of the truncated Gumbel series
#'
#' @param log_s0 Log of \eqn{s_0 = e^{\mu / \beta}}, that is `mu / beta`.
#'
#' @return The smallest integer `N`, at least 2, for which
#'   \eqn{s_0^n / n!} is below `1e-18 * max(s0, 1)` after `N` terms.
#'
#' @keywords internal
.gumbel_n_terms <- function(log_s0) {
  n <- seq_len(200L)
  log_bound <- n * log_s0 - lgamma(n + 1) - max(log_s0, 0)
  below <- which(log_bound < log(1e-18))
  max(2L, if (length(below) > 0L) below[[1L]] else length(n))
}

#' Location and scale of a truncated Gumbel primary
#'
#' @inheritParams pcens_cdf
#'
#' @return List with `mu` and `beta`, after checking that they are valid
#'   single numbers.
#'
#' @keywords internal
.gumbel_primary_args <- function(object) {
  primary <- object$primary_args[c("mu", "beta")]
  if (is.null(primary$mu) || is.null(primary$beta)) {
    stop(
      "mu and beta parameters are required for the truncated Gumbel ",
      "primary distribution",
      call. = FALSE
    )
  }
  .check_tgumbel(0, 1, primary$mu, primary$beta)
  if (length(primary$mu) != 1L || length(primary$beta) != 1L) {
    stop(
      "mu and beta must be single numbers for the truncated Gumbel ",
      "primary distribution",
      call. = FALSE
    )
  }
  primary
}

#' Whether the truncated Gumbel series can be used
#'
#' The choice depends on the parameters only.
#' The rounding error is about the machine precision times the size of the
#' log transform at the largest tilt times the terms \eqn{e^{s_0}}, so the
#' series is skipped where that is above `.gumbel_screen_tol`.
#'
#' @inheritParams pcens_cdf
#' @inheritParams tgumbel
#'
#' @return `TRUE` if the series can be used.
#'
#' @keywords internal
.gumbel_available <- function(object, mu, beta) {
  log_s0 <- mu / beta
  if (log_s0 > log(.gumbel_max_s0)) {
    return(FALSE)
  }
  xi <- .gumbel_n_terms(log_s0) / beta
  p <- .norm_mean_sd(object)
  size <- xi * abs(p$mean) + 0.5 * (xi * p$sd)^2
  log(.Machine$double.eps) + log1p(size) + exp(log_s0) <=
    log(.gumbel_screen_tol)
}

#' Truncated Gumbel primary censored CDF for a finite positive window
#'
#' Uses the series where `series` is `TRUE` and its estimated error is below
#' `.gumbel_tol`, and the numerical method otherwise.
#'
#' @inheritParams pcens_cdf
#' @inheritParams tgumbel
#'
#' @param series Logical; whether to try the series.
#'
#' @return Vector of CDFs.
#'
#' @keywords internal
.gumbel_cdf <- function(object, q, pwindow, mu, beta, series = FALSE) {
  result <- rep(NA_real_, length(q))
  result[!is.na(q) & q == Inf] <- 1
  result[!is.na(q) & q == -Inf] <- 0
  finite <- which(is.finite(q))
  if (length(finite) == 0L) {
    return(result)
  }
  q <- q[finite]
  numeric_needed <- rep(TRUE, length(q))
  if (series) {
    fit <- .gumbel_lcdf(object, q, pwindow, mu, beta)
    result[finite] <- pmin(1, exp(fit$log_cdf))
    numeric_needed <- is.na(fit$error) | fit$error > .gumbel_tol
  }
  if (any(numeric_needed)) {
    result[finite[numeric_needed]] <- .gumbel_numeric(
      object, q[numeric_needed], pwindow, mu, beta
    )
  }
  result
}

#' Numerical truncated Gumbel primary censored CDF
#'
#' The integral is taken in \eqn{u = s(z) - s(w)}, with
#' \eqn{s(z) = \exp\{(\mu - z) / \beta\}}, for \eqn{\mu} at or above the
#' window width \eqn{w}, and in \eqn{z} otherwise, so a narrow spike of the
#' density is resolved.
#'
#' @inheritParams .gumbel_cdf
#'
#' @return Vector of CDFs at finite `q`.
#'
#' @keywords internal
.gumbel_numeric <- function(object, q, pwindow, mu, beta) {
  log_sw <- (mu - pwindow) / beta
  cdf <- function(x) do.call(object$pdist, c(list(q = x), object$args))
  integrate_pieces <- function(integrand, breaks) {
    breaks <- sort(unique(breaks))
    sum(vapply(
      seq_len(length(breaks) - 1L),
      function(i) {
        fit <- stats::integrate(
          integrand, breaks[i], breaks[i + 1L],
          rel.tol = .gumbel_num_tol, abs.tol = 0, subdivisions = 1000L,
          stop.on.error = FALSE
        )
        if (fit$message != "OK") {
          warning(
            "Truncated Gumbel integration: ", fit$message, call. = FALSE
          )
        }
        fit$value
      },
      numeric(1)
    ))
  }
  if (log_sw >= 0) {
    log_delta <- log_sw + .log_diff_exp(pwindow / beta, 0)
    top <- min(exp(log_delta), .gumbel_u_max)
    scale <- -expm1(-exp(log_delta))
    breaks_u <- c(0, top, c(1, 5, 20, 60)[c(1, 5, 20, 60) < top])
    return(vapply(
      q,
      function(d) {
        # Keeps a small d - pwindow when u / s_w is tiny
        integrand <- function(u) {
          x <- (d - pwindow) + beta * log1p(exp(log(u) - log_sw))
          cdf(pmin(x, d)) * exp(-u)
        }
        kink <- if (d > 0 && d < pwindow) {
          log_sw + .log_diff_exp((pwindow - d) / beta, 0)
        } else {
          Inf
        }
        breaks <- if (kink < log(top)) c(breaks_u, exp(kink)) else breaks_u
        min(1, max(0, integrate_pieces(integrand, breaks) / scale))
      },
      numeric(1)
    ))
  }
  z_low <- max(0, mu - beta * log(.gumbel_u_max + exp(log_sw)))
  z_high <- min(pwindow, max(mu, 0) + beta * .gumbel_z_cut)
  breaks_z <- c(z_low, z_high, mu + beta * c(-3, -1, 0, 1, 3, 8, 20))
  breaks_z <- breaks_z[breaks_z >= z_low & breaks_z <= z_high]
  vapply(
    q,
    function(d) {
      integrand <- function(z) {
        cdf(d - z) * dtgumbel(z, 0, pwindow, mu, beta)
      }
      breaks <- c(breaks_z, if (d > z_low && d < z_high) d)
      min(1, max(0, integrate_pieces(integrand, breaks)))
    },
    numeric(1)
  )
}

#' Truncated Gumbel log CDF from the series
#'
#' @inheritParams .gumbel_cdf
#'
#' @return List with the log CDFs `log_cdf` and their estimated relative
#'   errors `error`.
#'
#' @keywords internal
.gumbel_lcdf <- function(object, q, pwindow, mu, beta) {
  n_terms <- .gumbel_n_terms(mu / beta)
  pos <- .exptilt_layout(q, pwindow, .pcens_tilt_lower(object))
  # Log transforms with a row per endpoint and a column per tilt n / beta,
  # n = 0, ..., N
  n_end <- length(pos$endpoints)
  ends <- rep(pos$endpoints, n_terms + 1L)
  xi <- rep(seq.int(0L, n_terms) / beta, each = n_end)
  lower_terms <- matrix(.pcens_tilt_transform(object, ends, xi), n_end)
  upper_terms <- matrix(
    .pcens_tilt_transform(object, ends, xi, upper = TRUE), n_end
  )
  at <- function(terms, index) terms[index, , drop = FALSE]
  log_delta <- matrix(
    .exptilt_tail_diff(
      c(at(lower_terms, pos$q)), c(at(lower_terms, pos$y)),
      c(at(upper_terms, pos$q)), c(at(upper_terms, pos$y))
    ),
    nrow = length(q)
  )
  finite_abs <- function(x) {
    x <- abs(x)
    x[!is.finite(x)] <- 0
    x
  }
  scale_t <- pmax(
    finite_abs(at(lower_terms, pos$q)), finite_abs(at(upper_terms, pos$q)),
    finite_abs(at(lower_terms, pos$y)), finite_abs(at(upper_terms, pos$y))
  )
  .gumbel_combine(
    log_delta, lower_terms[pos$y, 1L], scale_t, q, pwindow, mu, beta
  )
}

#' Sum the truncated Gumbel series
#'
#' @inheritParams .gumbel_cdf
#'
#' @param log_delta Matrix of the log of \eqn{\Delta_n(q)}, with a row per
#'   `q` and a column per tilt, \eqn{n = 0, \dots, N}.
#'
#' @param log_f_y Log delay CDF at `q - pwindow`.
#'
#' @param scale_t Matrix, shaped as `log_delta`, of the size of the log
#'   transforms, for the rounding estimate.
#'
#' @return List with the log CDFs `log_cdf` and their estimated relative
#'   errors `error`.
#'
#' @keywords internal
.gumbel_combine <- function(log_delta, log_f_y, scale_t, q, pwindow, mu,
                            beta) {
  n_terms <- ncol(log_delta) - 1L
  n <- seq_len(n_terms)
  log_s0 <- mu / beta
  log_c <- -(q - mu) / beta
  log_a <- sweep(log_delta[, n + 1L, drop = FALSE], 2L, lgamma(n + 1), "-") +
    outer(log_c, n)
  # Positive terms are the even n and (1 - G(0)) Delta_0, negative the odd
  is_odd <- n %% 2L == 1L
  log_pos <- .log_sum_exp_rows(cbind(
    .log1m_exp_neg_exp(log_s0) + log_delta[, 1L],
    log_a[, !is_odd, drop = FALSE]
  ))
  log_neg <- .log_sum_exp_rows(log_a[, is_odd, drop = FALSE])
  log_bracket <- .log_diff_exp(log_pos, log_neg)
  # log D = -s_w + log(1 - exp(-s_w expm1(w / beta)))
  log_s_w <- -(pwindow - mu) / beta
  log_d <- -exp(log_s_w) +
    .log1m_exp_neg_exp(log_s_w + .log_diff_exp(pwindow / beta, 0))
  log_cdf <- .log_sum_exp(log_f_y, log_bracket - log_d)
  # No delay mass in the window, so all terms vanish together
  no_mass <- is.infinite(log_delta[, 1L]) & log_delta[, 1L] < 0
  log_cdf[no_mass] <- log_f_y[no_mass]
  error <- .gumbel_error_bound(
    log_pos, log_neg, log_bracket, log_delta[, 1L], log_c, log_s0, n_terms,
    scale_t, n
  )
  error[no_mass] <- 0
  list(log_cdf = log_cdf, error = error)
}

#' Estimated relative error of the truncated Gumbel series
#'
#' The machine precision times the sum of the absolute terms and the size of
#' the logs added, plus the bound \eqn{s_0^{N + 1} / (N + 1)!} on the
#' truncation times the delay mass.
#'
#' @inheritParams .gumbel_combine
#' @inheritParams .gumbel_n_terms
#'
#' @param log_pos,log_neg Logs of the sums of the positive and negative terms.
#'
#' @param log_bracket Log of the bracket, `log_pos` minus `log_neg`.
#'
#' @param log_mass Log of the delay mass in the window.
#'
#' @param log_c Log of \eqn{c(q) = \exp\{-(q - \mu) / \beta\}}.
#'
#' @param n_terms Number of terms, \eqn{N}.
#'
#' @param n The term indices, `seq_len(n_terms)`.
#'
#' @return Vector of estimated relative errors.
#'
#' @keywords internal
.gumbel_error_bound <- function(log_pos, log_neg, log_bracket, log_mass,
                                log_c, log_s0, n_terms, scale_t, n) {
  magnitude <- 1 + .row_max(
    abs(outer(log_c, n)) + scale_t[, n + 1L, drop = FALSE]
  )
  log_rounding <- log(.Machine$double.eps) + log(magnitude) +
    .log_sum_exp(log_pos, log_neg) - log_bracket
  log_trunc <- log_mass + (n_terms + 1) * log_s0 - lgamma(n_terms + 2) -
    log1p(-min(exp(log_s0) / (n_terms + 2), 0.5)) - log_bracket
  error <- exp(log_rounding) + exp(log_trunc)
  error[!is.finite(log_bracket)] <- Inf
  error
}

#' Row maxima and row log-sum-exp of a matrix
#'
#' @param x Numeric matrix.
#'
#' @return Vector with one value per row.
#'
#' @keywords internal
.row_max <- function(x) {
  x[cbind(seq_len(nrow(x)), max.col(x, ties.method = "first"))]
}

#' @rdname dot-row_max
.log_sum_exp_rows <- function(x) {
  top <- .row_max(x)
  shifted <- x - ifelse(is.finite(top), top, 0)
  top + log(rowSums(exp(shifted)))
}
