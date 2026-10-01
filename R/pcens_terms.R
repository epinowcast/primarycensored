#' Evaluate a CDF from terms shared between endpoints
#'
#' The CDF at `q` combines terms at `q` and at `max(q - pwindow, 0)`. Terms
#' at `q` are reused for lower endpoints that equal another point. This is
#' done when it saves at least a quarter of the evaluations. Otherwise, and
#' for a vector `pwindow`, the terms are evaluated directly at both
#' endpoints. Terms are only shared when the delay parameters are scalar.
#'
#' @inheritParams pcens_cdf
#'
#' @param spec A terms specification from [.uniform_terms()].
#'
#' @return Vector of CDF values in \[0, 1\]. Missing values in `q` are an
#'   error.
#'
#' @keywords internal
.pcens_cdf_shared <- function(spec, q, pwindow) {
  if (anyNA(q)) {
    stop("q must not contain missing values.", call. = FALSE)
  }
  lower <- pmax.int(q - pwindow, 0)
  n <- length(q)
  shared <- n > 1L && length(pwindow) == 1L && isTRUE(spec$scalar)
  if (shared) {
    # Probing a few points skips the full match when nothing coincides
    probe <- if (n <= 6L) seq_len(n) else round(seq.int(1, n, length.out = 6))
    hits <- 0L
    for (i in probe) {
      hits <- hits + !is.na(match(lower[i], q))
    }
    shared <- 2L * hits >= length(probe)
  }
  if (shared) {
    idx <- match(lower, q)
    unmatched <- is.na(idx)
    extra <- unique(lower[unmatched])
    shared <- 2 * length(extra) <= n
  }
  if (shared) {
    idx[unmatched] <- n + match(lower[unmatched], extra)
    terms_d <- spec$terms(q)
    terms_q <- Map(
      function(d, e) c(d, e)[idx], terms_d, spec$terms(extra)
    )
  } else {
    terms_d <- spec$terms(q)
    terms_q <- spec$terms(lower)
  }
  res <- pmin.int(1, pmax.int(0, spec$combine(terms_d, terms_q, pwindow)))
  res[q <= 0] <- 0
  res
}

#' Uniform primary terms specification
#'
#' The CDF is
#' `(d F(d) - q F(q) - E (G(d) - G(q))) / pwindow` with `E` the delay mean
#' and `G` the partial expectation CDF, so `terms` returns `t F(t)` and
#' `G(t)` and is zero for `t <= 0`.
#'
#' @param terms Function of time returning a list of two vectors.
#'
#' @param mean Mean of the delay distribution.
#'
#' @param params List of the delay distribution parameters.
#'
#' @return A list with `terms`, `combine`, a function of the terms at the
#'   upper and lower endpoints and `pwindow` returning the CDF, and
#'   `scalar`, whether all of `params` have length one.
#'
#' @keywords internal
.uniform_terms <- function(terms, mean, params) {
  list(
    terms = terms,
    scalar = all(lengths(params) == 1L),
    combine = function(terms_d, terms_q, pwindow) {
      (terms_d[[1]] - terms_q[[1]] - mean * (terms_d[[2]] - terms_q[[2]])) /
        pwindow
    }
  )
}

#' Uniform primary terms for each delay distribution
#'
#' @param dist_args Delay distribution parameters.
#'
#' @inherit .uniform_terms return
#'
#' @keywords internal
.uniform_terms_gamma <- function(dist_args) {
  shape <- dist_args$shape
  scale <- dist_args$scale
  if (is.null(scale) && !is.null(dist_args$rate)) {
    scale <- 1 / dist_args$rate
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
  .uniform_terms(
    function(t) {
      t <- pmax.int(t, 0)
      list(
        t * pgamma(t, shape = shape, scale = scale),
        pgamma(t, shape = shape + 1, scale = scale)
      )
    },
    shape * scale,
    list(shape, scale)
  )
}

#' @rdname dot-uniform_terms_gamma
#' @keywords internal
.uniform_terms_lnorm <- function(dist_args) {
  mu <- dist_args$meanlog
  sigma <- dist_args$sdlog
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
  .uniform_terms(
    function(t) {
      t <- pmax.int(t, 0)
      list(
        t * stats::plnorm(t, meanlog = mu, sdlog = sigma),
        stats::plnorm(t, meanlog = mu + sigma^2, sdlog = sigma)
      )
    },
    exp(mu + 0.5 * sigma^2),
    list(mu, sigma)
  )
}

#' @rdname dot-uniform_terms_gamma
#' @keywords internal
.uniform_terms_weibull <- function(dist_args) {
  shape <- dist_args$shape
  scale <- dist_args$scale
  if (is.null(shape)) {
    stop("shape parameter is required for Weibull distribution", call. = FALSE)
  }
  if (is.null(scale)) {
    stop("scale parameter is required for Weibull distribution", call. = FALSE)
  }
  inv_scale <- 1 / scale
  a <- 1 + 1 / shape
  lgamma_a <- lgamma(a)
  # Regularised form avoids overflow in the unregularised series
  g <- function(t) {
    exp(pgamma((t * inv_scale)^shape, shape = a, log.p = TRUE) + lgamma_a)
  }
  .uniform_terms(
    function(t) {
      t <- pmax.int(t, 0)
      list(t * stats::pweibull(t, shape = shape, scale = scale), g(t))
    },
    scale,
    list(shape, scale)
  )
}

#' @rdname dot-uniform_terms_gamma
#' @param shape,scale,k Generalised gamma parameters in the Stacy
#'   parameterisation of `flexsurv::pgengamma.orig()`.
#' @keywords internal
.uniform_terms_gengamma <- function(shape, scale, k) {
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
  k_shift <- k + 1 / shape
  .uniform_terms(
    function(t) {
      t <- pmax.int(t, 0)
      x <- (t / scale)^shape
      list(t * pgamma(x, shape = k), pgamma(x, shape = k_shift))
    },
    scale * exp(lgamma(k_shift) - lgamma(k)),
    list(shape, scale, k)
  )
}

#' @rdname dot-uniform_terms_gamma
#' @keywords internal
.uniform_terms_prentice <- function(dist_args) {
  mu <- if (is.null(dist_args$mu)) 0 else dist_args$mu
  sigma <- if (is.null(dist_args$sigma)) 1 else dist_args$sigma
  Q <- dist_args$Q
  if (is.null(Q) || Q <= 0) {
    return(NULL)
  }
  .uniform_terms_gengamma(
    shape = Q / sigma, scale = exp(mu) * Q^(2 * sigma / Q), k = Q^-2
  )
}
