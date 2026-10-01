#' Primary event censored CDF from per-endpoint terms
#'
#' Evaluates the analytical CDF of a `pcens` object when its class has a
#' terms specification, see [.pcens_terms_spec()], and otherwise falls back
#' to [pcens_cdf.default()].
#'
#' @inheritParams pcens_cdf
#'
#' @return Vector of computed primary event censored CDFs.
#'
#' @keywords internal
.pcens_cdf_analytic <- function(object, q, pwindow, use_numeric) {
  spec <- if (!isTRUE(use_numeric)) .pcens_terms_spec(object)
  if (is.null(spec)) {
    return(pcens_cdf.default(object, q, pwindow, use_numeric))
  }
  .pcens_cdf_shared(spec, q, pwindow)
}

#' Evaluate a CDF from terms shared between endpoints
#'
#' The CDF at `q` combines terms at `q` and at `max(q - pwindow, 0)`. Terms
#' at `q` are reused for lower endpoints that equal another point. This is
#' done when it saves at least a quarter of the evaluations. Otherwise, and
#' for a vector `pwindow`, the terms are evaluated directly at both
#' endpoints.
#'
#' @inheritParams pcens_cdf
#'
#' @param spec A terms specification from [.pcens_terms_spec()].
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
  shared <- n > 1L && length(pwindow) == 1L
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

#' Terms specification for an analytical primary event censored CDF
#'
#' This is the registry of analytical solutions built from terms that each
#' depend on one endpoint. A solution is added by returning its
#' specification for the class of the `pcens` object. Objects with no
#' specification return `NULL` and use the numerical or class specific
#' methods.
#'
#' @param object A `pcens` object as created by [new_pcens()].
#'
#' @return `NULL`, or a list with `terms`, a function of time returning a
#'   list of the terms, and `combine`, a function of the terms at the upper
#'   and lower endpoints and `pwindow` returning the CDF.
#'
#' @keywords internal
.pcens_terms_spec <- function(object) {
  dist_args <- object$args
  switch(class(object)[1],
    pcens_pgamma_dunif = .uniform_terms_gamma(dist_args),
    pcens_plnorm_dunif = .uniform_terms_lnorm(dist_args),
    pcens_pweibull_dunif = .uniform_terms_weibull(dist_args),
    pcens_pgengamma.orig_dunif = .uniform_terms_gengamma(
      dist_args$shape, dist_args$scale, dist_args$k
    ),
    pcens_pgengamma_dunif = .uniform_terms_prentice(dist_args),
    NULL
  )
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
#' @return A terms specification, see [.pcens_terms_spec()].
#'
#' @keywords internal
.uniform_terms <- function(terms, mean) {
  list(
    terms = terms,
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
    shape * scale
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
    exp(mu + 0.5 * sigma^2)
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
    scale
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
    scale * exp(lgamma(k_shift) - lgamma(k))
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
