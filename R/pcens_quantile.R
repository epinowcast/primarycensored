#' Compute primary event censored quantiles
#'
#' This function inverts the primary event censored CDF to compute quantiles.
#' It uses numerical optimisation via optim to find the value q such that
#' [pcens_cdf()] is close to the specified probability. Currently, only the
#' default numerical inversion method is implemented. Future analytical
#' solutions may be added.
#'
#' @inheritParams pprimarycensored
#'
#' @param object A `primarycensored` object as created by [new_pcens()].
#'
#' @param p A vector of probabilities at which to compute the quantiles.
#'
#' @param use_numeric Logical; if TRUE forces the use of numeric inversion even
#' if an analytical solution is available (not yet implemented).
#'
#' @return Vector of primary event censored quantiles.
#'
#' @family pcens
#'
#' @export
pcens_quantile <- function(
    object,
    p,
    pwindow,
    D = Inf,
    L = 0,
    use_numeric = FALSE,
    ...) {
  UseMethod("pcens_quantile")
}

#' Default method for computing primary event censored quantiles
#'
#' This method inverts the primary event censored CDF using numerical
#' optimisation via optim. For each probability value, it searches for the
#' delay such that the CDF computed by [pcens_cdf()] approximates the target
#' probability.
#'
#' @param init Initial guess for the delay. By default, 5.
#'
#' @param tol Numeric tolerance for the convergence criterion in the
#'   optimisation routine.
#'
#' @param max_iter Integer specifying the maximum number of iterations allowed
#'   during optimisation.
#'
#' @param ... Additional arguments passed to underlying functions.
#'
#' @inheritParams pcens_quantile
#'
#' @details
#' The quantile is computed by minimising the squared difference between the
#' computed CDF and the target probability.
#'
#' @family pcens
#'
#' @return A numeric vector containing the computed primary event censored
#'   quantiles.
#'
#' @export
#' @examples
#' # Create a primarycensored object with gamma delay and uniform primary
#' pcens_obj <- new_pcens(
#'   pdist = pgamma,
#'   dprimary = dunif,
#'   dprimary_args = list(min = 0, max = 1),
#'   shape = 3,
#'   scale = 2
#' )
#'
#' # Compute quantile for a single probability
#' pcens_quantile(pcens_obj, p = 0.8, pwindow = 1)
#'
#' # Compute quantiles for multiple probabilities
#' pcens_quantile(pcens_obj, p = c(0.25, 0.5, 0.75), pwindow = 1)
#'
#' # Compute quantiles for multiple probabilities with truncation
#' pcens_quantile(pcens_obj, p = c(0.25, 0.5, 0.75), pwindow = 1, D = 10)
#'
#' # Compute quantiles with left truncation
#' pcens_quantile(pcens_obj, p = c(0.25, 0.5, 0.75), pwindow = 1, D = 10, L = 1)
pcens_quantile.default <- function(
    object,
    p,
    pwindow,
    D = Inf,
    L = 0,
    use_numeric = FALSE,
    init = 5,
    tol = 1e-8,
    max_iter = 10000,
    ...) {
  sapply(p, function(prob) {
    # Handle boundary cases.
    if (prob <= 0) {
      return(L)
    }
    if (prob >= 1) {
      return(NA_real_)
    }

    # Objective function: squared difference between the CDF value and prob.
    objective <- function(q) {
      cdf_val <- pcens_cdf(object, q, pwindow, use_numeric)
      if (!is.infinite(D) || L > 0) {
        cdf_val <- .normalise_cdf(cdf_val, q, D, L, object, pwindow)
      }
      (cdf_val - prob)^2
    }

    # Lower bound is L (minimum truncation point)
    lower_bound <- L

    opt_result <- stats::optim(
      par = max(init, L),
      fn = objective,
      method = "L-BFGS-B",
      lower = lower_bound,
      control = list(fnscale = 1, maxit = max_iter, factr = tol)
    )

    opt_result$par
  })
}
