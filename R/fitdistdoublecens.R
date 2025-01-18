#' Fit a distribution to doubly censored data
#'
#' This function wraps the custom approach for fitting distributions to doubly
#' censored data using fitdistrplus and primarycensored.
#'
#' @details
#' This function temporarily assigns and then removes functions from the global
#' environment in order to work with fitdistr. Users should be aware of this
#' behaviour, especially if they have existing functions with the same names in
#' their global environment.
#'
#' @param censdata A data frame with columns 'left' and 'right' representing
#' the lower and upper bounds of the censored observations. Unlike
#' [fitdistrplus::fitdistcens()] `NA` is not supported for either the
#' upper or lower bounds.
#'
#' @param distr A character string naming the distribution to be fitted.
#'
#' @inheritParams pprimarycensored
#'
#' @param ... Additional arguments to be passed to [fitdistrplus::fitdist()].
#'
#' @param truncation_check_multiplier Numeric multiplier to use for checking
#'   if the truncation time D is appropriate relative to the maximum delay.
#'   Set to NULL to skip the check. Default is 2.
#'
#' @return An object of class "fitdist" as returned by fitdistrplus::fitdist.
#'
#' @export
#' @family modelhelpers
#' @examplesIf requireNamespace("fitdistrplus", quietly = TRUE)
#' # Example with normal distribution
#' set.seed(123)
#' n <- 1000
#' true_mean <- 5
#' true_sd <- 2
#' pwindow <- 2
#' swindow <- 2
#' D <- 10
#' samples <- rprimarycensored(
#'   n, rnorm,
#'   mean = true_mean, sd = true_sd,
#'   pwindow = pwindow, swindow = swindow, D = D
#' )
#'
#' delay_data <- data.frame(
#'   left = samples,
#'   right = samples + swindow
#' )
#'
#' fit_norm <- fitdistdoublecens(
#'   delay_data,
#'   distr = "norm",
#'   start = list(mean = 0, sd = 1),
#'   D = D, pwindow = pwindow
#' )
#'
#' summary(fit_norm)
fitdistdoublecens <- function(censdata, distr,
                              pwindow = 1, D = Inf,
                              dprimary = stats::dunif,
                              dprimary_name = deprecated(),
                              dprimary_args = list(),
                              truncation_check_multiplier = 2,
                              ...) {
  .name_deprecation(pdist_name, dprimary_name)
  # Check if fitdistrplus is available
  if (!requireNamespace("fitdistrplus", quietly = TRUE)) {
    stop(
      "Package 'fitdistrplus' is required but not installed for this function."
    )
  }

  if (!requireNamespace("withr", quietly = TRUE)) {
    stop(
      "Package 'withr' is required but not installed for this function."
    )
  }

  if (!all(c("left", "right") %in% names(censdata))) {
    stop("censdata must contain 'left' and 'right' columns")
  }

  if (!is.null(truncation_check_multiplier)) {
    check_truncation(
      delays = censdata$left,
      D = D,
      multiplier = truncation_check_multiplier
    )
  }

  # Get the distribution functions
  pdist_name <- paste0("p", distr)
  pdist <- get(pdist_name)
  swindows <- censdata$right - censdata$left

  # Create the function definition with named arguments for dpcens
  dpcens_dist <- function() {
    args <- as.list(environment())
    do.call(.dpcens, c(
      args,
      list(
        swindows = swindows,
        pdist = pdist,
        pwindow = pwindow,
        D = D,
        dprimary = dprimary,
        dprimary_args = dprimary_args
      )
    ))
  }
  formals(dpcens_dist) <- formals(get(paste0("d", distr)))

  # Create the function definition with named arguments for ppcens
  ppcens_dist <- function() {
    args <- as.list(environment())
    do.call(.ppcens, c(
      args,
      list(
        pdist = pdist,
        pwindow = pwindow,
        D = D,
        dprimary = dprimary,
        dprimary_args = dprimary_args
      )
    ))
  }
  formals(ppcens_dist) <- formals(pdist)

  # Fit the distribution
  fit <- withr::with_environment(environment(), fitdistrplus::fitdist(
    censdata$left,
    distr = "pcens_dist",
    ...
  ))
  return(fit)
}

#' Define a fitdistrplus compatible wrapper around dprimarycensored
#' @inheritParams dprimarycensored
#'
#' @param swindows A numeric vector of secondary window sizes corresponding to
#' each element in x
#' @keywords internal
.dpcens <- function(x, swindows, pdist, pwindow, D, dprimary,
                    dprimary_args, pdist_name, dprimary_name, ...) {
  .name_deprecation(pdist_name, dprimary_name)
  tryCatch(
    {
      if (length(unique(swindows)) == 1) {
        dprimarycensored(
          x, pdist,
          pwindow = pwindow, swindow = swindows[1], D = D, dprimary = dprimary,
          dprimary_args = dprimary_args, ...
        )
      } else {
        # Group x and swindows by unique swindow values
        unique_swindows <- unique(swindows)
        result <- numeric(length(x))

        for (sw in unique_swindows) {
          mask <- swindows == sw
          result[mask] <- dprimarycensored(
            x[mask], pdist,
            pwindow = pwindow, swindow = sw, D = D,
            dprimary = dprimary, dprimary_args = dprimary_args,
            ...
          )
        }

        result
      }
    },
    error = function(e) {
      rep(NaN, length(x))
    }
  )
}

#' Define a fitdistrplus compatible wrapper around pprimarycensored
#' @inheritParams pprimarycensored
#' @keywords internal
.ppcens <- function(q, pdist, pwindow, D, dprimary, dprimary_args,
                    ...) {
  tryCatch(
    {
      pprimarycensored(
        q, pdist,
        pwindow = pwindow,
        D = D, dprimary = dprimary, dprimary_args = dprimary_args,
        ...
      )
    },
    error = function(e) {
      rep(NaN, length(q))
    }
  )
}
