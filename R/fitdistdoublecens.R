#' Fit a distribution to doubly censored data
#'
#' This function wraps the custom approach for fitting distributions to doubly
#' censored data using fitdistrplus and primarycensoreddist.
#'
#' @param censdata A data frame with columns 'left' and 'right' representing
#'   the lower and upper bounds of the censored observations.
#' @param distr A character string naming the distribution to be fitted.
#' @param start A named list giving the initial values for the parameters of
#'   the named distribution.
#' @param pwindow Numeric, primary event window (default: 1).
#' @param D Numeric, maximum delay for truncation (default: Inf).
#' @param dprimary Function for the primary event distribution
#'   (default: stats::dunif).
#' @param dprimary_args List of additional arguments for dprimary.
#' @param ... Additional arguments to be passed to fitdist.
#'
#' @return An object of class "fitdist" as returned by fitdistrplus::fitdist.
#'
#' @importFrom stats get
#'
#' @export
#'
#' @examples
#' # Example usage:
#' # censdata <- data.frame(left = c(1, 2, 3), right = c(2, 3, 4))
#' # fit <- fitdistdoublecens(censdata, "gamma", start = list(shape = 1, rate = 1))
fitdistdoublecens <- function(censdata, distr, start,
                              pwindow = 1, D = Inf,
                              dprimary = stats::dunif,
                              dprimary_args = list(), ...) {
  # Check if fitdistrplus is available
  if (!requireNamespace("fitdistrplus", quietly = TRUE)) {
    stop(
      "Package 'fitdistrplus' is required but not installed for this function."
    )
  }

  # Validate input
  if (!all(c("left", "right") %in% names(censdata))) {
    stop("censdata must contain 'left' and 'right' columns")
  }

  # Get the distribution functions
  pdist <- get(paste0("p", distr))

  # Calculate unique secondary intervals and their frequencies
  swindows <- censdata$right - censdata$left
  unique_swindows <- unique(swindows)
  swindow_freq <- table(swindows)

  # Fit the distribution
  fit <- fitdistrplus::fitdist(
    censdata$left,
    distr = paste0("pcens_", distr),
    start = start,
    custom.family = list(
      dname = paste0("dpcens_", distr),
      pname = paste0("ppcens_", distr)
    ),
    densfun = function(x, ...) dpcens(x, pdist, pwindow, D, dprimary, dprimary_args, unique_swindows, swindow_freq, ...),
    probfun = function(q, ...) ppcens_(q, pdist, pwindow, D, dprimary, dprimary_args, ...),
    ...
  )

  return(fit)
}

#' Define custom density function
#' @keywords internal
dpcens <- function(x, pdist, pwindow, D, dprimary, dprimary_args, unique_swindows, swindow_freq, ...) {
  result <- tryCatch(
    {
      densities <- vapply(unique_swindows, function(sw) {
        dprimarycensoreddist(
          x, pdist,
          pwindow = pwindow, swindow = sw,
          D = D, dprimary = dprimary, dprimary_args = dprimary_args, ...
        )
      }, numeric(length(x)))

      # Scale densities by frequency and sum
      rowSums(densities * rep(swindow_freq, each = length(x)))
    },
    error = function(e) {
      rep(NaN, length(x))
    }
  )
  return(result)
}

#' Define custom distribution function
#' @keywords internal
ppcens_ <- function(q, pdist, pwindow, D, dprimary, dprimary_args, ...) {
  result <- tryCatch(
    {
      pprimarycensoreddist(
        q, pdist,
        pwindow = pwindow,
        D = D, dprimary = dprimary, dprimary_args = dprimary_args, ...
      )
    },
    error = function(e) {
      rep(NaN, length(q))
    }
  )
  return(result)
}
