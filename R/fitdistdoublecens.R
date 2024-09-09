#' Fit a distribution to doubly censored data
#'
#' This function wraps the custom approach for fitting distributions to doubly
#' censored data using fitdistrplus and primarycensoreddist.
#'
#' @param censdata A data frame with columns 'left' and 'right' representing
#'   the lower and upper bounds of the censored observations. Unlike
#' [fitdistrplus::fitdistcens()] `NA` is not supported for either the
#' upper or lower bounds.
#'
#' @param distr A character string naming the distribution to be fitted.
#'
#' @inheritParams pprimarycensoreddist
#'
#' @param ... Additional arguments to be passed to [fitdistrplus::fitdist()].
#'
#' @return An object of class "fitdist" as returned by fitdistrplus::fitdist.
#'
#' @export
#'
#' @examples
#' censdata <- data.frame(left = c(1, 2, 3), right = c(2, 3, 4))
#' fit <- fitdistdoublecens(
#'   censdata, "gamma",
#'   start = list(shape = 1, rate = 1)
#' )
fitdistdoublecens <- function(censdata, distr,
                              pwindow = 1, D = Inf,
                              dprimary = stats::dunif,
                              dprimary_args = list(), ...) {
  # Check if fitdistrplus is available
  if (!requireNamespace("fitdistrplus", quietly = TRUE)) {
    stop(
      "Package 'fitdistrplus' is required but not installed for this function."
    )
  }

  if (!all(c("left", "right") %in% names(censdata))) {
    stop("censdata must contain 'left' and 'right' columns")
  }

  # Get the distribution functions
  pdist <- get(paste0("p", distr))

  swindows <- censdata$right - censdata$left

  assign(paste0("dpcens_", distr), function(x, ...) {
    .dpcens(
      x, swindows, pdist, pwindow, D, dprimary, dprimary_args, ...
    )
  })

  assign(paste0("ppcens_", distr), function(q, ...) {
    .ppcens(q, pdist, pwindow, D, dprimary, dprimary_args, ...)
  })


  # Fit the distribution
  fit <- fitdistrplus::fitdist(
    censdata$left,
    distr = paste0("pcens_", distr),
    ...
  )

  return(fit)
}

#' Define a fitdistrplus compatible wrapper around dprimarycensoreddist
#' @inheritParams dprimarycensoreddist
#'
#' @param swindows A numeric vector of secondary window sizes corresponding to
#' each element in x
#' @keywords internal
.dpcens <- function(x, swindows, pdist, pwindow, D, dprimary,
                    dprimary_args, ...) {
  tryCatch(
    {
      vapply(seq_along(x), function(i) {
        dprimarycensoreddist(
          x[i], pdist,
          pwindow = pwindow, swindow = swindows[i], D = D, dprimary = dprimary,
          dprimary_args = dprimary_args, ...
        )
      }, numeric(length(x)))
    },
    error = function(e) {
      rep(NaN, length(x))
    }
  )
}

#' Define a fitdistrplus compatible wrapper around pprimarycensoreddist
#' @inheritParams pprimarycensoreddist
#' @keywords internal
.ppcens <- function(q, pdist, pwindow, D, dprimary, dprimary_args, ...) {
  tryCatch(
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
}
