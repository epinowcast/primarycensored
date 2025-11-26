#' Fit a distribution to doubly censored data
#'
#' This function wraps the custom approach for fitting distributions to doubly
#' censored data using fitdistrplus and primarycensored. It handles primary
#' censoring (when the primary event time is not known exactly), secondary
#' censoring (when the secondary event time is interval-censored), and right
#' truncation (when events are only observed up to a maximum delay).
#'
#' @details
#' ## How distribution functions are resolved
#'
#' The `distr` parameter specifies the base name of the distribution. The
#' function automatically looks up the corresponding density (`d`) and
#' cumulative distribution (`p`) functions by prepending these prefixes to the
#' distribution name. For example:
#' - `distr = "gamma"` uses `dgamma()` and `pgamma()`
#' - `distr = "lnorm"` uses `dlnorm()` and `plnorm()`
#' - `distr = "weibull"` uses `dweibull()` and `pweibull()`
#'
#' Any distribution available in base R or loaded packages can be used, as long
#' as the corresponding `d<distr>` and `p<distr>` functions exist and follow
#' standard R distribution function conventions (first argument is `x` for
#' density, `q` for CDF).
#'
#' ## What this function does internally
#'
#' This function creates custom density and CDF functions that account for
#' primary censoring, secondary censoring, and truncation using
#' [dprimarycensored()] and [pprimarycensored()]. These custom functions are
#' then passed to [fitdistrplus::fitdist()] for maximum likelihood estimation.
#'
#' The function handles varying observation windows across observations,
#' making it suitable for real-world data where truncation times or censoring
#' windows may differ between observations.
#'
#' ## Environment handling
#'
#' This function temporarily assigns and then removes functions from the global
#' environment in order to work with fitdistr. Users should be aware of this
#' behaviour, especially if they have existing functions with the same names in
#' their global environment.
#'
#' @param censdata A data frame with columns 'left' and 'right' representing
#'  the lower and upper bounds of the censored observations. Unlike
#'  [fitdistrplus::fitdistcens()] `NA` is not supported for either the
#'  upper or lower bounds.
#'
#' @param distr A character string naming the distribution to be fitted. This
#'  should be the base name of a distribution with corresponding `d` (density)
#'  and `p` (cumulative distribution) functions available. For example, use
#'  `"gamma"` (which will use `dgamma` and `pgamma`), `"lnorm"` (for `dlnorm`
#'  and `plnorm`), `"weibull"`, `"norm"`, etc. Custom distributions can also
#'  be used as long as the corresponding `d<distr>()` and `p<distr>()`
#'  functions are defined and loaded.
#'
#' @param left Column name for lower bound of observed values (default: "left").
#'
#' @param right Column name for upper bound of observed values (default:
#'  "right").
#'
#' @param pwindow Column name for primary window (default: "pwindow").
#'
#' @param D Column name for maximum delay (truncation point). If finite, the
#'  distribution is truncated at D. If set to Inf, no truncation is applied.
#'  (default: "D").
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
#'   right = samples + swindow,
#'   pwindow = rep(pwindow, n),
#'   D = rep(D, n)
#' )
#'
#' fit_norm <- fitdistdoublecens(
#'   delay_data,
#'   distr = "norm",
#'   start = list(mean = 0, sd = 1)
#' )
#'
#' summary(fit_norm)
fitdistdoublecens <- function(
  censdata,
  distr,
  left = "left",
  right = "right",
  pwindow = "pwindow",
  D = "D",
  dprimary = stats::dunif,
  dprimary_name = lifecycle::deprecated(),
  dprimary_args = list(),
  truncation_check_multiplier = 2,
  ...
) {
  nms <- .name_deprecation(lifecycle::deprecated(), dprimary_name)
  if (!is.null(nms$dprimary)) {
    dprimary <- add_name_attribute(dprimary, nms$dprimary)
  }

  # Check if fitdistrplus is available
  if (!requireNamespace("fitdistrplus", quietly = TRUE)) {
    stop(
      "Package 'fitdistrplus' is required but not installed for this function.",
      call. = FALSE
    )
  }

  if (!requireNamespace("withr", quietly = TRUE)) {
    stop(
      "Package 'withr' is required but not installed for this function.",
      call. = FALSE
    )
  }

  # Deprecation handling for pwindow and D
  if (is.numeric(pwindow)) {
    lifecycle::deprecate_warn(
      "1.1.0",
      "fitdistdoublecens(pwindow)",
      details = "Use pwindow column in censdata instead."
    )
    censdata[["pwindow"]] <- pwindow
    pwindow <- "pwindow"
  }
  if (is.numeric(D)) {
    lifecycle::deprecate_warn(
      "1.1.0",
      "fitdistdoublecens(D)",
      details = "Use D column in censdata instead."
    )
    censdata[["D"]] <- D
    D <- "D"
  }

  required_cols <- c(left, right, pwindow, D)
  missing_cols <- setdiff(required_cols, names(censdata))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns: ",
      toString(missing_cols),
      call. = FALSE
    )
  }

  if (!is.null(truncation_check_multiplier)) {
    unique_D <- unique(censdata[[D]])
    for (d in unique_D) {
      delays_subset <- censdata[[left]][censdata[[D]] == d]
      check_truncation(
        delays = delays_subset,
        D = d,
        multiplier = truncation_check_multiplier
      )
    }
  }

  # Get the distribution functions
  pdist_name <- paste0("p", distr)
  pdist <- add_name_attribute(get(pdist_name), pdist_name)
  params <- data.frame(
    swindow = censdata[[right]] - censdata[[left]],
    pwindow = censdata[[pwindow]],
    D = censdata[[D]]
  )

  # Create the function definition with named arguments for dpcens
  dpcens_dist <- function() {
    env_args <- as.list(environment())
    do.call(
      .dpcens,
      c(
        env_args,
        list(
          params = params,
          pdist = pdist,
          dprimary = dprimary,
          dprimary_args = dprimary_args
        )
      )
    )
  }
  formals(dpcens_dist) <- formals(get(paste0("d", distr)))

  # Create the function definition with named arguments for ppcens
  ppcens_dist <- function() {
    env_args <- as.list(environment())
    do.call(
      .ppcens,
      c(
        env_args,
        list(
          params = params,
          pdist = pdist,
          dprimary = dprimary,
          dprimary_args = dprimary_args
        )
      )
    )
  }
  formals(ppcens_dist) <- formals(pdist)

  delays <- censdata[[left]]
  # Create a clean environment with only the necessary objects
  fit_env <- new.env(parent = emptyenv())

  # Copy only the required objects to the clean environment
  fit_env$delays <- delays
  fit_env$ppcens_dist <- ppcens_dist
  fit_env$dpcens_dist <- dpcens_dist

  # Perform the fitting in the clean environment
  fit <- withr::with_environment(
    fit_env,
    fitdistrplus::fitdist(
      delays,
      distr = "pcens_dist",
      ...
    )
  )
  return(fit)
}

#' Define a fitdistrplus compatible wrapper around dprimarycensored
#' @inheritParams dprimarycensored
#'
#' @param params A data frame with columns 'swindow', 'pwindow', and 'D'
#' corresponding to the secondary window sizes, primary window sizes, and
#' truncation times for each element in x.
#' @keywords internal
.dpcens <- function(
  x,
  params,
  pdist,
  dprimary,
  dprimary_args,
  ...
) {
  tryCatch(
    {
      unique_params <- unique(params)
      # Check if all parameters are constant
      if (nrow(unique_params) == 1) {
        dprimarycensored(
          x,
          pdist,
          pwindow = unique_params$pwindow[1],
          swindow = unique_params$swindow[1],
          D = unique_params$D[1],
          dprimary = dprimary,
          dprimary_args = dprimary_args,
          ...
        )
      } else {
        # Group by unique combinations of parameters
        result <- numeric(length(x))

        for (i in seq_len(nrow(unique_params))) {
          sw <- unique_params$swindow[i]
          pw <- unique_params$pwindow[i]
          Ds <- unique_params$D[i] # nolint
          mask <- params$swindow == sw &
            params$pwindow == pw &
            params$D == Ds

          result[mask] <- dprimarycensored(
            x[mask],
            pdist,
            pwindow = pw,
            swindow = sw,
            D = Ds,
            dprimary = dprimary,
            dprimary_args = dprimary_args,
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
.ppcens <- function(q, params, pdist, dprimary, dprimary_args, ...) {
  tryCatch(
    {
      # Vectorize the CDF calculation
      mapply(
        function(q_i, pw, D_i) {
          pprimarycensored(
            q_i,
            pdist,
            pwindow = pw,
            D = D_i,
            dprimary = dprimary,
            dprimary_args = dprimary_args,
            ...
          )
        },
        q,
        params$pwindow,
        params$D,
        SIMPLIFY = TRUE
      )
    },
    error = function(e) {
      rep(NaN, length(q))
    }
  )
}
