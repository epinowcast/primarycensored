#' Generate random samples from a primary event censored distribution
#'
#' This function generates random samples from a primary event censored
#' distribution. It adjusts the distribution by accounting for the primary
#' event distribution and potential truncation at a maximum delay (D). The
#' function allows for custom primary event distributions and delay
#' distributions.
#'
#' @inheritParams dprimarycensored
#'
#' @param rdist Function to generate random samples from the delay distribution
#' for example [stats::rlnorm()] for lognormal distribution.
#'
#' @param swindow Integer specifying the window size for rounding the delay
#' (default is 1). If `swindow = 0` then no rounding is applied.
#'
#' @param n Number of random samples to generate.
#'
#' @param rprimary Function to generate random samples from the primary
#' distribution (default is [stats::runif()]).
#'
#' @param rprimary_args List of additional arguments to be passed to rprimary.
#'
#' @param oversampling_factor Factor by which to oversample the number of
#' samples to account for truncation (default is 1.2).
#'
#' @param ... Additional arguments to be passed to the distribution function.
#'
#' @return Vector of random samples from the primary event censored
#' distribution censored by the secondary event window.
#'
#' @aliases rpcens
#'
#' @importFrom stats runif
#'
#' @export
#'
#' @details
#' The mathematical formulation for generating random samples from a primary
#' event censored distribution is as follows:
#'
#' 1. Generate primary event times (p) from the specified primary event
#'    distribution (f_p) with parameters phi, defined between 0 and the primary
#'    event window (pwindow):
#'    \deqn{p \sim f_p(\phi), \quad p \in [0, pwindow]}
#'
#' 2. Generate delays (d) from the specified delay distribution (f_d) with
#'    parameters theta:
#'    \deqn{d \sim f_d(\theta)}
#'
#' 3. Calculate the total delays (t) by adding the primary event times and
#'    the delays:
#'    \deqn{t = p + d}
#'
#' 4. Apply truncation (i.e. remove any delays that fall outside the observation
#'    window) to ensure that the delays are within the specified range \[0, D\],
#'    where D is the maximum observable delay:
#'    \deqn{t_{truncated} = \{t \mid 0 \leq t < D\}}
#'
#' 5. Round the truncated delays to the nearest secondary event window
#'    (swindow):
#'    \deqn{t_{valid} = \lfloor \frac{t_{truncated}}{swindow} \rfloor
#'      \times swindow}
#'
#' The function oversamples to account for potential truncation and generates
#' additional samples if needed to reach the desired number of valid samples.
#'
#' @family primarycensored
#'
#' @examples
#' # Example: Lognormal distribution with uniform primary events
#' rprimarycensored(10, rlnorm, meanlog = 0, sdlog = 1)
#'
#' # Example: Lognormal distribution with exponential growth primary events
#' rprimarycensored(
#'   10, rlnorm,
#'   rprimary = rexpgrowth, rprimary_args = list(r = 0.2),
#'   meanlog = 0, sdlog = 1
#' )
rprimarycensored <- function(
    n,
    rdist,
    pwindow = 1,
    swindow = 1,
    D = Inf,
    rprimary = stats::runif,
    rprimary_args = list(),
    oversampling_factor = 1.2,
    ...) {
  # Generate more samples than needed to account for truncation
  n_generate <- ceiling(n * oversampling_factor)

  # Generate primary event times
  p <- do.call(
    rprimary,
    c(list(n_generate), rprimary_args, list(min = 0, max = pwindow))
  )

  # Generate delays from the specified distribution
  delay <- rdist(n_generate, ...)

  # Calculate total delays
  total_delay <- p + delay

  # Apply truncation and select valid samples
  valid_samples <- total_delay[total_delay >= 0 & total_delay < D]

  # If we don't have enough samples, generate more
  while (length(valid_samples) < n) {
    additional_samples <- rprimarycensored(
      n - length(valid_samples),
      rdist,
      pwindow,
      swindow,
      D,
      rprimary,
      rprimary_args,
      ...
    )
    valid_samples <- c(valid_samples, additional_samples)
  }

  # Return exactly n samples
  samples <- valid_samples[1:n]

  # Round to the nearest swindow
  if (swindow > 0) {
    rounded_samples <- floor(samples / swindow) * swindow
  } else {
    rounded_samples <- samples
  }

  return(rounded_samples)
}

#' @rdname rprimarycensored
#' @export
rpcens <- rprimarycensored
