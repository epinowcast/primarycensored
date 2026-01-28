#' Generate random samples from a primary event censored distribution
#'
#' This function generates random samples from a primary event censored
#' distribution. It adjusts the distribution by accounting for the primary
#' event distribution and potential truncation at a maximum delay (D) and
#' minimum delay (L). The function allows for custom primary event
#' distributions and delay distributions.
#'
#' @inheritParams dprimarycensored
#'
#' @param rdist Function to generate random samples from the delay distribution
#'  for example [stats::rlnorm()] for lognormal distribution.
#'
#' @param swindow Integer specifying the window size for rounding the delay
#'  (default is 1). If `swindow = 0` then no rounding is applied.
#'
#' @param n Number of random samples to generate.
#'
#' @param rprimary Function to generate random samples from the primary
#'  distribution (default is [stats::runif()]).
#'
#' @param rprimary_args List of additional arguments to be passed to rprimary.
#'
#' @param oversampling_factor Factor by which to oversample the number of
#'  samples to account for truncation (default is 1.2).
#'
#' @param ... Additional arguments to be passed to the distribution function.
#'
#' @return Vector of random samples from the primary event censored
#'  distribution censored by the secondary event window.
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
#' 4. Apply upper truncation to remove delays >= D:
#'    \deqn{t_{upper} = \{t \mid t < D\}}
#'
#' 5. Round the delays to the nearest secondary event window (swindow):
#'    \deqn{t_{rounded} = \lfloor \frac{t_{upper}}{swindow} \rfloor
#'      \times swindow}
#'
#' 6. Apply lower truncation on the rounded values to ensure observed delays
#'    are >= L:
#'    \deqn{t_{valid} = \{t_{rounded} \mid t_{rounded} \geq L\}}
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
#'
#' # Example: Left-truncated distribution (e.g., for generation intervals)
#' rprimarycensored(10, rlnorm, L = 1, D = 10, meanlog = 0, sdlog = 1)
rprimarycensored <- function(
    n,
    rdist,
    pwindow = 1,
    swindow = 1,
    L = 0,
    D = Inf,
    rprimary = stats::runif,
    rprimary_args = list(),
    oversampling_factor = 1.2,
    ...) {
  .check_truncation_bounds(L, D)

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

  # Apply upper truncation (continuous delays must be < D)
  valid_continuous <- total_delay[total_delay < D]

  # Round to the nearest swindow
  if (swindow > 0) {
    rounded <- floor(valid_continuous / swindow) * swindow
  } else {
    rounded <- valid_continuous
  }

  # Apply lower truncation on observed values (must be >= L)
  valid_samples <- rounded[rounded >= L]

  # If we don't have enough samples, generate more
  while (length(valid_samples) < n) {
    additional_samples <- rprimarycensored(
      n - length(valid_samples),
      rdist,
      pwindow,
      swindow,
      L,
      D,
      rprimary,
      rprimary_args,
      ...
    )
    valid_samples <- c(valid_samples, additional_samples)
  }

  # Return exactly n samples
  return(valid_samples[1:n])
}

#' @rdname rprimarycensored
#' @export
rpcens <- rprimarycensored
