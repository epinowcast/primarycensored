% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitdistdoublecens.R
\name{fitdistdoublecens}
\alias{fitdistdoublecens}
\title{Fit a distribution to doubly censored data}
\usage{
fitdistdoublecens(
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
)
}
\arguments{
\item{censdata}{A data frame with columns 'left' and 'right' representing
the lower and upper bounds of the censored observations. Unlike
\code{\link[fitdistrplus:fitdistcens]{fitdistrplus::fitdistcens()}} \code{NA} is not supported for either the
upper or lower bounds.}

\item{distr}{A character string naming the distribution to be fitted.}

\item{left}{Column name for lower bound of observed values (default: "left").}

\item{right}{Column name for upper bound of observed values (default:
"right").}

\item{pwindow}{Column name for primary window (default: "pwindow").}

\item{D}{Column name for maximum delay (truncation point). If finite, the
distribution is truncated at D. If set to Inf, no truncation is applied.
(default: "D").}

\item{dprimary}{Function to generate the probability density function
(PDF) of primary event times. This function should take a value \code{x} and a
\code{pwindow} parameter, and return a probability density. It should be
normalized to integrate to 1 over [0, pwindow]. Defaults to a uniform
distribution over [0, pwindow]. Users can provide custom functions or use
helper functions like \code{dexpgrowth} for an exponential growth distribution.
See \code{\link[=pcd_primary_distributions]{pcd_primary_distributions()}} for examples. The package can identify
base R distributions for potential analytical solutions. For non-base R
functions, users can apply \code{\link[=add_name_attribute]{add_name_attribute()}} to yield properly tagged
functions if they wish to leverage analytical solutions.}

\item{dprimary_name}{\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} this argument will be
ignored in future versions; use \code{\link[=add_name_attribute]{add_name_attribute()}} on \code{dprimary}
instead}

\item{dprimary_args}{List of additional arguments to be passed to
dprimary. For example, when using \code{dexpgrowth}, you would
pass \code{list(min = 0, max = pwindow, r = 0.2)} to set the minimum, maximum,
and rate parameters}

\item{truncation_check_multiplier}{Numeric multiplier to use for checking
if the truncation time D is appropriate relative to the maximum delay.
Set to NULL to skip the check. Default is 2.}

\item{...}{Additional arguments to be passed to \code{\link[fitdistrplus:fitdist]{fitdistrplus::fitdist()}}.}
}
\value{
An object of class "fitdist" as returned by fitdistrplus::fitdist.
}
\description{
This function wraps the custom approach for fitting distributions to doubly
censored data using fitdistrplus and primarycensored.
}
\details{
This function temporarily assigns and then removes functions from the global
environment in order to work with fitdistr. Users should be aware of this
behaviour, especially if they have existing functions with the same names in
their global environment.
}
\examples{
\dontshow{if (requireNamespace("fitdistrplus", quietly = TRUE)) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
# Example with normal distribution
set.seed(123)
n <- 1000
true_mean <- 5
true_sd <- 2
pwindow <- 2
swindow <- 2
D <- 10
samples <- rprimarycensored(
  n, rnorm,
  mean = true_mean, sd = true_sd,
  pwindow = pwindow, swindow = swindow, D = D
)

delay_data <- data.frame(
  left = samples,
  right = samples + swindow,
  pwindow = rep(pwindow, n),
  D = rep(D, n)
)

fit_norm <- fitdistdoublecens(
  delay_data,
  distr = "norm",
  start = list(mean = 0, sd = 1)
)

summary(fit_norm)
\dontshow{\}) # examplesIf}
}
\seealso{
Modelling wrappers for external fitting packages
\code{\link{pcd_as_stan_data}()},
\code{\link{pcd_cmdstan_model}()}
}
\concept{modelhelpers}
