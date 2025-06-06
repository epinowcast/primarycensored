% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rprimarycensored.R
\name{rprimarycensored}
\alias{rprimarycensored}
\alias{rpcens}
\title{Generate random samples from a primary event censored distribution}
\usage{
rprimarycensored(
  n,
  rdist,
  pwindow = 1,
  swindow = 1,
  D = Inf,
  rprimary = stats::runif,
  rprimary_args = list(),
  oversampling_factor = 1.2,
  ...
)

rpcens(
  n,
  rdist,
  pwindow = 1,
  swindow = 1,
  D = Inf,
  rprimary = stats::runif,
  rprimary_args = list(),
  oversampling_factor = 1.2,
  ...
)
}
\arguments{
\item{n}{Number of random samples to generate.}

\item{rdist}{Function to generate random samples from the delay distribution
for example \code{\link[stats:Lognormal]{stats::rlnorm()}} for lognormal distribution.}

\item{pwindow}{Primary event window}

\item{swindow}{Integer specifying the window size for rounding the delay
(default is 1). If \code{swindow = 0} then no rounding is applied.}

\item{D}{Maximum delay (truncation point). If finite, the distribution is
truncated at D. If set to Inf, no truncation is applied. Defaults to Inf.}

\item{rprimary}{Function to generate random samples from the primary
distribution (default is \code{\link[stats:Uniform]{stats::runif()}}).}

\item{rprimary_args}{List of additional arguments to be passed to rprimary.}

\item{oversampling_factor}{Factor by which to oversample the number of
samples to account for truncation (default is 1.2).}

\item{...}{Additional arguments to be passed to the distribution function.}
}
\value{
Vector of random samples from the primary event censored
distribution censored by the secondary event window.
}
\description{
This function generates random samples from a primary event censored
distribution. It adjusts the distribution by accounting for the primary
event distribution and potential truncation at a maximum delay (D). The
function allows for custom primary event distributions and delay
distributions.
}
\details{
The mathematical formulation for generating random samples from a primary
event censored distribution is as follows:
\enumerate{
\item Generate primary event times (p) from the specified primary event
distribution (f_p) with parameters phi, defined between 0 and the primary
event window (pwindow):
\deqn{p \sim f_p(\phi), \quad p \in [0, pwindow]}
\item Generate delays (d) from the specified delay distribution (f_d) with
parameters theta:
\deqn{d \sim f_d(\theta)}
\item Calculate the total delays (t) by adding the primary event times and
the delays:
\deqn{t = p + d}
\item Apply truncation (i.e. remove any delays that fall outside the observation
window) to ensure that the delays are within the specified range [0, D],
where D is the maximum observable delay:
\deqn{t_{truncated} = \{t \mid 0 \leq t < D\}}
\item Round the truncated delays to the nearest secondary event window
(swindow):
\deqn{t_{valid} = \lfloor \frac{t_{truncated}}{swindow} \rfloor
     \times swindow}
}

The function oversamples to account for potential truncation and generates
additional samples if needed to reach the desired number of valid samples.
}
\examples{
# Example: Lognormal distribution with uniform primary events
rprimarycensored(10, rlnorm, meanlog = 0, sdlog = 1)

# Example: Lognormal distribution with exponential growth primary events
rprimarycensored(
  10, rlnorm,
  rprimary = rexpgrowth, rprimary_args = list(r = 0.2),
  meanlog = 0, sdlog = 1
)
}
\seealso{
Primary event censored distribution functions
\code{\link{dprimarycensored}()},
\code{\link{pprimarycensored}()},
\code{\link{qprimarycensored}()}
}
\concept{primarycensored}
