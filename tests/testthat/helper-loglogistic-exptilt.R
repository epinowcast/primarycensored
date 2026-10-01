# Helpers for the log-logistic delay tests with a tilted primary.

# Density of the log-logistic with a shape and a scale
dllogis_test <- function(u, shape, scale) {
  z <- shape * (log(u) - log(scale))
  exp(log(shape / u) + z - 2 * log1p(exp(z)))
}

# Reference transform int_0^t exp(xi u) f(u) du by integrating over the CDF
# scale p = F(u), where the integrand is smooth
loglogistic_transform_ref <- function(t, xi, shape, scale) {
  cdf_t <- pllogis_test(t, shape, scale)
  stats::integrate(
    function(p) exp(xi * scale * (p / (1 - p))^(1 / shape)),
    0, cdf_t,
    rel.tol = 1e-13, abs.tol = 0, subdivisions = 5000L
  )$value
}

# Cases where the combination of the lower transforms is ill conditioned
# near the series limit, for a positive tilt. Each is (shape, scale, rho,
# pwindow).
conditioning_cases <- list(
  c(0.325, 94.8, 0.0103, 0.3),
  c(8, 10, 0.2, 7),
  c(5, 5, 0.3, 3),
  c(4, 8, 0.25, 1),
  c(10.3, 31.9, 0.3026, 0.05),
  c(2, 5, 0.4, 2),
  c(8, 10, -0.2, 7),
  c(4, 8, -1, 1)
)

# Reference PMF of the integer delays 0:(n - 1) for a case of
# `conditioning_cases`, from differences of the smaller tail
loglogistic_pmf_reference <- function(case, n) {
  ref <- vapply(
    0:n, loglogistic_censored_reference, numeric(2),
    shape = case[1], scale = case[2], rho = case[3], pwindow = case[4]
  )
  ifelse(
    ref["cdf", -1] > 0.5,
    -diff(ref["survival", ]),
    diff(ref["cdf", ])
  )
}

# Cases in the upper tail of the tilted CDF, where the smaller tail is below
# 1e-7, each (shape, scale, pwindow, rho, q). They are beyond the series
# limit or ill conditioned, so they use the numerical CDF.
upper_tail_cases <- list(
  c(99, 9.38, 0.162, -0.00155, 11.49),
  c(14.6, 1.07, 7.5, 0.667, 10.64),
  c(3.4, 0.4, 14.6, -0.18, 58.3),
  c(1.57, 1.2, 0.31, -2.16, 3e4),
  c(1.5, 5, 1, 0.2, 1e4)
)

# Cases of the PMF in the upper tail at a large delay, each (shape, scale,
# pwindow, rho, d)
upper_tail_pmf_cases <- list(
  c(0.5, 5, 1, 0.3, 3e4),
  c(1.5, 5, 1, 0.2, 3000),
  c(3, 5, 1, -0.5, 500),
  c(8, 5, 1, -0.3, 40),
  c(1.5, 5, 1, 0.2, 1e5)
)

# Reference PMF at the integer delay d, from the survival so the upper tail
# keeps its relative precision
loglogistic_pmf_at <- function(case) {
  ref <- vapply(
    case[5] + 0:1, loglogistic_censored_reference, numeric(2),
    shape = case[1], scale = case[2], rho = case[4], pwindow = case[3]
  )
  -diff(ref["survival", ])
}
