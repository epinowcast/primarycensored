# High accuracy references for the uniform primary event censored CDF
#
# F_{S+}(d) = int_{d - pwindow}^{d} F_T(t) dt / pwindow is integrated with
# stats::integrate at a tolerance near double precision. The integrand has a
# kink wherever F_T changes form (t = 0 for non-negative delays, t = 1 for the
# Beta and t = y_min for the Pareto), so the integral is split there. The
# default stats::integrate tolerance of the package numerical path is not
# enough to check the analytical solutions to a relative 1e-9.
reference_uniform_cdf <- function(pdist, d, pwindow, kinks = numeric(0)) {
  vapply(d, function(di) {
    lower <- di - pwindow
    inside <- kinks[kinks > lower & kinks < di]
    breaks <- sort(unique(c(lower, inside, di)))
    total <- 0
    for (i in seq_len(length(breaks) - 1L)) {
      total <- total + stats::integrate(
        function(t) vapply(t, pdist, numeric(1)),
        breaks[i], breaks[i + 1L],
        rel.tol = 1.2e-14, abs.tol = 0, subdivisions = 500L,
        stop.on.error = FALSE
      )$value
    }
    total / pwindow
  }, numeric(1))
}

# Checks the relative error elementwise. Values of `expected` at or below
# `floor` are in the far tail, where the reference underflows, and need
# `actual` to be at or below `floor` too.
expect_rel_equal <- function(actual, expected, tolerance = 1e-9,
                             floor = 1e-250, info = NULL) {
  testthat::expect_length(actual, length(expected))
  keep <- expected > floor
  rel <- abs(actual[keep] - expected[keep]) / expected[keep]
  testthat::expect_true(all(is.finite(rel)), info = info)
  testthat::expect_true(
    all(rel <= tolerance),
    info = paste(
      info, "max relative error", format(suppressWarnings(max(rel, 0)))
    )
  )
  testthat::expect_true(
    all(actual[!keep] <= floor * 1e3),
    info = paste(info, "far tail above the floor")
  )
}
