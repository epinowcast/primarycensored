# Quadrature of the delay CDF over the primary window, split at the kinks
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

# Elementwise relative error, with values at or below `floor` compared as tails
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

# Log quadrature relative to F_T(d), so it holds where the CDF underflows
reference_uniform_lcdf <- function(lp, d, pwindow, lower = 0) {
  vapply(d, function(di) {
    q_lower <- max(di - pwindow, lower)
    lp_d <- lp(di)
    breaks <- sort(unique(c(q_lower, di - (di - q_lower) * 2^-(0:40), di)))
    breaks <- breaks[breaks >= q_lower & breaks <= di]
    total <- 0
    for (i in seq_len(length(breaks) - 1L)) {
      total <- total + stats::integrate(
        function(t) exp(vapply(t, lp, numeric(1)) - lp_d),
        breaks[i], breaks[i + 1L],
        rel.tol = 1e-13, abs.tol = 0, subdivisions = 500L,
        stop.on.error = FALSE
      )$value
    }
    lp_d + log(total) - log(pwindow)
  }, numeric(1))
}
