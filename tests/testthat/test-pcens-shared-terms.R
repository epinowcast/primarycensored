shared_delays <- list(
  list(pdist = pgamma, args = list(shape = 3, scale = 2)),
  list(pdist = pgamma, args = list(shape = 0.6, rate = 0.3)),
  list(pdist = plnorm, args = list(meanlog = 1.5, sdlog = 0.5)),
  list(pdist = plnorm, args = list(meanlog = 0.2, sdlog = 1.2)),
  list(pdist = pweibull, args = list(shape = 1.5, scale = 5)),
  list(pdist = pweibull, args = list(shape = 4, scale = 2))
)

flexsurv_delays <- function() {
  skip_if_not_installed("flexsurv")
  list(
    list(
      pdist = flexsurv::pgengamma.orig,
      args = list(shape = 1.5, scale = 3, k = 2)
    ),
    list(
      pdist = flexsurv::pgengamma,
      args = list(mu = 1, sigma = 0.6, Q = 0.8)
    )
  )
}

shared_obj <- function(delay) {
  do.call(new_pcens, c(list(delay$pdist, dunif), delay$args))
}

# Delay CDF evaluated for each point on its own
per_point_cdf <- function(obj, q, pwindow) {
  vapply(q, function(x) pcens_cdf(obj, x, pwindow), numeric(1))
}

test_that(".pcens_terms_spec gives terms for registered pairs only", {
  for (delay in c(shared_delays, flexsurv_delays())) {
    spec <- .pcens_terms_spec(shared_obj(delay))
    expect_type(spec$terms, "closure")
    expect_type(spec$combine, "closure")
    expect_identical(dim(spec$terms(c(0, 1, 2.5))), c(3L, 2L))
    expect_identical(spec$terms(c(0, -1)), matrix(0, 2, 2))
  }
  expect_null(.pcens_terms_spec(new_pcens(
    pgamma, dexpgrowth,
    primary_args = list(r = 0.2), shape = 2, scale = 1
  )))
  expect_null(.pcens_terms_spec(new_pcens(pexp, dunif, rate = 1)))
  skip_if_not_installed("flexsurv")
  expect_null(.pcens_terms_spec(new_pcens(
    flexsurv::pgengamma, dunif, mu = 1, sigma = 0.6, Q = -0.5
  )))
})

test_that(".pcens_cdf_shared evaluates each distinct endpoint once", {
  evaluated <- 0
  spec <- list(
    terms = function(t) {
      evaluated <<- evaluated + length(t)
      cbind(t, t^2)
    },
    combine = function(td, tq, pwindow) (td[, 1] - tq[, 1]) / pwindow
  )
  q <- 0:10
  .pcens_cdf_shared(spec, q, 2)
  expect_identical(evaluated, length(unique(c(q, pmax(q - 2, 0)))))

  evaluated <- 0
  .pcens_cdf_shared(spec, q, 1.5)
  expect_identical(evaluated, length(unique(c(q, pmax(q - 1.5, 0)))))

  # Nothing overlaps, so each endpoint is evaluated directly
  evaluated <- 0
  .pcens_cdf_shared(spec, c(2.3, 7.9), 1)
  expect_identical(evaluated, 4)
})

test_that("shared terms match single-point evaluation", {
  grids <- list(
    list(q = 0:40, pwindow = 1),
    list(q = 0:40, pwindow = 3),
    list(q = seq(0, 20, by = 0.5), pwindow = 1.5),
    list(q = c(-2, -0.5, 0, 0.2, 0.5, 1, 1.5, 4, 4, 30), pwindow = 1),
    list(q = seq(0, 6, by = 0.25), pwindow = 0.25)
  )
  for (delay in c(shared_delays, flexsurv_delays())) {
    obj <- shared_obj(delay)
    for (grid in grids) {
      shared <- pcens_cdf(obj, grid$q, grid$pwindow)
      expect_equal(
        shared, per_point_cdf(obj, grid$q, grid$pwindow),
        tolerance = 1e-14
      )
    }
  }
})

# Integrate the delay CDF against the uniform primary at high accuracy
reference_cdf <- function(delay, q, pwindow) {
  vapply(q, function(d) {
    stats::integrate(
      function(p) do.call(delay$pdist, c(list(d - p), delay$args)),
      lower = 0, upper = pwindow, rel.tol = 1e-12, abs.tol = 0,
      subdivisions = 1000L
    )$value / pwindow
  }, numeric(1))
}

test_that("shared terms match high accuracy integration", {
  q <- c(0, 0.3, 1, 2, 5, 12, 40)
  for (delay in c(shared_delays, flexsurv_delays())) {
    obj <- shared_obj(delay)
    for (pwindow in c(1, 2.5)) {
      expect_equal(
        pcens_cdf(obj, q, pwindow),
        reference_cdf(delay, q, pwindow),
        tolerance = 1e-9
      )
    }
  }
})

test_that("pprimarycensored and dprimarycensored share terms on a lattice", {
  x <- 0:25
  for (delay in c(shared_delays, flexsurv_delays())) {
    for (pwindow in c(1, 2, 3)) {
      args <- c(list(pdist = delay$pdist, pwindow = pwindow), delay$args)
      p <- do.call(pprimarycensored, c(list(q = x), args))
      expect_equal(
        p,
        vapply(
          x, function(i) do.call(pprimarycensored, c(list(q = i), args)),
          numeric(1)
        ),
        tolerance = 1e-14
      )
      d <- do.call(dprimarycensored, c(list(x = x, D = 26), args))
      p_d <- do.call(pprimarycensored, c(list(q = x, L = 0, D = 26), args))
      p_up <- do.call(pprimarycensored, c(list(q = x + 1, L = 0, D = 26), args))
      expect_equal(d, p_up - p_d, tolerance = 1e-12)
    }
  }
})

test_that("shared terms handle empty input", {
  obj <- shared_obj(shared_delays[[1]])
  expect_identical(pcens_cdf(obj, numeric(0), 1), numeric(0))
})
