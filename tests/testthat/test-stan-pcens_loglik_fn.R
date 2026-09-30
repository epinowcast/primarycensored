skip_on_cran()

# `pars` are the R parameters of `pdist` and `params` the Stan ones. The
# non-uniform primary uses the Stan ODE solver, so has a looser tolerance.
stan_cases <- list(
  gamma_unif = list(
    pdist = pgamma, dprimary = dunif, primary_args = list(),
    pars = list(shape = 2.5, rate = 0.8),
    dist_id = 2L, params = c(2.5, 0.8),
    primary_id = 1L, primary_params = numeric(0), tol = 1e-6
  ),
  lnorm_unif = list(
    pdist = plnorm, dprimary = dunif, primary_args = list(),
    pars = list(meanlog = 1.2, sdlog = 0.6),
    dist_id = 1L, params = c(1.2, 0.6),
    primary_id = 1L, primary_params = numeric(0), tol = 1e-6
  ),
  weibull_unif = list(
    pdist = pweibull, dprimary = dunif, primary_args = list(),
    pars = list(shape = 1.8, scale = 4),
    dist_id = 3L, params = c(1.8, 4),
    primary_id = 1L, primary_params = numeric(0), tol = 1e-6
  ),
  gamma_expgrowth = list(
    pdist = pgamma, dprimary = dexpgrowth, primary_args = list(r = 0.2),
    pars = list(shape = 2.5, rate = 0.8),
    dist_id = 2L, params = c(2.5, 0.8),
    primary_id = 2L, primary_params = 0.2, tol = 1e-3
  )
)

test_that("pcens_loglik_fn agrees with the Stan log PMF", {
  x <- c(0, 1, 2, 3, 5, 8)
  for (sc in stan_cases) {
    ll <- pcens_loglik_fn(
      x, sc$pdist,
      pwindow = 2, swindow = 1, L = 0, D = 10,
      dprimary = sc$dprimary, primary_args = sc$primary_args
    )
    stan_ll <- vapply(x, function(d) {
      primarycensored_lpmf( # nolint: object_usage_linter.
        as.integer(d), sc$dist_id, sc$params, 2, d + 1, 0, 10,
        sc$primary_id, sc$primary_params
      )
    }, numeric(1))
    expect_equal(do.call(ll, sc$pars), stan_ll, tolerance = sc$tol)
  }
})
