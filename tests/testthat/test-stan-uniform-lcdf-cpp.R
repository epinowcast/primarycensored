skip_on_cran()

# The C++ primarycensored_uniform_lcdf_vectorized() from pcd_stan_header()
# must return the same values as its Stan body and the same gradient as
# Stan's autodiff through that body. The probe model puts the parameters on
# the unconstrained scale, so there is no Jacobian term, and takes a
# weighted sum of the finite log CDFs as the target. The weights are
# positive, so the target does not cancel: far in the lower tail each
# gradient is large, and cancelling weights would turn their small relative
# errors into a large relative error of the sum.

cpp_probe_code <- "
data {
  int dist_id;
  int start;
  int n;
  real pwindow;
  vector[n] w;
}
parameters {
  vector[2] theta;
}
transformed parameters {
  array[2] real params;
  params[1] = dist_id == 1 ? theta[1] : exp(theta[1]);
  params[2] = exp(theta[2]);
}
model {
  vector[n] lcdf = primarycensored_uniform_lcdf_vectorized(
    start, n, dist_id, params, pwindow
  );
  for (d in start:n) {
    if (!is_inf(lcdf[d])) {
      target += w[d] * lcdf[d];
    }
  }
}
generated quantities {
  vector[n] lcdf = primarycensored_uniform_lcdf_vectorized(
    start, n, dist_id, params, pwindow
  );
}
"

cpp_cases <- list(
  list(dist_id = 1L, params = list(
    c(1.6, 0.42), c(1.5, 0.5), c(0.2, 1.2), c(-1, 0.3), c(3.5, 0.05),
    c(4, 2)
  )),
  list(dist_id = 2L, params = list(
    c(1.4, 0.38), c(2, 0.5), c(0.5, 0.2), c(20, 4), c(1.2, 0.02),
    c(50, 0.5), c(5, 1.2)
  )),
  list(dist_id = 3L, params = list(
    c(1.5, 5), c(0.7, 2), c(3, 30), c(5, 0.8), c(0.8, 3)
  ))
)

to_theta <- function(dist_id, params) {
  if (dist_id == 1L) c(params[1], log(params[2])) else log(params)
}

cpp_models <- function() {
  compile_cpp_and_stan(cpp_probe_code, "pcdcppprobe")
}

test_that("C++ uniform lcdf matches the Stan values", {
  models <- cpp_models()
  for (case in cpp_cases) {
    for (params in case$params) {
      for (setting in list(c(1, 40, 1), c(5, 40, 7), c(1, 70, 3))) {
        data <- list(
          dist_id = case$dist_id, start = setting[1], n = setting[2],
          pwindow = setting[3], w = rep(1, setting[2])
        )
        init <- list(theta = to_theta(case$dist_id, params))
        stan <- stan_gq(models$stan, data, init, "lcdf")
        cpp <- stan_gq(models$cpp, data, init, "lcdf")
        computed <- seq(setting[1], setting[2])
        info <- sprintf(
          "dist %d, params (%s), start %d, n %d, pwindow %d",
          case$dist_id, toString(params), setting[1], setting[2], setting[3]
        )
        expect_identical(
          is.infinite(cpp[computed]), is.infinite(stan[computed]),
          info = info
        )
        finite <- computed[is.finite(stan[computed])]
        expect_equal(cpp[finite], stan[finite], tolerance = 1e-13, info = info)
        expect_true(all(cpp[finite] <= 0), info = info)
      }
    }
  }
})

# Below a log CDF of about -700 the two terms are close to the smallest
# doubles and their difference has few significant digits, so both
# gradients are rounding error there. Those delays get zero weight.
test_that("C++ uniform lcdf gradient matches Stan autodiff", {
  skip_if_no_log_prob()
  models <- cpp_models()
  set.seed(1)
  for (case in cpp_cases) {
    for (params in case$params) {
      for (setting in list(c(1, 40, 1), c(5, 40, 7), c(1, 70, 3))) {
        data <- list(
          dist_id = case$dist_id, start = setting[1], n = setting[2],
          pwindow = setting[3], w = stats::runif(setting[2], 0.5, 1.5)
        )
        upars <- list(theta = to_theta(case$dist_id, params))
        lcdf <- stan_gq(models$stan, data, upars, "lcdf")
        data$w[!is.na(lcdf) & lcdf < -700] <- 0
        stan <- stan_log_prob(models$stan, data, upars)
        cpp <- stan_log_prob(models$cpp, data, upars)
        info <- sprintf(
          "dist %d, params (%s), start %d, n %d, pwindow %d",
          case$dist_id, toString(params), setting[1], setting[2], setting[3]
        )
        expect_equal(cpp$lp, stan$lp, tolerance = 1e-12, info = info)
        expect_equal(cpp$grad, stan$grad, tolerance = 1e-8, info = info)
      }
    }
  }
})

test_that("C++ uniform lcdf gradient matches finite differences", {
  skip_if_no_log_prob()
  models <- cpp_models()
  set.seed(2)
  for (case in cpp_cases) {
    params <- case$params[[1]]
    data <- list(
      dist_id = case$dist_id, start = 1, n = 30, pwindow = 1,
      w = stats::runif(30, 0.5, 1.5)
    )
    theta <- to_theta(case$dist_id, params)
    cpp <- stan_log_prob(models$cpp, data, list(theta = theta))
    h <- 1e-5
    fd <- vapply(1:2, function(i) {
      shift <- replace(numeric(2), i, h)
      (stan_log_prob(models$cpp, data, list(theta = theta + shift))$lp -
        stan_log_prob(models$cpp, data, list(theta = theta - shift))$lp) /
        (2 * h)
    }, numeric(1))
    expect_equal(cpp$grad, fd, tolerance = 1e-6, info = toString(params))
  }
})

# Lognormal log CDF and its gradient with respect to (mu, log(sigma)) at d,
# by integrating the delay CDF and its derivatives over the primary window.
# This does not cancel in the lower tail.
lognormal_uniform_lcdf_ref <- function(d, mu, sigma, pwindow) {
  z <- function(t) (log(pmax(t, 1e-300)) - mu) / sigma
  cdf <- function(t) ifelse(t > 0, stats::plnorm(t, mu, sigma), 0)
  dmu <- function(t) ifelse(t > 0, -stats::dnorm(z(t)) / sigma, 0)
  dsigma <- function(t) ifelse(t > 0, -z(t) * stats::dnorm(z(t)), 0)
  over_window <- function(f) {
    stats::integrate(
      function(u) f(d - u), 0, pwindow,
      rel.tol = 1e-13, abs.tol = 0, subdivisions = 1000L
    )$value
  }
  den <- over_window(cdf)
  c(log(den / pwindow), over_window(dmu) / den, over_window(dsigma) / den)
}

test_that("C++ lognormal uniform lcdf gradient is accurate in the lower
  tail", {
  skip_if_no_log_prob()
  models <- cpp_models()
  mu <- 3.5
  sigma <- 0.05
  # The log CDF runs from about -590 at d = 6 to about -2e-4 at d = 40
  for (d in c(6, 10, 20, 30, 40)) {
    data <- list(
      dist_id = 1L, start = 1, n = 40, pwindow = 1,
      w = replace(numeric(40), d, 1)
    )
    cpp <- stan_log_prob(
      models$cpp, data, list(theta = c(mu, log(sigma)))
    )
    ref <- lognormal_uniform_lcdf_ref(d, mu, sigma, 1)
    expect_equal(cpp$lp, ref[1], tolerance = 1e-10, info = paste("d =", d))
    expect_equal(cpp$grad, ref[2:3], tolerance = 1e-9, info = paste("d =", d))
  }
})

test_that("C++ uniform lcdf rejects a gamma rate that is not finite", {
  skip_if_no_log_prob()
  models <- cpp_models()
  data <- list(dist_id = 2L, start = 1, n = 10, pwindow = 1, w = rep(1, 10))
  upars <- list(theta = c(log(2), 800))
  for (model in models) {
    err <- tryCatch(
      stan_log_prob(model, data, upars),
      error = conditionMessage
    )
    expect_type(err, "character")
    expect_match(err, "rate must be positive finite")
  }
})
