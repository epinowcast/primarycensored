skip_on_cran()
skip_on_os("windows")
skip_on_os("mac")

test_that("Stan primary_censored_dist_cdf matches R pprimarycensoreddist", {
  d <- seq(0, 10, by = 0.5)
  dist_id <- 1 # Lognormal
  params <- c(0, 1) # meanlog, sdlog
  pwindow <- 1
  D <- Inf
  primary_dist_id <- 1 # Uniform
  primary_params <- array(numeric(0))

  stan_cdf <- sapply(
    d, primary_censored_dist_cdf, dist_id, params, pwindow, D,
    primary_dist_id, primary_params
  )
  r_cdf <- pprimarycensoreddist(
    d, plnorm,
    pwindow = pwindow, D = D, meanlog = params[1], sdlog = params[2]
  )

  expect_equal(stan_cdf, r_cdf, tolerance = 1e-6)
})

test_that(
  "Stan primary_censored_dist_lcdf matches R pprimarycensoreddist with
   log.p = TRUE",
  {
    d <- seq(0, 10, by = 0.5)
    dist_id <- 1 # Lognormal
    params <- c(0, 1) # meanlog, sdlog
    pwindow <- 1
    D <- Inf
    primary_dist_id <- 1 # Uniform
    primary_params <- numeric(0)

    stan_lcdf <- sapply(
      d, primary_censored_dist_lcdf, dist_id, params, pwindow, D,
      primary_dist_id, primary_params
    )
    r_lcdf <- log(
      pprimarycensoreddist(
        d, plnorm,
        pwindow = pwindow, D = D, meanlog = params[1],
        sdlog = params[2]
      )
    )

    expect_equal(stan_lcdf, r_lcdf, tolerance = 1e-6)
  }
)

test_that("Stan primary_censored_dist_pmf matches R dprimarycensoreddist", {
  d <- 0:10
  dist_id <- 1 # Lognormal
  params <- c(0, 1) # meanlog, sdlog
  pwindow <- 1
  swindow <- 1
  D <- Inf
  primary_dist_id <- 1 # Uniform
  primary_params <- numeric(0)

  stan_pmf <- sapply(
    d, primary_censored_dist_pmf, dist_id, params, pwindow, swindow, D,
    primary_dist_id, primary_params
  )
  r_pmf <- dprimarycensoreddist(
    d, plnorm,
    pwindow = pwindow, swindow = swindow, D = D,
    meanlog = params[1], sdlog = params[2]
  )

  expect_equal(stan_pmf, r_pmf, tolerance = 1e-6)
})

test_that(
  "Stan primary_censored_dist_lpmf matches R dprimarycensoreddist with
   log = TRUE",
  {
    d <- 0:10
    dist_id <- 1 # Lognormal
    params <- c(0, 1) # meanlog, sdlog
    pwindow <- 1
    swindow <- 1
    D <- Inf
    primary_dist_id <- 1 # Uniform
    primary_params <- numeric(0)

    stan_lpmf <- sapply(
      d, primary_censored_dist_lpmf, dist_id, params, pwindow, swindow, D,
      primary_dist_id, primary_params
    )
    r_lpmf <- log(
      dprimarycensoreddist(
        d, plnorm,
        pwindow = pwindow, swindow = swindow, D = D,
        meanlog = params[1], sdlog = params[2]
      )
    )

    expect_equal(stan_lpmf, r_lpmf, tolerance = 1e-6)
  }
)

test_that(
  "Stan primary_censored_sone_pmf_vectorized matches R dprimarycensoreddist",
  {
    max_delay <- 10
    dist_id <- 1 # Lognormal
    params <- c(0, 1) # meanlog, sdlog
    pwindow <- 1
    D <- Inf
    primary_dist_id <- 1 # Uniform
    primary_params <- numeric(0)

    stan_pmf_approx <- primary_censored_sone_pmf_vectorized(
      max_delay, D, dist_id, params, pwindow, primary_dist_id, primary_params, 1
    )
    stan_pmf_exact <- primary_censored_sone_pmf_vectorized(
      max_delay, D, dist_id, params, pwindow, primary_dist_id, primary_params, 0
    )
    r_pmf <- dprimarycensoreddist(
      0:max_delay, plnorm,
      pwindow = pwindow, swindow = 1, D = D,
      meanlog = params[1], sdlog = params[2]
    )

    expect_equal(stan_pmf_approx, r_pmf, tolerance = 1e-6)
    expect_equal(stan_pmf_exact, r_pmf, tolerance = 1e-6)
  }
)

test_that(
  "Stan primary_censored_sone_pmf_vectorized matches R dprimarycensoreddist
   with finite D",
  {
    max_delay <- 10
    dist_id <- 1 # Lognormal
    params <- c(0, 1) # meanlog, sdlog
    pwindow <- 1
    D <- 15
    primary_dist_id <- 1 # Uniform
    primary_params <- numeric(0)

    stan_pmf_approx <- primary_censored_sone_pmf_vectorized(
      max_delay, D, dist_id, params, pwindow, primary_dist_id, primary_params, 1
    )
    stan_pmf_exact <- primary_censored_sone_pmf_vectorized(
      max_delay, D, dist_id, params, pwindow, primary_dist_id, primary_params, 0
    )
    r_pmf <- dprimarycensoreddist(
      0:max_delay, plnorm,
      pwindow = pwindow, swindow = 1, D = D,
      meanlog = params[1], sdlog = params[2]
    )

    expect_equal(stan_pmf_approx, r_pmf, tolerance = 1e-2)
    expect_equal(stan_pmf_exact, r_pmf, tolerance = 1e-6)
    expect_true(all(abs(stan_pmf_exact - r_pmf) < abs(stan_pmf_approx - r_pmf)))
  }
)

test_that(
  "Stan primary_censored_sone_pmf_vectorized matches R dprimarycensoreddist
   with D equal to max_delay + 1",
  {
    max_delay <- 10
    dist_id <- 1 # Lognormal
    params <- c(0, 1) # meanlog, sdlog
    pwindow <- 1
    D <- max_delay + 1
    primary_dist_id <- 1 # Uniform
    primary_params <- numeric(0)

    stan_pmf_approx <- primary_censored_sone_pmf_vectorized(
      max_delay, D, dist_id, params, pwindow, primary_dist_id, primary_params, 1
    )
    stan_pmf_exact <- primary_censored_sone_pmf_vectorized(
      max_delay, D, dist_id, params, pwindow, primary_dist_id, primary_params, 0
    )
    r_pmf <- dprimarycensoreddist(
      0:max_delay, plnorm,
      pwindow = pwindow, swindow = 1, D = D,
      meanlog = params[1], sdlog = params[2]
    )

    expect_equal(stan_pmf_approx, r_pmf, tolerance = 1e-3)
    expect_equal(stan_pmf_exact, r_pmf, tolerance = 1e-6)
    expect_true(all(abs(stan_pmf_exact - r_pmf) < abs(stan_pmf_approx - r_pmf)))
  }
)

test_that(
  "Stan primary_censored_sone_lpmf_vectorized matches R dprimarycensoreddist
   with log = TRUE",
  {
    max_delay <- 10
    dist_id <- 1 # Lognormal
    params <- c(0, 1) # meanlog, sdlog
    pwindow <- 1
    D <- Inf
    primary_dist_id <- 1 # Uniform
    primary_params <- numeric(0)

    stan_lpmf_approx <- primary_censored_sone_lpmf_vectorized(
      max_delay, D, dist_id, params, pwindow, primary_dist_id, primary_params, 1
    )
    stan_lpmf_exact <- primary_censored_sone_lpmf_vectorized(
      max_delay, D, dist_id, params, pwindow, primary_dist_id, primary_params, 0
    )
    r_lpmf <- dprimarycensoreddist(
      0:max_delay, plnorm,
      pwindow = pwindow, swindow = 1, D = D,
      meanlog = params[1], sdlog = params[2], log = TRUE
    )

    expect_equal(stan_lpmf_approx, r_lpmf, tolerance = 1e-6)
    expect_equal(stan_lpmf_exact, r_lpmf, tolerance = 1e-8)
    expect_true(all(abs(stan_lpmf_exact - r_lpmf) < abs(stan_lpmf_approx - r_lpmf)))
  }
)
