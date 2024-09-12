skip_on_cran()
if (on_ci()) {
  skip_on_os("windows")
  skip_on_os("mac")
}


test_that("Stan primary_censored_dist_cdf matches R pprimarycensoreddist", {
  d_values <- list(
    seq(0, 10, by = 0.5),
    seq(0, 5, by = 0.25),
    seq(0, 15, by = 1)
  )
  D_values <- list(Inf, 10, 15)

  for (d in d_values) {
    for (D in D_values) {
      dist_id <- 1 # Lognormal
      params <- c(0, 1) # meanlog, sdlog
      pwindow <- 1
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

      expect_equal(
        stan_cdf, r_cdf,
        tolerance = 1e-6,
        info = paste(
          "Failed for d =", paste(range(d), collapse = "-"),
          "and D =", D
        )
      )
    }
  }
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
    r_cdf <- pprimarycensoreddist(
      d, plnorm,
      pwindow = pwindow, D = D, meanlog = params[1],
      sdlog = params[2]
    )

    expect_equal(exp(stan_lcdf), r_cdf, tolerance = 1e-6)
  }
)

test_that(
  "Stan primary_censored_dist_lpmf throws an error for invalid upper truncation
   point",
  {
    d <- 10
    dist_id <- 1 # Lognormal
    params <- c(0, 1) # meanlog, sdlog
    pwindow <- 1
    d_upper <- 11
    D <- 10
    primary_dist_id <- 1 # Uniform
    primary_params <- numeric(0)

    expect_error(
      primary_censored_dist_lpmf(
        d, dist_id, params, pwindow, d_upper, D, primary_dist_id, primary_params
      ),
      "Upper truncation point is greater than D"
    )
  }
)

test_that(
  "Stan primary_censored_dist matches R primarycensoreddist when d is the same
   as D - 1",
  {
    dist_id <- 1 # Lognormal
    params <- c(0, 1) # meanlog, sdlog
    pwindow <- 1
    primary_dist_id <- 1 # Uniform
    primary_params <- numeric(0)

    d_values <- 1:10
    D_values <- c(Inf, 11) # Test without and with truncation

    for (D in D_values) {
      truncation_label <- if (is.infinite(D)) {
        "without truncation"
      } else {
        "with truncation"
      }
      for (d in d_values) {
        stan_pmf <- primary_censored_dist_pmf(
          d, dist_id, params, pwindow, d + 1, D, primary_dist_id,
          primary_params
        )
        r_pmf <- dprimarycensoreddist(
          d, plnorm,
          pwindow = pwindow, swindow = 1, D = D,
          meanlog = params[1], sdlog = params[2]
        )
        expect_equal(
          stan_pmf, r_pmf,
          tolerance = 1e-6,
          info = paste(
            "Mismatch for d =", d, truncation_label, "\n",
            "Stan PMF:", stan_pmf, "\n",
            "R PMF:", r_pmf, "\n",
            "Difference:", abs(stan_pmf - r_pmf)
          )
        )
      }
    }
  }
)

test_that("Stan primary_censored_dist_pmf matches R dprimarycensoreddist", {
  d <- 0:10
  dist_id <- 1 # Lognormal
  params <- c(0, 1) # meanlog, sdlog
  pwindow <- 1
  d_upper <- d + 1
  D <- Inf
  primary_dist_id <- 1 # Uniform
  primary_params <- numeric(0)

  stan_pmf <- mapply(
    primary_censored_dist_pmf,
    d = d,
    d_upper = d + 1,
    MoreArgs = list(
      dist_id = dist_id,
      params = params,
      pwindow = pwindow,
      D = D,
      primary_dist_id = primary_dist_id,
      primary_params = primary_params
    )
  )
  r_pmf <- dprimarycensoreddist(
    d, plnorm,
    pwindow = pwindow, swindow = 1, D = D,
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
    d_upper <- d + 1
    D <- Inf
    primary_dist_id <- 1 # Uniform
    primary_params <- numeric(0)

    stan_lpmf <- mapply(
      primary_censored_dist_lpmf,
      d = d,
      d_upper = d_upper,
      MoreArgs = list(
        dist_id = dist_id,
        params = params,
        pwindow = pwindow,
        D = D,
        primary_dist_id = primary_dist_id,
        primary_params = primary_params
      )
    )
    r_lpmf <- log(
      dprimarycensoreddist(
        d, plnorm,
        pwindow = pwindow, swindow = 1, D = D,
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

    stan_pmf <- primary_censored_sone_pmf_vectorized(
      max_delay, D, dist_id, params, pwindow, primary_dist_id, primary_params
    )
    r_pmf <- dprimarycensoreddist(
      0:max_delay, plnorm,
      pwindow = pwindow, swindow = 1, D = D,
      meanlog = params[1], sdlog = params[2]
    )

    expect_equal(stan_pmf, r_pmf, tolerance = 1e-6)
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

    stan_pmf <- primary_censored_sone_pmf_vectorized(
      max_delay, D, dist_id, params, pwindow, primary_dist_id, primary_params
    )
    r_pmf <- dprimarycensoreddist(
      0:max_delay, plnorm,
      pwindow = pwindow, swindow = 1, D = D,
      meanlog = params[1], sdlog = params[2]
    )

    expect_equal(stan_pmf, r_pmf, tolerance = 1e-6)
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

    stan_pmf <- primary_censored_sone_pmf_vectorized(
      max_delay, D, dist_id, params, pwindow, primary_dist_id, primary_params
    )
    r_pmf <- dprimarycensoreddist(
      0:max_delay, plnorm,
      pwindow = pwindow, swindow = 1, D = D,
      meanlog = params[1], sdlog = params[2]
    )

    expect_equal(stan_pmf, r_pmf, tolerance = 1e-6)
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

    stan_lpmf <- primary_censored_sone_lpmf_vectorized(
      max_delay, D, dist_id, params, pwindow, primary_dist_id, primary_params
    )
    r_lpmf <- dprimarycensoreddist(
      0:max_delay, plnorm,
      pwindow = pwindow, swindow = 1, D = D,
      meanlog = params[1], sdlog = params[2], log = TRUE
    )

    expect_equal(stan_lpmf, r_lpmf, tolerance = 1e-6)
  }
)
