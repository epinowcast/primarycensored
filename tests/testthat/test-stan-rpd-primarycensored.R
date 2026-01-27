skip_on_cran()

# Left truncation (L parameter) tests

test_that("Stan primarycensored_cdf matches R pprimarycensored with L > 0", {
  d_values <- seq(0, 10, by = 0.5)
  D <- 12
  L <- 2

  dist_id <- 1 # Lognormal
  params <- c(1, 0.5) # meanlog, sdlog
  pwindow <- 1
  primary_id <- 1 # Uniform
  primary_params <- array(numeric(0))

  stan_cdf <- sapply(
    d_values, primarycensored_cdf, dist_id, params, pwindow, L, D,
    primary_id, primary_params
  )
  r_cdf <- pprimarycensored(
    d_values, plnorm,
    pwindow = pwindow, D = D, L = L, meanlog = params[1], sdlog = params[2]
  )

  expect_equal(stan_cdf, r_cdf, tolerance = 1e-6)
})

test_that("Stan primarycensored_lcdf matches R pprimarycensored with L > 0", {
  d_values <- seq(0, 10, by = 0.5)
  D <- 12
  L <- 2

  dist_id <- 1 # Lognormal
  params <- c(1, 0.5) # meanlog, sdlog
  pwindow <- 1
  primary_id <- 1 # Uniform
  primary_params <- array(numeric(0))

  stan_lcdf <- sapply(
    d_values, primarycensored_lcdf, dist_id, params, pwindow, L, D,
    primary_id, primary_params
  )
  r_cdf <- pprimarycensored(
    d_values, plnorm,
    pwindow = pwindow, D = D, L = L, meanlog = params[1], sdlog = params[2]
  )

  expect_equal(exp(stan_lcdf), r_cdf, tolerance = 1e-6)
})

test_that("Stan primarycensored_lpmf matches R dprimarycensored with L > 0", {
  D <- 12
  L <- 2
  d <- L:(D - 1)

  dist_id <- 1 # Lognormal
  params <- c(1, 0.5) # meanlog, sdlog
  pwindow <- 1
  primary_id <- 1 # Uniform
  primary_params <- numeric(0)

  stan_lpmf <- mapply(
    primarycensored_lpmf,
    d = d,
    d_upper = d + 1,
    MoreArgs = list(
      dist_id = dist_id,
      params = params,
      pwindow = pwindow,
      L = L,
      D = D,
      primary_id = primary_id,
      primary_params = primary_params
    )
  )
  r_lpmf <- dprimarycensored(
    d, plnorm,
    pwindow = pwindow, swindow = 1, D = D, L = L,
    meanlog = params[1], sdlog = params[2], log = TRUE
  )

  expect_equal(stan_lpmf, r_lpmf, tolerance = 1e-6)
})

test_that(
  "Stan primarycensored_sone_lpmf_vectorized matches R dprimarycensored
   with L > 0",
  {
    max_delay <- 10
    D <- 15
    L <- 2

    dist_id <- 1 # Lognormal
    params <- c(1, 0.5) # meanlog, sdlog
    pwindow <- 1
    primary_id <- 1 # Uniform
    primary_params <- numeric(0)

    stan_lpmf <- primarycensored_sone_lpmf_vectorized(
      max_delay, L, D, dist_id, params, pwindow, primary_id, primary_params
    )

    # Stan returns -Inf for indices <= L, R errors for x < L
    # Check that Stan correctly returns -Inf for excluded delays
    expect_true(all(is.infinite(stan_lpmf[1:L])))
    expect_true(all(stan_lpmf[1:L] < 0))

    # Compare valid range with R
    valid_indices <- (L + 1):(max_delay + 1)
    r_lpmf <- dprimarycensored(
      L:max_delay, plnorm,
      pwindow = pwindow, swindow = 1, D = D, L = L,
      meanlog = params[1], sdlog = params[2], log = TRUE
    )

    expect_equal(stan_lpmf[valid_indices], r_lpmf, tolerance = 1e-6)
  }
)

test_that(
  "Stan primarycensored functions handle L > 0 with D = Inf",
  {
    d_values <- seq(0, 10, by = 0.5)
    D <- Inf
    L <- 2

    dist_id <- 1 # Lognormal
    params <- c(1, 0.5) # meanlog, sdlog
    pwindow <- 1
    primary_id <- 1 # Uniform
    primary_params <- array(numeric(0))

    stan_cdf <- sapply(
      d_values, primarycensored_cdf, dist_id, params, pwindow, L, D,
      primary_id, primary_params
    )
    r_cdf <- pprimarycensored(
      d_values, plnorm,
      pwindow = pwindow, D = D, L = L, meanlog = params[1], sdlog = params[2]
    )

    expect_equal(stan_cdf, r_cdf, tolerance = 1e-6)
  }
)

test_that("Stan primarycensored_cdf matches R pprimarycensored", {
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
      primary_id <- 1 # Uniform
      primary_params <- array(numeric(0))

      stan_cdf <- sapply(
        d, primarycensored_cdf, dist_id, params, pwindow, 0, D,
        primary_id, primary_params
      )
      r_cdf <- pprimarycensored(
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
  "Stan primarycensored_lcdf matches R pprimarycensored with
   log.p = TRUE",
  {
    d <- seq(0, 10, by = 0.5)
    dist_id <- 1 # Lognormal
    params <- c(0, 1) # meanlog, sdlog
    pwindow <- 1
    D <- 12
    primary_id <- 1 # Uniform
    primary_params <- numeric(0)

    stan_lcdf <- sapply(
      d, primarycensored_lcdf, dist_id, params, pwindow, 0, D,
      primary_id, primary_params
    )
    r_cdf <- pprimarycensored(
      d, plnorm,
      pwindow = pwindow, D = D, meanlog = params[1],
      sdlog = params[2]
    )

    expect_equal(exp(stan_lcdf), r_cdf, tolerance = 1e-6)

    # Gamma
    dist_id <- 2 # Gamma
    params <- c(2, 1) # shape, scale
    pwindow <- 2
    D <- 12
    primary_id <- 1 # Uniform
    primary_params <- numeric(0)

    stan_lcdf <- sapply(
      d, primarycensored_lcdf, dist_id, params, pwindow, 0, D,
      primary_id, primary_params
    )
    r_cdf <- pprimarycensored(
      d, pgamma,
      pwindow = pwindow, D = D, shape = params[1], scale = params[2]
    )

    expect_equal(exp(stan_lcdf), r_cdf, tolerance = 1e-6)
  }
)

test_that(
  "Stan primarycensored_lpmf throws an error for invalid upper truncation
   point",
  {
    d <- 10
    dist_id <- 1 # Lognormal
    params <- c(0, 1) # meanlog, sdlog
    pwindow <- 1
    d_upper <- 11
    D <- 10
    primary_id <- 1 # Uniform
    primary_params <- numeric(0)

    expect_error(
      primarycensored_lpmf(
        d, dist_id, params, pwindow, d_upper, 0, D, primary_id, primary_params
      ),
      "Upper truncation point is greater than D"
    )
  }
)

test_that(
  "Stan primarycensored matches R primarycensored when d is the same
   as D - 1",
  {
    dist_id <- 1 # Lognormal
    params <- c(0, 1) # meanlog, sdlog
    pwindow <- 1
    primary_id <- 1 # Uniform
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
        stan_pmf <- primarycensored_pmf(
          d, dist_id, params, pwindow, d + 1, 0, D, primary_id,
          primary_params
        )
        r_pmf <- dprimarycensored(
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

test_that("Stan primarycensored_pmf matches R dprimarycensored", {
  d <- 0:10
  dist_id <- 1 # Lognormal
  params <- c(0, 1) # meanlog, sdlog
  pwindow <- 1
  d_upper <- d + 1
  D <- Inf
  primary_id <- 1 # Uniform
  primary_params <- numeric(0)

  stan_pmf <- mapply(
    primarycensored_pmf,
    d = d,
    d_upper = d + 1,
    MoreArgs = list(
      dist_id = dist_id,
      params = params,
      pwindow = pwindow,
      L = 0,
      D = D,
      primary_id = primary_id,
      primary_params = primary_params
    )
  )
  r_pmf <- dprimarycensored(
    d, plnorm,
    pwindow = pwindow, swindow = 1, D = D,
    meanlog = params[1], sdlog = params[2]
  )

  expect_equal(stan_pmf, r_pmf, tolerance = 1e-6)
})

test_that(
  "Stan primarycensored_lpmf matches R dprimarycensored with
   log = TRUE",
  {
    d <- 0:10
    dist_id <- 1 # Lognormal
    params <- c(0, 1) # meanlog, sdlog
    pwindow <- 1
    d_upper <- d + 1
    D <- Inf
    primary_id <- 1 # Uniform
    primary_params <- numeric(0)

    stan_lpmf <- mapply(
      primarycensored_lpmf,
      d = d,
      d_upper = d_upper,
      MoreArgs = list(
        dist_id = dist_id,
        params = params,
        pwindow = pwindow,
        L = 0,
        D = D,
        primary_id = primary_id,
        primary_params = primary_params
      )
    )
    r_lpmf <- log(
      dprimarycensored(
        d, plnorm,
        pwindow = pwindow, swindow = 1, D = D,
        meanlog = params[1], sdlog = params[2]
      )
    )

    expect_equal(stan_lpmf, r_lpmf, tolerance = 1e-6)
  }
)

test_that(
  "Stan primarycensored_sone_pmf_vectorized matches R dprimarycensored",
  {
    max_delay <- 10
    dist_id <- 1 # Lognormal
    params <- c(0, 1) # meanlog, sdlog
    pwindow <- 1
    D <- Inf
    primary_id <- 1 # Uniform
    primary_params <- numeric(0)

    stan_pmf <- primarycensored_sone_pmf_vectorized(
      max_delay, 0, D, dist_id, params, pwindow, primary_id, primary_params
    )
    r_pmf <- dprimarycensored(
      0:max_delay, plnorm,
      pwindow = pwindow, swindow = 1, D = D,
      meanlog = params[1], sdlog = params[2]
    )

    expect_equal(stan_pmf, r_pmf, tolerance = 1e-6)
  }
)

test_that(
  "Stan primarycensored_sone_pmf_vectorized matches R dprimarycensored
   with finite D",
  {
    max_delay <- 10
    dist_id <- 1 # Lognormal
    params <- c(0, 1) # meanlog, sdlog
    pwindow <- 1
    D <- 15
    primary_id <- 1 # Uniform
    primary_params <- numeric(0)

    stan_pmf <- primarycensored_sone_pmf_vectorized(
      max_delay, 0, D, dist_id, params, pwindow, primary_id, primary_params
    )
    r_pmf <- dprimarycensored(
      0:max_delay, plnorm,
      pwindow = pwindow, swindow = 1, D = D,
      meanlog = params[1], sdlog = params[2]
    )

    expect_equal(stan_pmf, r_pmf, tolerance = 1e-6)
  }
)

test_that(
  "Stan primarycensored_sone_pmf_vectorized matches R dprimarycensored
   with D equal to max_delay + 1",
  {
    max_delay <- 10
    dist_id <- 1 # Lognormal
    params <- c(0, 1) # meanlog, sdlog
    pwindow <- 1
    D <- max_delay + 1
    primary_id <- 1 # Uniform
    primary_params <- numeric(0)

    stan_pmf <- primarycensored_sone_pmf_vectorized(
      max_delay, 0, D, dist_id, params, pwindow, primary_id, primary_params
    )
    r_pmf <- dprimarycensored(
      0:max_delay, plnorm,
      pwindow = pwindow, swindow = 1, D = D,
      meanlog = params[1], sdlog = params[2]
    )

    expect_equal(stan_pmf, r_pmf, tolerance = 1e-6)
  }
)

test_that(
  "Stan primarycensored_sone_lpmf_vectorized matches R dprimarycensored
   with log = TRUE",
  {
    max_delay <- 10
    dist_id <- 1 # Lognormal
    params <- c(0, 1) # meanlog, sdlog
    pwindow <- 1
    D <- Inf
    primary_id <- 1 # Uniform
    primary_params <- numeric(0)

    stan_lpmf <- primarycensored_sone_lpmf_vectorized(
      max_delay, 0, D, dist_id, params, pwindow, primary_id, primary_params
    )
    r_lpmf <- dprimarycensored(
      0:max_delay, plnorm,
      pwindow = pwindow, swindow = 1, D = D,
      meanlog = params[1], sdlog = params[2], log = TRUE
    )

    expect_equal(stan_lpmf, r_lpmf, tolerance = 1e-6)
  }
)
