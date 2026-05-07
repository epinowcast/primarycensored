# fmt: skip file
# Skip tests if fitdistrplus is not installed
if (!requireNamespace("fitdistrplus", quietly = TRUE)) {
  skip("Package 'fitdistrplus' is required for these tests")
}

test_that("fitdistdoublecens works correctly with column names", {
  # Set seed for reproducibility
  set.seed(123)

  # Define true distribution parameters
  n <- 1000
  shape <- 1.77
  rate <- 0.44

  # Generate samples
  samples <- rprimarycensored(
    n,
    rgamma,
    shape = shape,
    rate = rate,
    pwindow = 1,
    swindow = 1,
    D = 8
  )

  # Create data frame with column names
  delay_data <- data.frame(
    left = samples,
    right = samples + 1,
    pwindow = rep(1, n),
    D = rep(8, n)
  )

  # Fit the model using fitdistdoublecens
  fit <- fitdistdoublecens(
    delay_data,
    distr = "gamma",
    start = list(shape = 1, rate = 1)
  )

  # Check that the function returns a fitdist object
  expect_s3_class(fit, "fitdist")

  # Check that the estimated parameters are close to the true values
  expect_equal(unname(fit$estimate["shape"]), shape, tolerance = 0.2)
  expect_equal(unname(fit$estimate["rate"]), rate, tolerance = 0.2)

  # Check that the log-likelihood is not NA or -Inf
  expect_false(is.na(fit$loglik))
  expect_false(is.infinite(fit$loglik))

  # Check that the AIC and BIC are calculated
  expect_false(is.na(fit$aic))
  expect_false(is.na(fit$bic))
})

test_that("fitdistdoublecens handles errors correctly", {
  # Test with missing columns
  expect_error(
    fitdistdoublecens(
      data.frame(x = 1:10), # Missing required columns
      distr = "gamma"
    ),
    "Missing required columns"
  )

  # Test with non-existent distribution
  expect_error(
    fitdistdoublecens(
      data.frame(
        left = 1:10,
        right = 2:11,
        pwindow = rep(1, 10),
        D = rep(10, 10)
      ),
      distr = "nonexistent_dist"
    )
  )

  # Test with observations below L
  expect_error(
    fitdistdoublecens(
      data.frame(
        left = 0:9,
        right = 1:10,
        pwindow = rep(1, 10),
        L = rep(2, 10),
        D = rep(15, 10)
      ),
      distr = "gamma",
      start = list(shape = 1, rate = 1)
    ),
    "Observations must be >= L"
  )

  # Test with L >= D
  expect_error(
    fitdistdoublecens(
      data.frame(
        left = 1:10,
        right = 2:11,
        pwindow = rep(1, 10),
        L = rep(10, 10),
        D = rep(10, 10)
      ),
      distr = "gamma",
      start = list(shape = 1, rate = 1)
    ),
    "L must be less than D"
  )
})

test_that("fitdistdoublecens works with different distributions", {
  set.seed(123)
  n <- 1000

  # Test with normal distribution. `fitdistdoublecens()` treats a missing
  # `L` column as `L = -Inf`, so the sampler uses the default upper-only
  # truncation as well.
  true_mean <- 5
  true_sd <- 2
  samples <- rprimarycensored(
    n,
    rnorm,
    mean = true_mean,
    sd = true_sd,
    pwindow = 2,
    swindow = 2,
    D = 10
  )

  delay_data <- data.frame(
    left = samples,
    right = samples + 2,
    pwindow = rep(2, n),
    D = rep(10, n)
  )

  fit_norm <- fitdistdoublecens(
    delay_data,
    distr = "norm",
    start = list(mean = 0, sd = 1)
  )

  expect_s3_class(fit_norm, "fitdist")
  expect_equal(unname(fit_norm$estimate["mean"]), true_mean, tolerance = 0.2)
  expect_equal(unname(fit_norm$estimate["sd"]), true_sd, tolerance = 0.2)
})

test_that("fitdistdoublecens works with mixed secondary windows", {
  set.seed(456)
  n <- 1000

  # True parameters for gamma distribution
  true_shape <- 3
  true_rate <- 0.5

  # Generate samples with mixed secondary windows
  generate_sample <- function(pwindow, swindow, obs_time) {
    rpcens(
      1,
      rgamma,
      shape = true_shape,
      rate = true_rate,
      pwindow = pwindow,
      swindow = swindow,
      D = obs_time
    )
  }

  pwindows <- rep(1, n)
  swindows <- sample(c(1, 2), n, replace = TRUE)
  obs_times <- rep(10, n)

  samples <- mapply(generate_sample, pwindows, swindows, obs_times)

  delay_data <- data.frame(
    left = samples,
    right = samples + swindows,
    pwindow = pwindows,
    D = obs_times
  )

  fit_gamma <- fitdistdoublecens(
    delay_data,
    distr = "gamma",
    start = list(shape = 2, rate = 1)
  )

  expect_s3_class(fit_gamma, "fitdist")
  expect_equal(unname(fit_gamma$estimate["shape"]), true_shape, tolerance = 0.3)
  expect_equal(unname(fit_gamma$estimate["rate"]), true_rate, tolerance = 0.2)
})

test_that("fitdistdoublecens works with mixed D and primary windows", {
  set.seed(789)
  n <- 1000

  # True parameters for gamma distribution
  true_shape <- 2.5
  true_rate <- 0.6

  # Generate samples with mixed D and primary windows
  generate_sample <- function(pwindow, swindow, obs_time) {
    rpcens(
      1,
      rgamma,
      shape = true_shape,
      rate = true_rate,
      pwindow = pwindow,
      swindow = swindow,
      D = obs_time
    )
  }

  # Create mixed pwindows and D values
  pwindows <- sample(c(1, 2), n, replace = TRUE)
  swindows <- rep(1, n)
  obs_times <- sample(c(8, 12), n, replace = TRUE)

  samples <- mapply(generate_sample, pwindows, swindows, obs_times)

  delay_data <- data.frame(
    left = samples,
    right = samples + swindows,
    pwindow = pwindows,
    D = obs_times
  )

  fit_gamma <- fitdistdoublecens(
    delay_data,
    distr = "gamma",
    start = list(shape = 2, rate = 1)
  )

  expect_s3_class(fit_gamma, "fitdist")
  expect_equal(unname(fit_gamma$estimate["shape"]), true_shape, tolerance = 0.3)
  expect_equal(unname(fit_gamma$estimate["rate"]), true_rate, tolerance = 0.3)
})

test_that("fitdistdoublecens works with custom column names", {
  set.seed(123)
  n <- 1000
  shape <- 1.77
  rate <- 0.44

  samples <- rprimarycensored(
    n,
    rgamma,
    shape = shape,
    rate = rate,
    pwindow = 1,
    swindow = 1,
    D = 8
  )

  # Create data frame with custom column names
  delay_data <- data.frame(
    lower_bound = samples,
    upper_bound = samples + 1,
    primary_window = rep(1, n),
    truncation_time = rep(8, n)
  )

  fit <- fitdistdoublecens(
    delay_data,
    distr = "gamma",
    start = list(shape = 1, rate = 1),
    left = "lower_bound",
    right = "upper_bound",
    pwindow = "primary_window",
    D = "truncation_time"
  )

  expect_s3_class(fit, "fitdist")
  expect_equal(unname(fit$estimate["shape"]), shape, tolerance = 0.2)
  expect_equal(unname(fit$estimate["rate"]), rate, tolerance = 0.2)
})

test_that("fitdistdoublecens handles additional columns correctly", {
  set.seed(123)
  n <- 1000
  shape <- 1.77
  rate <- 0.44

  samples <- rprimarycensored(
    n,
    rgamma,
    shape = shape,
    rate = rate,
    pwindow = 1,
    swindow = 1,
    D = 8
  )

  # Create data frame with extra columns that should be ignored
  delay_data <- data.frame(
    left = samples,
    right = samples + 1,
    pwindow = rep(1, n),
    D = rep(8, n),
    extra_column = rep("extra", n),
    another_extra = seq_len(n),
    metadata = rep(TRUE, n),
    stringsAsFactors = FALSE
  )

  # Should work without error despite extra columns
  fit <- fitdistdoublecens(
    delay_data,
    distr = "gamma",
    start = list(shape = 1, rate = 1)
  )

  expect_s3_class(fit, "fitdist")
  expect_equal(unname(fit$estimate["shape"]), shape, tolerance = 0.2)
  expect_equal(unname(fit$estimate["rate"]), rate, tolerance = 0.2)
})

test_that("fitdistdoublecens handles additional columns with custom names", {
  set.seed(123)
  n <- 1000
  shape <- 1.77
  rate <- 0.44

  samples <- rprimarycensored(
    n,
    rgamma,
    shape = shape,
    rate = rate,
    pwindow = 1,
    swindow = 1,
    D = 8
  )

  # Create data frame with custom column names AND extra columns
  delay_data <- data.frame(
    lower_bound = samples,
    upper_bound = samples + 1,
    primary_window = rep(1, n),
    truncation_time = rep(8, n),
    spurious_data = runif(n),
    id = seq_len(n),
    notes = rep("test", n),
    stringsAsFactors = FALSE
  )

  fit <- fitdistdoublecens(
    delay_data,
    distr = "gamma",
    start = list(shape = 1, rate = 1),
    left = "lower_bound",
    right = "upper_bound",
    pwindow = "primary_window",
    D = "truncation_time"
  )

  expect_s3_class(fit, "fitdist")
  expect_equal(unname(fit$estimate["shape"]), shape, tolerance = 0.2)
  expect_equal(unname(fit$estimate["rate"]), rate, tolerance = 0.2)
})

test_that("fitdistdoublecens handles truncation_check_multiplier correctly", {
  set.seed(123)
  n <- 100
  shape <- 1.77
  rate <- 0.44

  samples <- rprimarycensored(
    n,
    rgamma,
    shape = shape,
    rate = rate,
    pwindow = 1,
    swindow = 1,
    D = 100 # Very large D
  )

  delay_data <- data.frame(
    left = samples,
    right = samples + 1,
    pwindow = rep(1, n),
    D = rep(100, n)
  )

  # Should show a message about large D
  expect_message(
    fitdistdoublecens(
      delay_data,
      distr = "gamma",
      start = list(shape = 1, rate = 1),
      truncation_check_multiplier = 2
    ),
    "truncation time D"
  )

  # Should not show a message when check is disabled
  expect_no_message(
    fitdistdoublecens(
      delay_data,
      distr = "gamma",
      start = list(shape = 1, rate = 1),
      truncation_check_multiplier = NULL
    )
  )
})

test_that(
  "fitdistdoublecens throws error when required packages are not installed",
  {
    # Create dummy data
    dummy_data <- data.frame(
      left = 1:5,
      right = 2:6,
      pwindow = rep(1, 5),
      D = rep(10, 5)
    )

    # Test for fitdistrplus
    with_mocked_bindings(
      expect_error(
        fitdistdoublecens(dummy_data, "norm"),
        "Package 'fitdistrplus' is required but not installed for this",
        fixed = TRUE
      ),
      requireNamespace = function(pkg, ...) {
        if (pkg == "fitdistrplus") {
          return(FALSE)
        }
        TRUE
      },
      .package = "base"
    )

    # Test for withr
    with_mocked_bindings(
      expect_error(
        fitdistdoublecens(dummy_data, "norm"),
        "Package 'withr' is required but not installed for this function.",
        fixed = TRUE
      ),
      requireNamespace = function(pkg, ...) {
        if (pkg == "withr") {
          return(FALSE)
        }
        TRUE
      },
      .package = "base"
    )

    # Test when both packages are missing
    with_mocked_bindings(
      expect_error(
        fitdistdoublecens(dummy_data, "norm"),
        "Package 'fitdistrplus' is required but not installed",
        fixed = TRUE
      ),
      requireNamespace = function(...) FALSE,
      .package = "base"
    )
  }
)

# ---- non-parametric tests --------------------------------------------------

test_that("fitdistdoublecens discretestep: recovers approximate PMF", {
  skip_if_not(
    exists("pdiscretestep"),
    message = "pdiscretestep not yet available"
  )
  set.seed(101)
  true_pmf <- c(0.1, 0.3, 0.4, 0.15, 0.05)
  K <- 5L
  # D = K + 1 so the last bin (right edge = 5) is observable with swindow = 1
  D_val <- K + 1L
  n <- 800

  samples <- rprimarycensored(
    n, rdiscretestep,
    boundaries = 0:K, pmf = true_pmf,
    pwindow = 1, swindow = 1, D = D_val
  )

  step_data <- data.frame(
    left = samples,
    right = samples + 1,
    pwindow = rep(1, n),
    D = rep(D_val, n)
  )

  start <- as.list(stats::setNames(rep(1 / K, K - 1L), paste0("p", 1:(K - 1))))
  fit <- fitdistdoublecens(
    step_data,
    distr = "discretestep",
    start = start,
    boundaries = 0:K,
    truncation_check_multiplier = NULL
  )

  expect_s3_class(fit, "fitdist")
  expect_false(is.na(fit$loglik))
  expect_false(is.infinite(fit$loglik))

  # Reconstruct full PMF from estimated free parameters
  p_free <- unname(fit$estimate)
  p_last <- 1 - sum(p_free)
  est_pmf <- c(p_free, p_last)

  expect_equal(sum(est_pmf), 1, tolerance = 1e-6)
  expect_true(all(est_pmf >= -0.01)) # allow tiny numeric noise
  expect_equal(est_pmf, true_pmf, tolerance = 0.1)
})

test_that("fitdistdoublecens discretestep: infers K from start", {
  skip_if_not(
    exists("pdiscretestep"),
    message = "pdiscretestep not yet available"
  )
  set.seed(202)
  true_pmf <- c(0.2, 0.5, 0.3)
  K <- 3L
  D_val <- K + 1L
  n <- 500

  samples <- rprimarycensored(
    n, rdiscretestep,
    boundaries = 0:K, pmf = true_pmf,
    pwindow = 1, swindow = 1, D = D_val
  )

  step_data <- data.frame(
    left = samples,
    right = samples + 1,
    pwindow = rep(1, n),
    D = rep(D_val, n)
  )

  # K not supplied; inferred from start (length 2 -> K = 3)
  fit <- fitdistdoublecens(
    step_data,
    distr = "discretestep",
    start = list(p1 = 0.33, p2 = 0.33),
    boundaries = 0:K,
    truncation_check_multiplier = NULL
  )

  expect_s3_class(fit, "fitdist")
  p_free <- unname(fit$estimate)
  p_last <- 1 - sum(p_free)
  est_pmf <- c(p_free, p_last)
  expect_equal(sum(est_pmf), 1, tolerance = 1e-6)
})

test_that(
  "fitdistdoublecens discretehazard: produces valid PMF close to truth",
  {
    skip_if_not(
      exists("pdiscretestep"),
      message = "pdiscretestep not yet available"
    )
    set.seed(303)
    true_pmf <- c(0.1, 0.3, 0.4, 0.15, 0.05)
    K <- 5L
    D_val <- K + 1L
    n <- 800

    samples <- rprimarycensored(
      n, rdiscretestep,
      boundaries = 0:K, pmf = true_pmf,
      pwindow = 1, swindow = 1, D = D_val
    )

    haz_data <- data.frame(
      left = samples,
      right = samples + 1,
      pwindow = rep(1, n),
      D = rep(D_val, n)
    )

    haz_start <- c(
      list(alpha = -2, log_sigma = log(0.5)),
      as.list(stats::setNames(rep(0, K - 1L), paste0("eps_", seq_len(K - 1L))))
    )
    fit <- fitdistdoublecens(
      haz_data,
      distr = "discretehazard",
      start = haz_start,
      boundaries = 0:K,
      truncation_check_multiplier = NULL
    )

    expect_s3_class(fit, "fitdist")
    expect_false(is.na(fit$loglik))
    expect_false(is.infinite(fit$loglik))

    # Reconstruct PMF from estimated parameters
    est <- fit$estimate
    alpha_val <- est[["alpha"]]
    log_sigma_val <- est[["log_sigma"]]
    eps_vals <- unname(est[grep("^eps_", names(est))])
    sigma <- exp(log_sigma_val)
    logit_h <- alpha_val + sigma * cumsum(c(0, eps_vals))
    h <- 1 / (1 + exp(-logit_h))
    h[K] <- 1
    est_pmf <- hazards_to_pmf(h)

    expect_equal(sum(est_pmf), 1, tolerance = 1e-6)
    expect_true(all(est_pmf >= 0))
    expect_equal(est_pmf, true_pmf, tolerance = 0.1)
  }
)

test_that(
  "fitdistdoublecens discretestep: error when start is missing",
  {
    skip_if_not(
      exists("pdiscretestep"),
      message = "pdiscretestep not yet available"
    )
    dummy <- data.frame(
      left = 1:5,
      right = 2:6,
      pwindow = rep(1, 5),
      D = rep(6, 5)
    )
    expect_error(
      fitdistdoublecens(
        dummy,
        distr = "discretestep",
        truncation_check_multiplier = NULL
      ),
      "`start` must be supplied"
    )
  }
)

test_that("fitdistdoublecens discretehazard: error when start is missing", {
  skip_if_not(
    exists("pdiscretestep"),
    message = "pdiscretestep not yet available"
  )
  dummy <- data.frame(
    left = 1:5,
    right = 2:6,
    pwindow = rep(1, 5),
    D = rep(6, 5)
  )
  expect_error(
    fitdistdoublecens(dummy, distr = "discretehazard"),
    "`start` must be supplied"
  )
})

test_that(
  "fitdistdoublecens discretestep: structural zeros in PMF",
  {
    skip_if_not(
      exists("pdiscretestep"),
      message = "pdiscretestep not yet available"
    )
    set.seed(404)
    # PMF with a zero bin
    true_pmf <- c(0.0, 0.5, 0.5)
    K <- 3L
    D_val <- K + 1L
    n <- 500

    samples <- rprimarycensored(
      n, rdiscretestep,
      boundaries = 0:K, pmf = true_pmf,
      pwindow = 1, swindow = 1, D = D_val
    )

    step_data <- data.frame(
      left = samples,
      right = samples + 1,
      pwindow = rep(1, n),
      D = rep(D_val, n)
    )

    par_names <- paste0("p", 1:(K - 1))
    start <- as.list(stats::setNames(rep(1 / K, K - 1L), par_names))
    # Should complete without error
    fit <- fitdistdoublecens(
      step_data,
      distr = "discretestep",
      start = start,
      boundaries = 0:K,
      truncation_check_multiplier = NULL
    )
    expect_s3_class(fit, "fitdist")
    p_free <- unname(fit$estimate)
    p_last <- 1 - sum(p_free)
    est_pmf <- c(p_free, p_last)
    expect_equal(sum(est_pmf), 1, tolerance = 1e-6)
  }
)

test_that("fitdistdoublecens discretehazard: prior argument changes penalty", {
  skip_if_not(
    exists("pdiscretestep"),
    message = "pdiscretestep not yet available"
  )
  set.seed(606)
  true_pmf <- c(0.2, 0.5, 0.3)
  K <- 3L
  D_val <- K + 1L
  n <- 200
  samples <- rprimarycensored(
    n, rdiscretestep,
    boundaries = 0:K, pmf = true_pmf,
    pwindow = 1, swindow = 1, D = D_val
  )
  haz_data <- data.frame(
    left = samples,
    right = samples + 1,
    pwindow = rep(1, n),
    D = rep(D_val, n)
  )
  haz_start <- c(
    list(alpha = -2, log_sigma = log(0.5)),
    as.list(stats::setNames(rep(0, K - 1L), paste0("eps_", seq_len(K - 1L))))
  )
  fit_default <- suppressWarnings(fitdistdoublecens(
    haz_data,
    distr = "discretehazard",
    start = haz_start,
    boundaries = 0:K,
    truncation_check_multiplier = NULL
  ))
  fit_tight <- suppressWarnings(fitdistdoublecens(
    haz_data,
    distr = "discretehazard",
    start = haz_start,
    boundaries = 0:K,
    prior = list(
      alpha = list(mean = 0, sd = 0.1),
      log_sigma = list(mean = 0, sd = 0.05)
    ),
    truncation_check_multiplier = NULL
  ))
  expect_false(
    isTRUE(all.equal(
      unname(fit_default$estimate), unname(fit_tight$estimate)
    ))
  )
})

test_that(
  "fitdistdoublecens accepts deprecated dprimary_args with a warning",
  {
    set.seed(123)
    n <- 200
    samples <- rprimarycensored(
      n, rgamma,
      shape = 2, rate = 1,
      pwindow = 1, swindow = 1, D = 8
    )
    delay_data <- data.frame(
      left = samples,
      right = samples + 1,
      pwindow = rep(1, n),
      D = rep(8, n)
    )
    expect_warning(
      fitdistdoublecens(
        delay_data,
        distr = "gamma",
        start = list(shape = 1, rate = 1),
        dprimary_args = list()
      ),
      "deprecated"
    )
  }
)

test_that(
  "fitdistdoublecens discretestep: recovers PMF under mixed windows",
  {
    skip_if_not(
      exists("pdiscretestep"),
      message = "pdiscretestep not yet available"
    )
    set.seed(404)
    true_pmf <- c(0.05, 0.20, 0.35, 0.20, 0.10, 0.05, 0.03, 0.02)
    K <- 8L
    n <- 1500
    boundaries <- 0:K

    pwindows <- sample(c(1, 2), n, replace = TRUE)
    swindows <- sample(c(1, 2), n, replace = TRUE)
    obs_times <- sample(c(K + 1L, K + 2L), n, replace = TRUE)

    samples <- mapply(
      function(pw, sw, ot) {
        rprimarycensored(
          1, rdiscretestep,
          boundaries = boundaries, pmf = true_pmf,
          pwindow = pw, swindow = sw, D = ot
        )
      },
      pwindows, swindows, obs_times
    )

    step_data <- data.frame(
      left    = samples,
      right   = pmin(samples + swindows, obs_times),
      pwindow = pwindows,
      D       = obs_times
    )

    start <- as.list(stats::setNames(
      rep(1 / K, K - 1L), paste0("p", seq_len(K - 1L))
    ))
    fit <- fitdistdoublecens(
      step_data,
      distr      = "discretestep",
      start      = start,
      boundaries = boundaries
    )

    p_free <- unname(fit$estimate)
    est_pmf <- c(p_free, 1 - sum(p_free))
    expect_equal(sum(est_pmf), 1, tolerance = 1e-6)
    # Mixed windows + 1500 obs: the recovered PMF should be within
    # ~0.05 absolute on every bin
    expect_lt(max(abs(est_pmf - true_pmf)), 0.05)
  }
)
