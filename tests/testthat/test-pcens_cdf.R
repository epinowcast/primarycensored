test_that("pcens_cdf methods dispatch correctly to existing
   analytical solutions", {
  pdist <- pgamma
  dprimary <- dunif

  obj_gamma <- new_pcens(
    pdist,
    dprimary,
    list(),
    shape = 2,
    rate = 1
  )

  pdist <- plnorm
  dprimary <- dunif

  obj_lnorm <- new_pcens(
    pdist,
    dprimary,
    list(),
    meanlog = 0,
    sdlog = 1
  )

  pdist <- pweibull
  dprimary <- dunif

  obj_weibull <- new_pcens(
    pdist,
    dprimary,
    list(),
    shape = 2,
    scale = 1
  )

  expect_s3_class(obj_gamma, "pcens_pgamma_dunif")
  expect_s3_class(obj_lnorm, "pcens_plnorm_dunif")
  expect_s3_class(obj_weibull, "pcens_pweibull_dunif")

  q_values <- c(5, 10)
  pwindow <- 2

  expect_no_error(
    pcens_cdf(obj_gamma, q = q_values, pwindow = pwindow)
  )
  expect_no_error(
    pcens_cdf(obj_lnorm, q = q_values, pwindow = pwindow)
  )
  expect_no_error(
    pcens_cdf(obj_weibull, q = q_values, pwindow = pwindow)
  )
})

test_that("pcens_cdf errors as expected when the wrong distributional
   parameters are supplied", {
  pdist <- pgamma
  dprimary <- dunif

  obj_gamma <- new_pcens(
    pdist,
    dprimary,
    list(),
    rate = 1
  )

  expect_error(
    pcens_cdf(obj_gamma, q = 1, pwindow = 1),
    "shape parameter is required for Gamma distribution"
  )

  obj_gamma_no_rate <- new_pcens(
    pdist,
    dprimary,
    list(),
    shape = 2
  )

  expect_error(
    pcens_cdf(obj_gamma_no_rate, q = 1, pwindow = 1),
    "scale or rate parameter is required for Gamma distribution"
  )

  pdist <- plnorm

  obj_lnorm_no_meanlog <- new_pcens(
    pdist,
    dprimary,
    list(),
    sdlog = 1
  )

  expect_error(
    pcens_cdf(obj_lnorm_no_meanlog, q = 1, pwindow = 1),
    "meanlog parameter is required for Log-Normal distribution"
  )

  obj_lnorm_no_sdlog <- new_pcens(
    pdist,
    dprimary,
    list(),
    meanlog = 0
  )

  expect_error(
    pcens_cdf(obj_lnorm_no_sdlog, q = 1, pwindow = 1),
    "sdlog parameter is required for Log-Normal distribution"
  )

  pdist <- pweibull

  obj_weibull_no_shape <- new_pcens(
    pdist,
    dprimary,
    list(),
    scale = 1
  )

  expect_error(
    pcens_cdf(obj_weibull_no_shape, q = 1, pwindow = 1),
    "shape parameter is required for Weibull distribution"
  )

  obj_weibull_no_scale <- new_pcens(
    pdist,
    dprimary,
    list(),
    shape = 2
  )

  expect_error(
    pcens_cdf(obj_weibull_no_scale, q = 1, pwindow = 1),
    "scale parameter is required for Weibull distribution"
  )
})

test_that("pcens_cdf.default computes the same values as
   pcens_cdf.pcens_pgamma_dunif", {
  pdist <- pgamma
  dprimary <- dunif

  shapes <- c(0.5, 1, 2, 5)
  rates <- c(0.1, 0.5, 1, 2)
  pwindows <- c(1, 2, 5, 10)

  for (shape in shapes) {
    for (rate in rates) {
      for (pwindow in pwindows) {
        obj <- new_pcens(
          pdist,
          dprimary,
          list(),
          shape = shape,
          rate = rate
        )

        q_values <- seq(0, 30, by = 0.1)
        result_numeric <- pcens_cdf(
          obj,
          q = q_values,
          pwindow = pwindow,
          use_numeric = TRUE
        )
        result_analytical <- pcens_cdf(
          obj,
          q = q_values,
          pwindow = pwindow,
          use_numeric = FALSE
        )

        # Check properties of numeric result
        expect_type(result_numeric, "double")
        expect_length(result_numeric, length(q_values))
        expect_true(
          all(diff(result_numeric) >= 0)
        ) # Ensure CDF is non-decreasing

        # Check that analytical and numeric results are the same
        expect_equal(
          result_numeric,
          result_analytical,
          tolerance = 1e-5,
          info = sprintf(
            "Mismatch for shape = %s, rate = %s, pwindow = %s",
            shape,
            rate,
            pwindow
          )
        )
      }
    }
  }
})

test_that("pcens_cdf.default computes the same values as
   pcens_cdf.pcens_plnorm_dunif", {
  pdist <- plnorm
  dprimary <- dunif

  meanlogs <- c(-1, 0, 1, 2)
  sdlogs <- c(0.5, 1, 1.5)
  pwindows <- c(1, 2, 5, 8)

  for (meanlog in meanlogs) {
    for (sdlog in sdlogs) {
      for (pwindow in pwindows) {
        obj <- new_pcens(
          pdist,
          dprimary,
          list(),
          meanlog = meanlog,
          sdlog = sdlog
        )

        q_values <- seq(0, 30, by = 0.1)
        result_numeric <- pcens_cdf(
          obj,
          q = q_values,
          pwindow = pwindow,
          use_numeric = TRUE
        )
        result_analytical <- pcens_cdf(
          obj,
          q = q_values,
          pwindow = pwindow,
          use_numeric = FALSE
        )
        # Check properties of numeric result
        expect_type(result_numeric, "double")
        expect_length(result_numeric, length(q_values))
        expect_true(
          all(diff(result_numeric) >= 0)
        ) # Ensure CDF is non-decreasing

        # Check that analytical and numeric results are the same
        expect_equal(
          result_numeric,
          result_analytical,
          tolerance = 1e-5,
          info = sprintf(
            "Mismatch for meanlog = %s, sdlog = %s, pwindow = %s",
            meanlog,
            sdlog,
            pwindow
          )
        )
      }
    }
  }
})

test_that("pcens_cdf.default computes the same values as
   pcens_cdf.pcens_pweibull_dunif", {
  pdist <- pweibull
  dprimary <- dunif

  shapes <- c(0.5, 1, 2)
  scales <- c(0.5, 1, 2)
  pwindows <- c(1, 2, 3, 4, 5)

  for (shape in shapes) {
    for (scale in scales) {
      for (pwindow in pwindows) {
        obj <- new_pcens(
          pdist,
          dprimary,
          list(),
          shape = shape,
          scale = scale
        )

        q_values <- seq(0, 30, by = 0.1)
        result_numeric <- pcens_cdf(
          obj,
          q = q_values,
          pwindow = pwindow,
          use_numeric = TRUE
        )
        result_analytical <- pcens_cdf(
          obj,
          q = q_values,
          pwindow = pwindow,
          use_numeric = FALSE
        )

        # Check properties of numeric result
        expect_type(result_numeric, "double")
        expect_length(result_numeric, length(q_values))
        expect_true(
          all(diff(result_numeric) >= 0)
        ) # Ensure CDF is non-decreasing

        # Check that analytical and numeric results are the same
        expect_equal(
          result_numeric,
          result_analytical,
          tolerance = 1e-5,
          info = sprintf(
            "Mismatch for shape = %s, scale = %s, pwindow = %s",
            shape,
            scale,
            pwindow
          )
        )
      }
    }
  }
})

test_that("pcens_cdf returns values in [0, 1] for all analytical methods", {
  # Test case from issue #238: parameters that stress numerical precision
  # Lognormal with tight distribution where CDF approaches 1 quickly
  obj_lnorm <- new_pcens(
    plnorm,
    dunif,
    list(),
    meanlog = 0.55,
    sdlog = 0.27
  )

  q_values <- seq(0, 30, by = 1)
  cdf_lnorm <- pcens_cdf(obj_lnorm, q = q_values, pwindow = 1)
  expect_true(
    all(cdf_lnorm >= 0 & cdf_lnorm <= 1),
    info = "Lognormal CDF should be in [0, 1]"
  )

  # Gamma distribution
  obj_gamma <- new_pcens(
    pgamma,
    dunif,
    list(),
    shape = 2,
    rate = 0.5
  )

  cdf_gamma <- pcens_cdf(obj_gamma, q = q_values, pwindow = 1)
  expect_true(
    all(cdf_gamma >= 0 & cdf_gamma <= 1),
    info = "Gamma CDF should be in [0, 1]"
  )

  # Weibull distribution
  obj_weibull <- new_pcens(
    pweibull,
    dunif,
    list(),
    shape = 2,
    scale = 1
  )

  cdf_weibull <- pcens_cdf(obj_weibull, q = q_values, pwindow = 1)
  expect_true(
    all(cdf_weibull >= 0 & cdf_weibull <= 1),
    info = "Weibull CDF should be in [0, 1]"
  )

  # Test with exponential growth primary (uses numerical integration)
  obj_expgrowth <- new_pcens(
    plnorm,
    dexpgrowth,
    list(r = 0.2),
    meanlog = 0.55,
    sdlog = 0.27
  )

  cdf_expgrowth <- pcens_cdf(obj_expgrowth, q = q_values, pwindow = 1)
  expect_true(
    all(cdf_expgrowth >= 0 & cdf_expgrowth <= 1),
    info = "CDF with expgrowth primary should be in [0, 1]"
  )
})

test_that("new_pcens works with custom function with name attribute", {
  # Create custom functions with name attributes
  custom_pdist <- function(x, shape, rate) pgamma(x, shape, rate)
  custom_dprimary <- function(x, min = 0, max = 1) dunif(x, min, max)

  named_pdist <- add_name_attribute(custom_pdist, "pgamma")
  named_dprimary <- add_name_attribute(custom_dprimary, "dunif")

  # Create pcens object with custom named functions
  obj <- new_pcens(
    named_pdist,
    named_dprimary,
    list(),
    shape = 2,
    rate = 1
  )

  # Check class is set correctly using function names
  expect_s3_class(obj, "pcens_pgamma_dunif")

  # Check functions are preserved correctly
  expect_identical(body(obj$pdist), body(custom_pdist))
  expect_identical(formals(obj$pdist), formals(custom_pdist))
  expect_identical(body(obj$dprimary), body(custom_dprimary))
  expect_identical(formals(obj$dprimary), formals(custom_dprimary))

  # Check arguments are preserved
  expect_identical(obj$args, list(shape = 2, rate = 1))
})

test_that("pcens_cdf.pcens_pdiscretestep_dunif dispatches and returns [0,1]", {
  boundaries <- 0:3
  pmf <- c(0.2, 0.5, 0.3)
  obj <- new_pcens(
    pdiscretestep,
    dunif,
    list(),
    boundaries = boundaries,
    pmf = pmf
  )
  expect_s3_class(obj, "pcens_pdiscretestep_dunif")
  expect_s3_class(obj, "pcens_pdiscretestep")
  expect_s3_class(obj, "pcens")
  q_values <- seq(0, 4, by = 0.5)
  result <- pcens_cdf(obj, q = q_values, pwindow = 1)
  expect_true(all(result >= 0 & result <= 1))
  expect_true(all(diff(result) >= 0))
})

test_that(
  "pcens_cdf.pcens_pdiscretestep_dunif analytic matches numeric integration",
  {
    boundaries <- 0:3
    pmf <- c(0.2, 0.5, 0.3)
    obj <- new_pcens(
      pdiscretestep,
      dunif,
      list(),
      boundaries = boundaries,
      pmf = pmf
    )
    pwindow <- 1
    # Compute analytic result at points not on boundaries
    q_values <- c(0.3, 0.7, 1.2, 1.7, 2.2, 2.7)
    analytic <- pcens_cdf(obj, q = q_values, pwindow = pwindow)
    # Compute numeric reference by hand:
    # F_obs(q) = (1/pwindow) * integrate F_step(q-p) dp from 0 to pwindow
    numeric_ref <- vapply(q_values, function(qi) {
      stats::integrate(
        function(p) pdiscretestep(qi - p, boundaries, pmf),
        lower = 0,
        upper = pwindow
      )$value
    }, numeric(1))
    expect_equal(analytic, numeric_ref, tolerance = 1e-8)
  }
)

test_that(
  "pcens_cdf.pcens_pdiscretestep_dunif errors for step width > pwindow",
  {
    # boundaries 0, 3 — step width 3 > pwindow 1
    boundaries <- c(0, 3, 4)
    pmf <- c(0.6, 0.4)
    obj <- new_pcens(
      pdiscretestep,
      dunif,
      list(),
      boundaries = boundaries,
      pmf = pmf
    )
    expect_error(
      pcens_cdf(obj, q = 2, pwindow = 1),
      "pwindow"
    )
  }
)

test_that("pcens_cdf.pcens_pdiscretestep_dunif handles two-bin PMF correctly", {
  # Simplest case: two bins, pwindow = 1, hand-computed
  # boundaries = 0:2, pmf = c(0.4, 0.6)
  # right edges: 1, 2
  # F_step: 0 for q < 1, 0.4 for 1 <= q < 2, 1 for q >= 2
  # F_obs(q) = (1/1) * integral_0^1 F_step(q-p) dp
  # For q = 1.5: integral_0^1 F_step(1.5-p) dp
  #   F_step(1.5-p): p in [0, 0.5] -> 1.5-p in [1,1.5] -> F=0.4
  #                  p in (0.5, 1] -> 1.5-p in [0.5,1) -> F=0
  #   = 0.4 * 0.5 + 0 * 0.5 = 0.2
  boundaries <- 0:2
  pmf <- c(0.4, 0.6)
  obj <- new_pcens(
    pdiscretestep,
    dunif,
    list(),
    boundaries = boundaries,
    pmf = pmf
  )
  result <- pcens_cdf(obj, q = 1.5, pwindow = 1)
  expect_equal(result, 0.2, tolerance = 1e-10)
})

test_that(
  "pcens_cdf.pcens_pdiscretestep analytic matches numeric for expgrowth",
  {
    boundaries <- 0:3
    pmf <- c(0.2, 0.5, 0.3)
    pwindow <- 1
    r <- 0.5
    obj <- new_pcens(
      pdiscretestep,
      dexpgrowth,
      list(r = r),
      boundaries = boundaries,
      pmf = pmf
    )
    # Verify pprimary was found in registry
    expect_false(is.null(obj$pprimary))

    q_values <- c(0.3, 0.7, 1.2, 1.7, 2.2, 2.7)
    analytic <- pcens_cdf(obj, q = q_values, pwindow = pwindow)

    # Reference: partition integral using pexpgrowth directly
    ref <- vapply(q_values, function(qi) {
      K <- length(pmf)
      cum_pmf <- cumsum(pmf)
      right_edges <- boundaries[-1L]
      .fstep <- function(x) {
        idx <- findInterval(x, right_edges, left.open = FALSE)
        if (idx == 0L) 0 else if (idx >= K) 1 else cum_pmf[idx]
      }
      p_knots <- qi - right_edges
      inside <- p_knots[p_knots > 0 & p_knots < pwindow]
      breaks <- sort(unique(c(0, inside, pwindow)))
      total <- 0
      for (j in seq_len(length(breaks) - 1L)) {
        a <- breaks[j]
        b <- breaks[j + 1L]
        if ((b - a) <= 0) next
        mid <- 0.5 * (a + b)
        c_k <- .fstep(qi - mid)
        total <- total + c_k * (
          pexpgrowth(b, min = 0, max = pwindow, r = r) -
            pexpgrowth(a, min = 0, max = pwindow, r = r)
        )
      }
      total
    }, numeric(1))

    expect_equal(analytic, ref, tolerance = 1e-8)
  }
)

test_that(
  "pcens_cdf.pcens_pdiscretehazard dispatch chain is correct",
  {
    hazards <- c(0.2, 0.3, 1)
    boundaries <- 0:3
    pmf_equiv <- hazards_to_pmf(hazards)
    pwindow <- 1

    obj_haz <- new_pcens(
      pdiscretehazard,
      dunif,
      list(),
      boundaries = boundaries,
      hazards = hazards
    )
    expect_s3_class(obj_haz, "pcens_pdiscretehazard_dunif")
    expect_s3_class(obj_haz, "pcens_pdiscretehazard")
    expect_s3_class(obj_haz, "pcens")

    obj_step <- new_pcens(
      pdiscretestep,
      dunif,
      list(),
      boundaries = boundaries,
      pmf = pmf_equiv
    )

    q_values <- c(0.5, 1, 1.5, 2, 2.5)
    result_haz <- pcens_cdf(obj_haz, q = q_values, pwindow = pwindow)
    result_step <- pcens_cdf(obj_step, q = q_values, pwindow = pwindow)

    expect_equal(result_haz, result_step, tolerance = 1e-12)
  }
)
