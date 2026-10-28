# The analytical uniform primary CDFs share `.pcens_cdf_uniform()`.
# These tests check them against numerical integration at a tight tolerance.
# The agreement with Stan is in test-stan-primarycensored_analytical_cdf.R.

test_that("uniform primary analytical CDFs match tight numerical
  integration across parameters, windows and tails", {
  q <- c(1e-8, 1e-3, 0.05, 0.3, 0.75, 1, 1.5, 2.5, 6, 12, 40)
  for (case in unif_cases()) {
    for (args in case$grid) {
      for (pwindow in c(0.5, 1, 3)) {
        obj <- do.call(
          new_pcens, c(list(case$pdist, dunif), args)
        )
        info <- paste(
          case$name, toString(names(args)), toString(unlist(args)),
          "pwindow", pwindow
        )
        expect_close(
          pcens_cdf(obj, q, pwindow),
          do.call(
            unif_reference,
            c(list(case$pdist, q, pwindow), args, lower = case_lower(case))
          ),
          info = info
        )
      }
    }
  }
})

test_that("uniform primary analytical CDFs are 0 at and below the support and
  tend to 1", {
  for (case in unif_cases()) {
    obj <- do.call(
      new_pcens, c(list(case$pdist, dunif), case$grid[[1]])
    )
    expect_identical(
      pcens_cdf(obj, case_lower(case) - c(5, 0.1, 0), 1), c(0, 0, 0),
      info = case$name
    )
    expect_equal(pcens_cdf(obj, 1e4, 1), 1, tolerance = 1e-10)
  }
  if (requireNamespace("flexsurv", quietly = TRUE)) {
    # A small Q gives a large shape, where the terms at q < 0 are undefined
    obj <- new_pcens(
      flexsurv::pgengamma, dunif, mu = 1, sigma = 0.6, Q = 1e-3
    )
    expect_identical(pcens_cdf(obj, c(-2, -1e-12, 0, 0), 1), rep(0, 4))
  }
})

test_that("uniform primary analytical CDFs handle empty input and error for
  missing q and invalid pwindow", {
  for (case in unif_cases()) {
    obj <- do.call(
      new_pcens, c(list(case$pdist, dunif), case$grid[[1]])
    )
    expect_identical(pcens_cdf(obj, numeric(0), 1), numeric(0))
    expect_identical(pcens_cdf(obj, numeric(0), c(1, 2)), numeric(0))
    expect_error(pcens_cdf(obj, c(1, NA), 1), "missing")
    expect_error(pcens_cdf(obj, NaN, 1), "missing")
    expect_error(pcens_cdf(obj, c(1, 2), NA_real_), "pwindow")
    expect_error(pcens_cdf(obj, c(1, 2), c(1, NA)), "pwindow")
    expect_error(pcens_cdf(obj, c(1, 2), -1), "pwindow")
  }
})

test_that("uniform primary analytical CDFs with a vector pwindow match
  element-wise calls across non-positive, narrow and upper tail cases", {
  q <- c(-1, 0.4, 2, 5, 0, 9, 5, 4e3, 1e5, 3, 1e6, 1e9, 60, 25, Inf)
  pwindow <- c(
    1, 2, 0.5, 3, 1, 4, 1e-7, 1, 1e-2, 1e-9, 0.1, 5, 0.01, 1, 1
  )
  for (case in unif_cases()) {
    obj <- do.call(
      new_pcens, c(list(case$pdist, dunif), case$grid[[1]])
    )
    expected <- vapply(
      seq_along(q), function(i) pcens_cdf(obj, q[[i]], pwindow[[i]]),
      numeric(1)
    )
    expect_equal(
      pcens_cdf(obj, q, pwindow), expected,
      tolerance = 1e-14,
      info = case$name
    )
    expect_equal(
      pcens_cdf(obj, 2, c(1, 2)),
      vapply(c(1, 2), function(w) pcens_cdf(obj, 2, w), numeric(1)),
      tolerance = 1e-14
    )
  }
})

test_that(".pcens_cdf_uniform combines the terms at both ends of the window
  for a new delay family", {
  # Exponential delay, with partial expectation (1 - exp(-rt) (1 + rt)) / r
  rate <- 0.4
  terms_fn <- function(t, p) {
    t * -expm1(-p$rate * t) - p$mean * (1 - exp(-p$rate * t) * (1 + p$rate * t))
  }
  upper_fn <- function(t, p) -p$mean * exp(-p$rate * t)
  delay_cdf <- function(t, p) pexp(t, p$rate)
  q <- c(-1, 0, 0.2, 0.9, 1, 3, 10, 40)
  for (pwindow in c(0.5, 1, 4)) {
    expect_close(
      .pcens_cdf_uniform(
        q, pwindow, terms_fn, upper_fn, delay_cdf,
        list(rate = rate, mean = 1 / rate)
      ),
      unif_reference(pexp, q, pwindow, rate = rate)
    )
  }
})

test_that(".uniform_endpoint_diff evaluates each distinct window end once
  when enough lower ends match a delay", {
  counter <- new.env()
  fn <- function(t, p) {
    counter$n <- counter$n + length(t)
    t^2
  }
  count_evaluations <- function(q, pwindow, share = TRUE) {
    counter$n <- 0L
    lo <- pmax.int(q - pwindow, 0)
    expect_identical(
      .uniform_endpoint_diff(fn, q, lo, list(), share), q^2 - lo^2
    )
    counter$n
  }
  distinct_ends <- function(q, pwindow) {
    length(unique(c(q, pmax.int(q - pwindow, 0))))
  }
  q <- 0:80
  expect_identical(count_evaluations(q, 2), distinct_ends(q, 2))
  q_half <- seq(0, 40, by = 0.5)
  expect_identical(count_evaluations(q_half, 1.5), distinct_ends(q_half, 1.5))
  q_mixed <- c(seq(0.37, 9.9, length.out = 30), 0:50)
  expect_identical(count_evaluations(q_mixed, 1), distinct_ends(q_mixed, 1))
  # Not shared without the flag, for short inputs, or when few lower ends
  # match a delay
  expect_identical(count_evaluations(q, 2, share = FALSE), 2L * length(q))
  expect_identical(count_evaluations(0:30, 1), 62L)
  expect_identical(count_evaluations(q, 1.5), 2L * length(q))
  q_none <- seq(1.37, 40, length.out = 60)
  expect_identical(count_evaluations(q_none, 1), 2L * length(q_none))
  q_zeros <- sort(c(0, 0, 0, seq(1.37, 40, length.out = 200)))
  expect_identical(count_evaluations(q_zeros, 1), 2L * length(q_zeros))
})

test_that("uniform primary analytical CDFs that share terms match single
  point evaluation", {
  grids <- list(
    list(q = 0:100, pwindow = 1),
    list(q = 0:100, pwindow = 3),
    list(q = seq(0, 40, by = 0.5), pwindow = 1.5),
    list(q = seq(0, 20, by = 0.25), pwindow = 0.25),
    list(q = c(seq(0.37, 9.9, length.out = 30), 0:50), pwindow = 1),
    list(q = c(0:100, 2e3, 5e3), pwindow = 1)
  )
  for (case in unif_cases()) {
    for (args in case$grid[1:2]) {
      obj <- do.call(new_pcens, c(list(case$pdist, dunif), args))
      for (grid in grids) {
        expect_equal(
          pcens_cdf(obj, grid$q, grid$pwindow),
          vapply(
            grid$q, function(x) pcens_cdf(obj, x, grid$pwindow), numeric(1)
          ),
          tolerance = 1e-14,
          info = paste(case$name, toString(unlist(args)))
        )
      }
    }
  }
})

test_that("pprimarycensored and dprimarycensored agree on a lattice where
  uniform primary terms are shared", {
  x <- 0:80
  for (case in unif_cases()) {
    args <- c(list(pdist = case$pdist, pwindow = 2), case$grid[[3]])
    p <- do.call(pprimarycensored, c(list(q = x), args))
    expect_equal(
      p,
      vapply(
        x, function(i) do.call(pprimarycensored, c(list(q = i), args)),
        numeric(1)
      ),
      tolerance = 1e-14,
      info = case$name
    )
    d <- do.call(dprimarycensored, c(list(x = x, D = 81), args))
    p_lo <- do.call(pprimarycensored, c(list(q = x, L = 0, D = 81), args))
    p_up <- do.call(
      pprimarycensored, c(list(q = x + 1, L = 0, D = 81), args)
    )
    expect_equal(d, p_up - p_lo, tolerance = 1e-12, info = case$name)
  }
})

test_that("uniform primary analytical CDFs evaluate vector delay parameters
  per point", {
  q <- c(0:10, 2e3, 1e7)
  shape <- seq(1, 3, length.out = 13)
  obj <- new_pcens(pgamma, dunif, shape = shape, scale = 1)
  expected <- vapply(seq_along(q), function(i) {
    pcens_cdf(new_pcens(pgamma, dunif, shape = shape[i], scale = 1), q[i], 1)
  }, numeric(1))
  expect_no_warning(pcens_cdf(obj, q, 1))
  expect_equal(pcens_cdf(obj, q, 1), expected, tolerance = 1e-14)
  expect_equal(
    pcens_cdf(obj, 2, 1),
    vapply(shape, function(s) {
      pcens_cdf(new_pcens(pgamma, dunif, shape = s, scale = 1), 2, 1)
    }, numeric(1)),
    tolerance = 1e-14
  )
  expect_warning(pcens_cdf(obj, 1:5, 1), "multiple")
})

test_that("uniform primary analytical CDFs are 1 at Inf and 0 at -Inf", {
  for (case in unif_cases()) {
    for (args in case$grid) {
      obj <- do.call(new_pcens, c(list(case$pdist, dunif), args))
      expect_identical(
        pcens_cdf(obj, c(-Inf, Inf), 1), c(0, 1), info = case$name
      )
      expect_equal(
        pcens_cdf(obj, c(Inf, 2, -Inf, Inf), c(1, 1, 2, 0.5))[c(1, 3, 4)],
        c(1, 0, 1),
        info = case$name
      )
    }
  }
})

test_that("uniform primary analytical CDFs are 1 far into the upper tail,
  where the lower tail form cancels", {
  q <- c(1e6, 1e9, 1e12, 1e14, 1e15, 1e16, 1e17)
  for (case in unif_cases()) {
    for (pwindow in c(0.1, 1, 5)) {
      obj <- do.call(
        new_pcens, c(list(case$pdist, dunif), case$grid[[1]])
      )
      expect_equal(
        pcens_cdf(obj, q, pwindow), rep(1, length(q)),
        tolerance = 1e-9,
        info = paste(case$name, "pwindow", pwindow)
      )
    }
  }
})

test_that("uniform primary analytical CDFs match numerical integration where
  the upper tail form is used and the tail is not negligible", {
  q <- c(30, 60, 120, 240)
  for (case in unif_cases()) {
    for (args in case$grid) {
      obj <- do.call(new_pcens, c(list(case$pdist, dunif), args))
      expect_close(
        pcens_cdf(obj, q, 0.01),
        do.call(
          unif_reference,
          c(list(case$pdist, q, 0.01), args, lower = case_lower(case))
        ),
        info = paste(case$name, toString(unlist(args)))
      )
    }
  }
})

test_that("uniform primary analytical CDFs agree on either side of the
  switch to the narrow window form", {
  d <- 5
  for (case in unif_cases()) {
    for (args in case$grid[1:3]) {
      obj <- do.call(new_pcens, c(list(case$pdist, dunif), args))
      for (ratio in c(0.999, 1.001)) {
        w <- d / (ratio * .narrow_window_ratio)
        narrow <- pcens_cdf(obj, d, w)
        wide <- with_mocked_bindings(
          pcens_cdf(obj, d, w),
          .narrow_window_ratio = Inf
        )
        expect_close(
          narrow, wide, rtol = 1e-8,
          info = paste(case$name, toString(unlist(args)), "ratio", ratio)
        )
      }
    }
  }
})

test_that("uniform primary analytical CDFs match numerical integration
  across the switch to the upper tail form and with narrow windows", {
  q <- c(0.5, 2, 5, 8, 15, 30, 60, 150)
  for (case in unif_cases()) {
    for (args in case$grid[1:3]) {
      for (pwindow in c(1e-3, 0.05, 2)) {
        obj <- do.call(new_pcens, c(list(case$pdist, dunif), args))
        expect_close(
          pcens_cdf(obj, q, pwindow),
          do.call(
            unif_reference,
            c(list(case$pdist, q, pwindow), args, lower = case_lower(case))
          ),
          info = paste(
            case$name, toString(unlist(args)), "pwindow", pwindow
          )
        )
      }
    }
  }
})

test_that("uniform primary analytical CDFs match numerical integration for
  windows many orders of magnitude narrower than the delay", {
  cases <- list(
    list(pgamma, list(shape = 500, scale = 100), 49967, 1e-6),
    list(pgamma, list(shape = 500, scale = 100), c(4.9e4, 5.1e4), 1e-5),
    list(pweibull, list(shape = 0.1, scale = 1), 4.29e6, 1e-6),
    list(pweibull, list(shape = 0.1, scale = 100), 4.29e8, 1e-3),
    list(pweibull, list(shape = 0.1, scale = 100), 4.29e8, 1e-6),
    list(pweibull, list(shape = 1.5, scale = 5), c(3, 8, 20), 1e-8),
    list(plnorm, list(meanlog = 2, sdlog = 2), c(5, 4e3, 1e5), 1e-7),
    list(plnorm, list(meanlog = 0, sdlog = 1), c(1, 3, 40), 1e-9)
  )
  for (case in cases) {
    obj <- do.call(new_pcens, c(list(case[[1]], dunif), case[[2]]))
    expect_close(
      pcens_cdf(obj, case[[3]], case[[4]]),
      do.call(unif_reference, c(list(case[[1]], case[[3]], case[[4]]),
                                case[[2]])),
      info = paste(toString(unlist(case[[2]])), "pwindow", case[[4]])
    )
  }
})

test_that(".check_pwindow errors for missing q and invalid pwindow", {
  expect_null(.check_pwindow(c(1, 2), c(0, 1)))
  expect_error(.check_pwindow(NA_real_, 1), "q must not")
  expect_error(.check_pwindow(1, numeric(0)), "pwindow")
  expect_error(.check_pwindow(1, NA_real_), "pwindow")
  expect_error(.check_pwindow(1, -1), "pwindow")
})

test_that("uniform primary analytical CDFs warn when q and a vector pwindow
  have lengths that do not recycle evenly", {
  for (case in unif_cases()) {
    obj <- do.call(
      new_pcens, c(list(case$pdist, dunif), case$grid[[3]])
    )
    expect_warning(
      pcens_cdf(obj, c(1, 2, 3), c(1, 2)), "multiple"
    )
    expect_warning(
      pcens_cdf(obj, c(1, 2), c(1, 2, 3)), "multiple"
    )
    expect_no_warning(pcens_cdf(obj, c(1, 2, 3, 4), c(1, 2)))
    expect_no_warning(pcens_cdf(obj, 1, c(1, 2, 3)))
    expect_no_warning(pcens_cdf(obj, c(1, 2, 3), 1))
  }
})

test_that("uniform primary analytical CDFs give the delay CDF for the zero
  width elements of a vector pwindow", {
  q <- c(2, 2, 5, -1, 0)
  pwindow <- c(0, 1, 0, 0, 0)
  for (case in unif_cases()) {
    args <- case$grid[[3]]
    obj <- do.call(new_pcens, c(list(case$pdist, dunif), args))
    expected <- c(
      do.call(case$pdist, c(list(2), args)),
      pcens_cdf(obj, 2, 1),
      do.call(case$pdist, c(list(c(5, -1, 0)), args))
    )
    expect_equal(
      pcens_cdf(obj, q, pwindow), expected, tolerance = 1e-14,
      info = case$name
    )
  }
})

test_that("uniform primary analytical CDFs match the empirical CDF of
  rprimarycensored samples", {
  set.seed(378)
  n <- 2e5
  for (case in unif_cases()) {
    for (pwindow in c(0.5, 3)) {
      args <- case$grid[[3]]
      obj <- do.call(new_pcens, c(list(case$pdist, dunif), args))
      draws <- do.call(
        rprimarycensored,
        c(list(n, case$rdist, pwindow = pwindow, swindow = 0), args)
      )
      q <- stats::quantile(draws, c(0.02, 0.1, 0.5, 0.9, 0.98))
      # Within five standard errors of the empirical CDF
      expect_close(
        unname(pcens_cdf(obj, q, pwindow)),
        vapply(q, function(x) mean(draws <= x), numeric(1)),
        rtol = 0, atol = 6e-3,
        info = paste(case$name, "pwindow", pwindow)
      )
    }
  }
})

test_that("the exponential, normal, chi-square and beta delays dispatch to
  analytical pcens_cdf methods", {
  expect_s3_class(new_pcens(pexp, dunif, list()), "pcens_pexp_dunif")
  expect_s3_class(new_pcens(pnorm, dunif, list()), "pcens_pnorm_dunif")
  expect_s3_class(
    new_pcens(pchisq, dunif, list(), df = 1), "pcens_pchisq_dunif"
  )
  expect_s3_class(
    new_pcens(pbeta, dunif, list(), shape1 = 1, shape2 = 1),
    "pcens_pbeta_dunif"
  )
  for (class in c("pexp", "pnorm", "pchisq", "pbeta")) {
    method <- paste0("pcens_", class, "_dunif")
    expect_false(
      is.null(utils::getS3method("pcens_cdf", method, optional = TRUE)),
      info = method
    )
  }
})

test_that("the exponential and chi-square methods agree with the gamma", {
  q <- c(0.001, 0.5, 2, 10, 40)
  for (pwindow in c(0.5, 2, 7)) {
    expect_equal(
      pcens_cdf(new_pcens(pexp, dunif, list(), rate = 0.3), q, pwindow),
      pcens_cdf(
        new_pcens(pgamma, dunif, list(), shape = 1, rate = 0.3), q, pwindow
      ),
      tolerance = 1e-10
    )
  }
  for (df in c(0.5, 3, 40)) {
    expect_identical(
      pcens_cdf(new_pcens(pchisq, dunif, list(), df = df), q, 3),
      pcens_cdf(
        new_pcens(pgamma, dunif, list(), shape = df / 2, scale = 2), q, 3
      ),
      info = paste("df", df)
    )
  }
})

test_that("the exponential and normal methods use the pdist defaults", {
  q <- c(-2, 0.5, 3, 10)
  expect_identical(
    pcens_cdf(new_pcens(pexp, dunif, list()), q, 2),
    pcens_cdf(new_pcens(pexp, dunif, list(), rate = 1), q, 2)
  )
  expect_identical(
    pcens_cdf(new_pcens(pnorm, dunif, list()), q, 2),
    pcens_cdf(new_pcens(pnorm, dunif, list(), mean = 0, sd = 1), q, 2)
  )
})

test_that("the normal method handles infinite q and is continuous across the
  switch to the asymptotic series", {
  obj <- new_pcens(pnorm, dunif, list(), mean = 1, sd = 2)
  expect_identical(pcens_cdf(obj, c(-Inf, Inf), 1), c(0, 1))
  expect_identical(
    pcens_cdf(obj, c(-Inf, 2, Inf), c(1, 0, 0)), c(0, pnorm(2, 1, 2), 1)
  )
  obj <- new_pcens(pnorm, dunif, list(), mean = 0, sd = 1)
  expect_equal(
    pcens_cdf(obj, -10 + 1e-9 - 1e-7, 1), pcens_cdf(obj, -10 + 1e-9 + 1e-7, 1),
    tolerance = 1e-5
  )
  z <- c(-10.000001, -9.999999)
  expect_rel_equal(
    .norm_shortfall(z),
    vapply(z, function(zi) {
      stats::integrate(
        pnorm, -Inf, zi,
        rel.tol = 1.2e-14, abs.tol = 0, subdivisions = 500L
      )$value
    }, numeric(1)),
    tolerance = 1e-12
  )
})

test_that("the exponential, normal, chi-square and beta methods use the
  numerical path for arguments they do not handle, rather than ignoring them", {
  q <- c(0.5, 2, 6)
  # A non-zero ncp, lower.tail and log.p change what pdist returns
  objs <- list(
    new_pcens(pexp, dunif, list(), rate = 0.4, lower.tail = FALSE),
    new_pcens(pnorm, dunif, list(), mean = 2, log.p = TRUE),
    new_pcens(pchisq, dunif, list(), df = 3, lower.tail = FALSE),
    new_pcens(pchisq, dunif, list(), df = 3, ncp = 1),
    new_pcens(pbeta, dunif, list(), shape1 = 2, shape2 = 3, lower.tail = FALSE),
    new_pcens(pbeta, dunif, list(), shape1 = 2, shape2 = 3, ncp = 1)
  )
  for (obj in objs) {
    expect_identical(
      pcens_cdf(obj, q, 0.5), pcens_cdf(obj, q, 0.5, use_numeric = TRUE),
      info = toString(names(obj$args))
    )
  }
  expect_gt(
    max(abs(
      pcens_cdf(objs[[1]], q, 2) -
        pcens_cdf(new_pcens(pexp, dunif, list(), rate = 0.4), q, 2)
    )),
    0.1
  )
  # An unknown argument is an error in pdist
  exp_scale <- new_pcens(pexp, dunif, list(), scale = 2)
  expect_error(pcens_cdf(exp_scale, q, 2), "unused argument")
  # A zero ncp is the central distribution
  expect_identical(
    pcens_cdf(new_pcens(pchisq, dunif, list(), df = 3, ncp = 0), q, 1),
    pcens_cdf(new_pcens(pchisq, dunif, list(), df = 3), q, 1)
  )
})

test_that("the normal method evaluates vector means per point and is
  centred on the mean", {
  q <- c(-3, 0.5, 2, 6)
  mean <- c(-1, 0, 2, 5)
  expect_identical(
    pcens_cdf(new_pcens(pnorm, dunif, list(), mean = mean, sd = 2), q, 1),
    vapply(seq_along(q), function(i) {
      obj_i <- new_pcens(pnorm, dunif, list(), mean = mean[i], sd = 2)
      pcens_cdf(obj_i, q[i], 1)
    }, numeric(1))
  )
  # The window is far narrower than the delay but not than the sd
  far <- new_pcens(pnorm, dunif, list(), mean = 1e7, sd = 1)
  centred <- new_pcens(pnorm, dunif, list(), mean = 0, sd = 1)
  expect_equal(
    pcens_cdf(far, 1e7 + c(-1, 0, 2), 3), pcens_cdf(centred, c(-1, 0, 2), 3),
    tolerance = 1e-8
  )
})

test_that("the chi-square and beta methods error for missing parameters", {
  expect_error(
    pcens_cdf(new_pcens(pchisq, dunif, list()), 1, 1),
    "df parameter is required for Chi-square distribution"
  )
  for (obj in list(
    new_pcens(pbeta, dunif, list(), shape2 = 2),
    new_pcens(pbeta, dunif, list(), shape1 = 2)
  )) {
    expect_error(
      pcens_cdf(obj, 1, 1),
      "shape1 and shape2 parameters are required for Beta distribution"
    )
  }
})

test_that("pprimarycensored and dprimarycensored work end to end for the
  exponential and normal delays", {
  expect_identical(
    pprimarycensored(c(2, 5, 9), pexp, pwindow = 2, rate = 0.3),
    pcens_cdf(new_pcens(pexp, dunif, list(), rate = 0.3), c(2, 5, 9), 2)
  )
  pmf <- dprimarycensored(
    0:15, pnorm,
    pwindow = 2, swindow = 1, L = 0, D = 16, mean = 6, sd = 2
  )
  expect_equal(sum(pmf), 1, tolerance = 1e-10)
  # The window is not clipped at 0, so a normal delay has mass below 0
  lower <- pprimarycensored(-1, pnorm, pwindow = 2, mean = 2, sd = 2)
  expect_gt(lower, 0)
})

test_that("the Pareto method dispatches, uses the extraDistr defaults and
  uses the numerical path for arguments it does not handle", {
  skip_if_not_installed("extraDistr")
  obj <- new_pcens(extraDistr::ppareto, dunif, list(), a = 1, b = 1)
  expect_s3_class(obj, "pcens_ppareto_dunif")
  named <- add_name_attribute(extraDistr::ppareto, "pareto")
  expect_s3_class(new_pcens(named, dunif, list()), "pcens_ppareto_dunif")
  q <- c(0.5, 1.5, 3, 10)
  expect_identical(
    pcens_cdf(new_pcens(extraDistr::ppareto, dunif, list()), q, 2),
    pcens_cdf(obj, q, 2)
  )
  upper <- new_pcens(
    extraDistr::ppareto, dunif, list(), a = 2, b = 1, lower.tail = FALSE
  )
  expect_identical(
    pcens_cdf(upper, q, 0.5), pcens_cdf(upper, q, 0.5, use_numeric = TRUE)
  )
})

test_that("the Pareto method is accurate just above the minimum, for shapes
  near 1 and for vector shapes", {
  skip_if_not_installed("extraDistr")
  b <- 0.5
  d <- b * (1 + c(1e-6, 1e-3, 0.05, 0.5, 4))
  for (a in c(0.3, 1 - 1e-9, 1, 1 + 1e-9, 3)) {
    obj <- new_pcens(extraDistr::ppareto, dunif, list(), a = a, b = b)
    for (pwindow in c(1e-3, 0.5, 2)) {
      expect_rel_equal(
        pcens_cdf(obj, d, pwindow),
        reference_uniform_cdf(
          function(t) extraDistr::ppareto(t, a, b), d, pwindow, kinks = b
        ),
        info = paste("a", a, "pwindow", pwindow)
      )
    }
  }
  a <- c(0.5, 1, 2, 4, 0.9)
  q <- c(0.6, 1, 3, 8, 50)
  obj <- new_pcens(extraDistr::ppareto, dunif, list(), a = a, b = b)
  expect_identical(
    pcens_cdf(obj, q, 1),
    vapply(seq_along(q), function(i) {
      obj_i <- new_pcens(extraDistr::ppareto, dunif, list(), a = a[i], b = b)
      pcens_cdf(obj_i, q[i], 1)
    }, numeric(1))
  )
})

test_that("the inverse gamma and inverse chi-square methods dispatch, use
  the extraDistr defaults and map to the inverse gamma", {
  skip_if_not_installed("extraDistr")
  q <- c(0.2, 1, 3, 10)
  invgamma <- new_pcens(extraDistr::pinvgamma, dunif, list(), alpha = 3)
  expect_s3_class(invgamma, "pcens_pinvgamma_dunif")
  expect_identical(
    pcens_cdf(invgamma, q, 2),
    pcens_cdf(
      new_pcens(extraDistr::pinvgamma, dunif, list(), alpha = 3, beta = 1),
      q, 2
    )
  )
  invchisq <- new_pcens(extraDistr::pinvchisq, dunif, list(), nu = 5)
  expect_s3_class(invchisq, "pcens_pinvchisq_dunif")
  expect_identical(
    pcens_cdf(invchisq, q, 2),
    pcens_cdf(
      new_pcens(extraDistr::pinvgamma, dunif, list(), alpha = 2.5, beta = 0.5),
      q, 2
    )
  )
  scaled <- new_pcens(extraDistr::pinvchisq, dunif, list(), nu = 5, tau = 2)
  expect_identical(
    pcens_cdf(scaled, q, 2),
    pcens_cdf(
      new_pcens(extraDistr::pinvgamma, dunif, list(), alpha = 2.5, beta = 5),
      q, 2
    )
  )
  for (name in c("inverse gamma", "invgamma")) {
    named <- add_name_attribute(extraDistr::pinvgamma, name)
    expect_s3_class(
      new_pcens(named, dunif, list(), alpha = 3), "pcens_pinvgamma_dunif"
    )
  }
  for (name in c("invchisq", "scaleinvchisq")) {
    named <- add_name_attribute(extraDistr::pinvchisq, name)
    expect_s3_class(
      new_pcens(named, dunif, list(), nu = 3), "pcens_pinvchisq_dunif"
    )
  }
})

test_that("the inverse gamma and inverse chi-square methods use the
  numerical path without a finite mean or for arguments they do not handle", {
  skip_if_not_installed("extraDistr")
  q <- c(0.5, 2, 6)
  objs <- list(
    new_pcens(extraDistr::pinvgamma, dunif, list(), alpha = 1, beta = 2),
    new_pcens(extraDistr::pinvgamma, dunif, list(), alpha = 0.6),
    new_pcens(
      extraDistr::pinvgamma, dunif, list(), alpha = c(3, 0.8, 2), beta = 1
    ),
    new_pcens(
      extraDistr::pinvgamma, dunif, list(), alpha = 3, lower.tail = FALSE
    ),
    new_pcens(extraDistr::pinvchisq, dunif, list(), nu = 2),
    new_pcens(extraDistr::pinvchisq, dunif, list(), nu = 1.5, tau = 1)
  )
  for (obj in objs) {
    expect_identical(
      pcens_cdf(obj, q, 0.5), pcens_cdf(obj, q, 0.5, use_numeric = TRUE),
      info = toString(unlist(obj$args))
    )
  }
  expect_error(
    pcens_cdf(new_pcens(extraDistr::pinvgamma, dunif, list(), beta = 1), 1, 1),
    "alpha parameter is required for inverse gamma distribution"
  )
  expect_error(
    pcens_cdf(new_pcens(extraDistr::pinvchisq, dunif, list()), 1, 1),
    "nu parameter is required for inverse chi-square distribution"
  )
})

test_that("the inverse gamma method is accurate just above shape 1 and deep
  in the lower tail", {
  skip_if_not_installed("extraDistr")
  for (args in list(list(1.001, 1), list(3, 50), list(30, 5))) {
    obj <- new_pcens(
      extraDistr::pinvgamma, dunif, list(), alpha = args[[1]], beta = args[[2]]
    )
    cdf <- function(t) extraDistr::pinvgamma(t, args[[1]], args[[2]])
    for (pwindow in c(0.5, 2)) {
      d <- c(0.05, 0.3, 1, 10, 100)
      expect_rel_equal(
        pcens_cdf(obj, d, pwindow),
        reference_uniform_cdf(cdf, d, pwindow, kinks = 0),
        tolerance = 1e-8,
        info = paste(toString(unlist(args)), "pwindow", pwindow)
      )
    }
  }
})
