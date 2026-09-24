# Pin the results of dprimarycensored(), pprimarycensored(), new_pcens(),
# pcens_pmf() and the fitdistdoublecens() likelihood wrappers against
# reference implementations built from pcens_cdf(), so refactors that
# reduce overhead cannot change results.

equiv_delays <- list(
  gamma = list(pdist = pgamma, name = "pgamma", args = list(
    shape = 3, scale = 2
  )),
  lnorm = list(pdist = plnorm, name = "plnorm", args = list(
    meanlog = 1, sdlog = 0.5
  )),
  weibull = list(pdist = pweibull, name = "pweibull", args = list(
    shape = 1.5, scale = 3
  ))
)

equiv_primaries <- list(
  unif = list(
    dprimary = dunif, pprimary = punif, name = "dunif", args = list()
  ),
  expgrowth = list(
    dprimary = dexpgrowth, pprimary = pexpgrowth, name = "dexpgrowth",
    args = list(r = 0.2)
  )
)

equiv_bounds <- list(
  list(L = -Inf, D = Inf),
  list(L = -Inf, D = 15),
  list(L = 1, D = 15),
  list(L = 1, D = Inf)
)

# A pcens object built by hand, without any name lookups
ref_pcens <- function(delay, primary) {
  obj <- list(
    pdist = delay$pdist,
    dprimary = primary$dprimary,
    primary_args = primary$args,
    dprimary_args = primary$args,
    pprimary = primary$pprimary,
    args = delay$args
  )
  class(obj) <- c(
    sprintf("pcens_%s_%s", delay$name, primary$name),
    sprintf("pcens_%s", delay$name),
    "pcens"
  )
  obj
}

# Reference CDF at a single truncation point
ref_cdf_at <- function(obj, bound, pwindow, inf_value) {
  if (is.infinite(bound)) {
    return(inf_value)
  }
  pcens_cdf(obj, bound, pwindow)
}

# Reference PMF: difference the CDF at x and min(x + swindow, D) and
# normalise over [L, D]
ref_pmf <- function(obj, x, pwindow, swindow, L, D, log = FALSE) {
  upper <- pmin(x + swindow, D)
  pts <- sort(unique(c(x, upper)))
  cdfs <- pcens_cdf(obj, pts, pwindow)
  cdfs[pts == -Inf] <- 0
  cdfs[pts == Inf] <- 1
  result <- cdfs[match(upper, pts)] - cdfs[match(x, pts)]
  if (!(is.infinite(L) && is.infinite(D))) {
    normaliser <- ref_cdf_at(obj, D, pwindow, 1) -
      ref_cdf_at(obj, L, pwindow, 0)
    if (normaliser != 1) {
      result <- result / normaliser
    }
  }
  result <- pmax(0, result)
  if (log) log(result) else result
}

# Reference CDF normalised over [L, D] and clamped outside it
ref_cdf <- function(obj, q, pwindow, L, D) {
  result <- pcens_cdf(obj, q, pwindow)
  cdf_D <- ref_cdf_at(obj, D, pwindow, 1)
  cdf_L <- ref_cdf_at(obj, L, pwindow, 0)
  normaliser <- cdf_D - cdf_L
  if (!(cdf_L == 0 && normaliser == 1)) {
    result <- (result - cdf_L) / normaliser
  }
  result <- ifelse(q <= L, 0, result)
  ifelse(q >= D, 1, result)
}

# Loop over every delay, primary, truncation and window combination
for_each_equiv_case <- function(fn) {
  for (delay in equiv_delays) {
    for (primary in equiv_primaries) {
      for (bounds in equiv_bounds) {
        for (pwindow in c(1, 2)) {
          for (swindow in c(1, 2)) {
            lower <- if (is.finite(bounds$L)) bounds$L else 0
            upper <- if (is.finite(bounds$D)) bounds$D - 0.5 else 20
            fn(
              delay = delay, primary = primary, L = bounds$L, D = bounds$D,
              pwindow = pwindow, swindow = swindow,
              x = seq(lower, upper, by = 0.5)
            )
          }
        }
      }
    }
  }
}

test_that("new_pcens matches a hand-built object", {
  for (delay in equiv_delays) {
    for (primary in equiv_primaries) {
      obj <- do.call(
        new_pcens,
        c(
          list(delay$pdist, primary$dprimary, primary_args = primary$args),
          delay$args
        )
      )
      expect_identical(obj, ref_pcens(delay, primary))
    }
  }
})

test_that("dprimarycensored and pcens_pmf match the reference PMF", {
  for_each_equiv_case(function(delay, primary, L, D, pwindow, swindow, x) {
    obj <- ref_pcens(delay, primary)
    expected <- ref_pmf(obj, x, pwindow, swindow, L, D)
    actual <- suppressMessages(do.call(
      dprimarycensored,
      c(
        list(
          x, delay$pdist,
          pwindow = pwindow, swindow = swindow, L = L, D = D,
          dprimary = primary$dprimary, primary_args = primary$args
        ),
        delay$args
      )
    ))
    expect_identical(actual, expected)
    expect_identical(
      suppressMessages(pcens_pmf(obj, x, pwindow, swindow, L, D)),
      expected
    )
    expect_identical(
      suppressMessages(pcens_pmf(obj, x, pwindow, swindow, L, D, log = TRUE)),
      log(expected)
    )
  })
})

test_that("pprimarycensored matches the reference CDF", {
  for_each_equiv_case(function(delay, primary, L, D, pwindow, swindow, x) {
    if (swindow != 1) {
      return(invisible(NULL))
    }
    qs <- c(x, D)
    obj <- ref_pcens(delay, primary)
    actual <- do.call(
      pprimarycensored,
      c(
        list(
          qs, delay$pdist,
          pwindow = pwindow, L = L, D = D,
          dprimary = primary$dprimary, primary_args = primary$args
        ),
        delay$args
      )
    )
    expect_identical(actual, ref_cdf(obj, qs, pwindow, L, D))
  })
})

test_that("pcens_pmf handles non-integer x and duplicated points", {
  obj <- ref_pcens(equiv_delays$gamma, equiv_primaries$unif)
  x <- c(3.2, 0.1, 3.2, 7.75, 0.1, 12)
  for (swindow in c(0.3, 1, 2.5)) {
    for (D in c(Inf, 13)) {
      expect_identical(
        suppressMessages(pcens_pmf(obj, x, 1, swindow, D = D)),
        ref_pmf(obj, x, 1, swindow, -Inf, D)
      )
    }
  }
})

equiv_fit_params <- data.frame(
  swindow = c(1, 1, 2, 2, 1, 1),
  pwindow = c(1, 1, 1, 2, 2, 2),
  L = c(-Inf, -Inf, -Inf, 0, 0, 0),
  D = c(Inf, Inf, 15, 15, 20, 20)
)

test_that(".dpcens matches dprimarycensored per observation", {
  x <- c(1, 4, 2, 6, 3, 8)
  for (delay in equiv_delays) {
    for (primary in equiv_primaries) {
      expected <- vapply(seq_along(x), function(i) {
        do.call(
          dprimarycensored,
          c(
            list(
              x[i], delay$pdist,
              pwindow = equiv_fit_params$pwindow[i],
              swindow = equiv_fit_params$swindow[i],
              L = equiv_fit_params$L[i], D = equiv_fit_params$D[i],
              dprimary = primary$dprimary, primary_args = primary$args
            ),
            delay$args
          )
        )
      }, numeric(1))
      actual <- do.call(
        .dpcens,
        c(
          list(
            x, equiv_fit_params, delay$pdist, primary$dprimary,
            primary$args
          ),
          delay$args
        )
      )
      expect_identical(actual, expected)
    }
  }
})

test_that(".ppcens matches pprimarycensored per observation", {
  q <- c(1, 4, 2, 6, 3, 8)
  for (delay in equiv_delays) {
    for (primary in equiv_primaries) {
      expected <- vapply(seq_along(q), function(i) {
        do.call(
          pprimarycensored,
          c(
            list(
              q[i], delay$pdist,
              pwindow = equiv_fit_params$pwindow[i],
              L = equiv_fit_params$L[i], D = equiv_fit_params$D[i],
              dprimary = primary$dprimary, primary_args = primary$args
            ),
            delay$args
          )
        )
      }, numeric(1))
      actual <- do.call(
        .ppcens,
        c(
          list(
            q, equiv_fit_params, delay$pdist, primary$dprimary,
            primary$args
          ),
          delay$args
        )
      )
      expect_identical(actual, expected)
    }
  }
})

test_that("fitdistdoublecens matches a fit built on dprimarycensored", {
  skip_if_not_installed("fitdistrplus")
  withr::local_seed(101)
  n <- 100
  pwindows <- rep_len(c(1, 2), n)
  Ds <- rep(c(15, Inf), each = n / 2)
  delays <- vapply(seq_len(n), function(i) {
    rprimarycensored(
      1, rgamma,
      shape = 3, scale = 2, pwindow = pwindows[i], swindow = 1, D = Ds[i]
    )
  }, numeric(1))
  censdata <- data.frame(
    left = delays, right = delays + 1, pwindow = pwindows, D = Ds
  )
  start <- list(shape = 1, scale = 1)

  fit <- fitdistdoublecens(
    censdata, "gamma",
    start = start, truncation_check_multiplier = NULL
  )

  # Reference: one dprimarycensored() call per observation
  dref <- function(x, shape, scale) {
    vapply(seq_along(x), function(i) {
      tryCatch(
        dprimarycensored(
          x[i], pgamma,
          pwindow = pwindows[i], swindow = 1, D = Ds[i],
          shape = shape, scale = scale, check = FALSE
        ),
        error = function(e) NaN
      )
    }, numeric(1))
  }
  pref <- function(q, shape, scale) {
    pgamma(q, shape = shape, scale = scale)
  }
  ref_env <- new.env()
  ref_env$ddref <- dref
  ref_env$pdref <- pref
  ref <- withr::with_environment(
    ref_env,
    fitdistrplus::fitdist(delays, "dref", start = start)
  )

  expect_equal(fit$estimate, ref$estimate, tolerance = 1e-10)
  expect_equal(fit$loglik, ref$loglik, tolerance = 1e-10)
})

test_that("a primary CDF found by alias fails the name check", {
  dprim <- add_name_attribute(
    function(x, min, max) dunif(x, min, max), "uniform"
  )
  msg <- "refer to different distributions: 'uniform' vs 'punif'"
  # new_pcens() does not check a primary CDF it looks up itself
  obj <- new_pcens(pgamma, dprim, shape = 1, scale = 1)
  expect_identical(obj$pprimary, punif)
  expect_error(
    dprimarycensored(1:3, pgamma, dprimary = dprim, shape = 1, scale = 1),
    msg
  )
  expect_error(
    pprimarycensored(1:3, pgamma, dprimary = dprim, shape = 1, scale = 1),
    msg
  )
})

test_that(".dpcens and .ppcens give NaN when any input errors", {
  # fitdistrplus probes the fitted functions with inputs like these
  params <- data.frame(
    swindow = 1, pwindow = rep(c(1, 2), 3), L = -Inf,
    D = rep(c(Inf, 10), each = 3)
  )
  for (x in list(c(0, 1, NA), c(0, 1, Inf, NaN, -1))) {
    nans <- rep(NaN, length(x))
    expect_identical(
      .dpcens(x, params, pgamma, dunif, list(), shape = 2, scale = 1), nans
    )
    expect_identical(
      .ppcens(x, params, pgamma, dunif, list(), shape = 2, scale = 1), nans
    )
  }
  expect_length(
    .ppcens(numeric(0), params, pgamma, dunif, list(), shape = 2, scale = 1),
    0L
  )
})
