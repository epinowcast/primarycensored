loglik_cases <- list(
  gamma_unif = list(
    pdist = pgamma, dprimary = dunif, primary_args = list(),
    pars = list(shape = 2.5, scale = 1.5)
  ),
  lnorm_unif = list(
    pdist = plnorm, dprimary = dunif, primary_args = list(),
    pars = list(meanlog = 1.2, sdlog = 0.6)
  ),
  weibull_unif = list(
    pdist = pweibull, dprimary = dunif, primary_args = list(),
    pars = list(shape = 1.8, scale = 4)
  ),
  gamma_expgrowth = list(
    pdist = pgamma, dprimary = dexpgrowth, primary_args = list(r = 0.2),
    pars = list(shape = 2.5, scale = 1.5)
  )
)

# Per-row reference: one dprimarycensored() call for every observation
reference_loglik <- function(case, x, pwindow, swindow, L, D) {
  n <- length(x)
  pwindow <- rep_len(pwindow, n)
  swindow <- rep_len(swindow, n)
  L <- rep_len(L, n)
  D <- rep_len(D, n)
  vapply(seq_len(n), function(i) {
    suppressMessages(log(do.call(
      dprimarycensored,
      c(
        list(
          x[i], case$pdist,
          pwindow = pwindow[i], swindow = swindow[i], L = L[i], D = D[i],
          dprimary = case$dprimary, primary_args = case$primary_args
        ),
        case$pars
      )
    )))
  }, numeric(1))
}

make_loglik <- function(case, x, ...) {
  pcens_loglik_fn(
    x, case$pdist,
    dprimary = case$dprimary, primary_args = case$primary_args, ...
  )
}

# Agreement with dprimarycensored() is to a relative tolerance of 1e-10, as
# both paths evaluate the same CDFs.
tol <- 1e-10

test_that("pcens_loglik_fn returns a function of the delay parameters", {
  ll <- pcens_loglik_fn(0:5, pgamma)
  expect_type(ll, "closure")
  out <- ll(shape = 2, scale = 1)
  expect_type(out, "double")
  expect_length(out, 6)
})

test_that("pcens_loglik_fn matches log(dprimarycensored) per row", {
  x <- c(0, 0.5, 1, 2, 3.25, 7)
  for (case in loglik_cases) {
    ll <- make_loglik(case, x, pwindow = 1, swindow = 1)
    expect_equal(
      do.call(ll, case$pars),
      reference_loglik(case, x, 1, 1, -Inf, Inf),
      tolerance = tol
    )
  }
})

test_that("pcens_loglik_fn handles truncation", {
  x <- c(1, 2, 3, 4, 6, 9)
  for (case in loglik_cases) {
    ll <- make_loglik(case, x, pwindow = 2, swindow = 1, L = 1, D = 10)
    expect_equal(
      do.call(ll, case$pars),
      reference_loglik(case, x, 2, 1, 1, 10),
      tolerance = tol
    )
  }
})

test_that("pcens_loglik_fn supports per-row windows and truncation", {
  x <- c(0, 1, 1, 2, 3, 3, 5, 0.5)
  pwindow <- c(1, 1, 2, 2, 1, 2, 1, 1)
  swindow <- c(1, 2, 1, 1, 1, 2, 2, 1)
  L <- c(-Inf, -Inf, 0, 0, -Inf, 0, -Inf, -Inf)
  D <- c(Inf, 10, 10, Inf, 10, Inf, 10, Inf)
  for (case in loglik_cases) {
    ll <- make_loglik(
      case, x,
      pwindow = pwindow, swindow = swindow, L = L, D = D
    )
    expect_equal(
      do.call(ll, case$pars),
      reference_loglik(case, x, pwindow, swindow, L, D),
      tolerance = tol
    )
  }
})

test_that("pcens_loglik_fn returns values in row order with repeated x", {
  set.seed(1)
  x <- sample(0:8, 60, replace = TRUE)
  pwindow <- sample.int(2, 60, replace = TRUE)
  case <- loglik_cases$gamma_unif
  ll <- make_loglik(case, x, pwindow = pwindow, swindow = 1, D = 10)
  out <- do.call(ll, case$pars)
  expect_equal(
    out, reference_loglik(case, x, pwindow, 1, -Inf, 10),
    tolerance = tol
  )
  # Rows with the same x and settings get the same value
  expect_identical(out[x == 3 & pwindow == 1], rep(
    out[x == 3 & pwindow == 1][1], sum(x == 3 & pwindow == 1)
  ))
})

test_that("pcens_loglik_fn evaluates each group at its unique x", {
  seen <- new.env(parent = emptyenv())
  seen$n <- integer(0)
  pgamma_counted <- add_name_attribute(
    function(q, shape, scale) {
      seen$n <- c(seen$n, length(q))
      pgamma(q, shape = shape, scale = scale)
    },
    "pgamma"
  )
  x <- rep(0:4, times = 200)
  ll <- pcens_loglik_fn(x, pgamma_counted, pwindow = 1, swindow = 1)
  seen$n <- integer(0)
  out <- ll(shape = 2, scale = 1)
  expect_length(out, 1000)
  # 5 unique x give 6 unique CDF points, not 2000
  expect_gt(length(seen$n), 0)
  expect_lte(max(seen$n), 6)
})

test_that("pcens_loglik_fn handles exact windows and mixed swindow", {
  x <- c(1, 2, 3, 4)
  swindow <- c(1, 0, 1, 0)
  for (case in loglik_cases) {
    ll <- make_loglik(case, x, pwindow = 1, swindow = swindow)
    expect_equal(
      do.call(ll, case$pars),
      reference_loglik(case, x, 1, swindow, -Inf, Inf),
      tolerance = tol
    )
  }
  case <- loglik_cases$lnorm_unif
  ll <- make_loglik(case, x, pwindow = c(0, 1, 0, 2), swindow = 1)
  expect_equal(
    do.call(ll, case$pars),
    reference_loglik(case, x, c(0, 1, 0, 2), 1, -Inf, Inf),
    tolerance = tol
  )
})

test_that("pcens_loglik_fn returns -Inf for zero probability rows", {
  ll <- pcens_loglik_fn(c(0, 50), pgamma, pwindow = 1, swindow = 1, D = 100)
  out <- suppressWarnings(ll(shape = 2, scale = 0.01))
  expect_identical(out[2], -Inf)
  expect_true(is.finite(out[1]))
})

test_that("pcens_loglik_fn takes fixed delay parameters through ...", {
  x <- c(0, 1, 2, 3)
  ll <- pcens_loglik_fn(x, plnorm, sdlog = 0.6, pwindow = 1, swindow = 1)
  expect_identical(
    ll(meanlog = 1.2),
    pcens_loglik_fn(x, plnorm, pwindow = 1, swindow = 1)(
      meanlog = 1.2, sdlog = 0.6
    )
  )
  # A parameter given to the function overrides the fixed one
  expect_identical(
    ll(meanlog = 1.2, sdlog = 0.3),
    pcens_loglik_fn(x, plnorm, pwindow = 1, swindow = 1)(
      meanlog = 1.2, sdlog = 0.3
    )
  )
})

test_that("pcens_loglik_fn calls are independent of each other", {
  x <- c(0, 1, 2, 3, 4)
  ll <- pcens_loglik_fn(x, pgamma, pwindow = 1, swindow = 1)
  a <- ll(shape = 2, rate = 1)
  b <- ll(shape = 3, rate = 0.5)
  expect_identical(ll(shape = 2, rate = 1), a)
  expect_false(identical(a, b))
  # Changing the parameterisation does not leave the old parameters behind
  expect_equal(
    ll(shape = 2, scale = 1),
    pcens_loglik_fn(x, pgamma, pwindow = 1, swindow = 1)(shape = 2, scale = 1),
    tolerance = tol
  )
})

test_that("pcens_loglik_fn checks parameter names on every new name set", {
  ll <- pcens_loglik_fn(0:4, pgamma, pwindow = 1, swindow = 1)
  expect_no_error(ll(shape = 2, rate = 1))
  expect_error(ll(shape = 2, rte = 1), "Unknown delay parameter")
  expect_error(ll(2, 1), "must be named")
  # A valid call still works after the failures
  expect_no_error(ll(shape = 3, rate = 1))
})

test_that("pcens_loglik_fn validates pdist on first use", {
  decreasing <- function(q, shape, rate) exp(-q / 1000)
  ll <- pcens_loglik_fn(0:4, decreasing, pwindow = 1, swindow = 1)
  expect_error(
    ll(shape = 2, rate = 1),
    "pdist is not a valid cumulative distribution function"
  )
  ll <- pcens_loglik_fn(
    0:4, decreasing,
    pwindow = 1, swindow = 1, check = FALSE
  )
  expect_no_error(suppressWarnings(ll(shape = 2, rate = 1)))
})

test_that("pcens_loglik_fn validates dprimary at construction", {
  bad_dprimary <- function(x, min, max) dunif(x, min, max) * 2
  expect_error(
    pcens_loglik_fn(0:4, pgamma, dprimary = bad_dprimary),
    "dprimary is not a valid probability density function"
  )
  expect_error(
    pcens_loglik_fn(0:4, pgamma, dprimary = function(x) x),
    "dprimary must take min and max"
  )
  expect_no_error(
    pcens_loglik_fn(0:4, pgamma, dprimary = bad_dprimary, check = FALSE)
  )
})

test_that("pcens_loglik_fn validates its inputs at construction", {
  expect_error(pcens_loglik_fn("a", pgamma), "x must be numeric")
  expect_error(pcens_loglik_fn(c(1, NA), pgamma), "missing values")
  expect_error(
    pcens_loglik_fn(0:4, pgamma, pwindow = c(1, 2)),
    "length 1 or"
  )
  expect_error(
    pcens_loglik_fn(0:4, pgamma, D = c(10, 11, 12)),
    "length 1 or"
  )
  expect_error(pcens_loglik_fn(0:4, pgamma, pwindow = -1), "non-negative")
  expect_error(pcens_loglik_fn(0:4, pgamma, swindow = -1), "non-negative")
  expect_error(
    pcens_loglik_fn(0:4, pgamma, L = 5, D = 5), "L must be less than D"
  )
  expect_error(
    pcens_loglik_fn(0:4, pgamma, L = 1),
    "below L"
  )
  expect_error(
    pcens_loglik_fn(0:4, pgamma, D = 4),
    "Upper truncation point is greater than D"
  )
})

test_that("pcens_loglik_fn rejects non-numeric settings", {
  for (nm in c("pwindow", "swindow", "L", "D")) {
    args <- list(x = 0:3, pdist = plnorm)
    args[[nm]] <- "a"
    expect_error(do.call(pcens_loglik_fn, args), paste(nm, "must be numeric"))
  }
  # A function passed by position lands in pwindow
  expect_error(pcens_loglik_fn(0:3, plnorm, dunif), "pwindow must be numeric")
})

test_that("pcens_loglik_fn supports non-parametric delays", {
  boundaries <- 0:5
  pmf1 <- c(0.1, 0.2, 0.3, 0.25, 0.15)
  pmf2 <- rev(pmf1)
  x <- c(0, 1, 1, 2, 3, 4)
  ll <- pcens_loglik_fn(x, pdiscretestep, pwindow = 1, swindow = 1)
  reference <- function(pmf) {
    vapply(x, function(xi) {
      log(dprimarycensored(
        xi, pdiscretestep,
        pwindow = 1, swindow = 1, boundaries = boundaries, pmf = pmf
      ))
    }, numeric(1))
  }
  expect_equal(
    ll(pmf = pmf1, boundaries = boundaries), reference(pmf1),
    tolerance = tol
  )
  # Same names, new values: the unchecked update path
  expect_equal(
    ll(pmf = pmf2, boundaries = boundaries), reference(pmf2),
    tolerance = tol
  )
  expect_equal(
    ll(pmf = pmf1, boundaries = boundaries), reference(pmf1),
    tolerance = tol
  )
})

test_that("pcens_loglik_fn evaluates the CDF once per unique pwindow", {
  calls <- new.env(parent = emptyenv())
  calls$cdf <- 0L
  local_mocked_bindings(
    pcens_cdf = function(object, q, pwindow, ...) {
      calls$cdf <- calls$cdf + 1L
      pgamma(q, shape = 2, scale = 1)
    }
  )
  x <- rep(0:4, times = 20)
  # Truncation points that are not delays are evaluated in the same call
  ll <- suppressMessages(
    pcens_loglik_fn(x, pgamma, pwindow = 1, swindow = 1, L = -1, D = 10)
  )
  out <- ll(shape = 2, scale = 1)
  expect_length(out, 100)
  expect_identical(calls$cdf, 1L)
  calls$cdf <- 0L
  ll <- pcens_loglik_fn(
    x, pgamma,
    pwindow = rep(1:2, each = 50), swindow = 1, D = Inf
  )
  ll(shape = 2, scale = 1)
  expect_identical(calls$cdf, 2L)
})

test_that("pcens_loglik_fn shares CDF points across groups", {
  calls <- new.env(parent = emptyenv())
  calls$cdf <- 0L
  calls$q <- list()
  calls$pwindow <- numeric(0)
  local_mocked_bindings(
    pcens_cdf = function(object, q, pwindow, ...) {
      calls$cdf <- calls$cdf + 1L
      calls$q[[calls$cdf]] <- q
      calls$pwindow[[calls$cdf]] <- pwindow
      pgamma(q, shape = 2, scale = 1)
    }
  )
  set.seed(4)
  n <- 600
  x <- sample(0:20, n, replace = TRUE)
  pwindow <- sample.int(2, n, replace = TRUE)
  swindow <- sample(c(0.5, 1, 2), n, replace = TRUE)
  L <- sample(c(-Inf, 0), n, replace = TRUE)
  D <- sample(c(Inf, 25, 40), n, replace = TRUE)
  ll <- suppressMessages(
    pcens_loglik_fn(x, pgamma,
      pwindow = pwindow, swindow = swindow,
      L = L, D = D
    )
  )
  ll(shape = 2, scale = 1)
  # One call per unique pwindow, however many groups there are
  expect_identical(calls$cdf, 2L)
  expect_setequal(calls$pwindow, 1:2)
  for (i in seq_len(calls$cdf)) {
    q <- calls$q[[i]]
    # Each point is evaluated once
    expect_false(anyDuplicated(q) > 0)
    expect_false(is.unsorted(q))
    # The points are those needed by the rows with this pwindow
    rows <- pwindow == calls$pwindow[[i]]
    upper <- pmin(x[rows] + swindow[rows], D[rows])
    needed <- c(x[rows], upper, L[rows], D[rows])
    expect_setequal(q, unique(needed[is.finite(needed)]))
  }
})

test_that("pcens_loglik_fn matches dprimarycensored with many groups", {
  set.seed(5)
  n <- 300
  x <- sample(0:12, n, replace = TRUE)
  pwindow <- sample.int(2, n, replace = TRUE)
  swindow <- sample(c(0, 0.5, 1, 2), n, replace = TRUE)
  L <- sample(c(-Inf, 0, -0.5), n, replace = TRUE)
  D <- sample(c(Inf, 25, 22.5, 40), n, replace = TRUE)
  for (case in loglik_cases[c("gamma_unif", "lnorm_unif")]) {
    ll <- suppressMessages(make_loglik(
      case, x,
      pwindow = pwindow, swindow = swindow, L = L, D = D
    ))
    expect_equal(
      do.call(ll, case$pars),
      reference_loglik(case, x, pwindow, swindow, L, D),
      tolerance = tol
    )
  }
})

test_that("pcens_loglik_fn returns NaN for invalid parameters", {
  x <- c(0, 1, 2)
  for (D in c(Inf, 10)) {
    ll <- pcens_loglik_fn(x, pgamma, D = D)
    valid <- ll(shape = 2, rate = 1)
    expect_true(all(is.finite(valid)))
    invalid <- suppressWarnings(ll(shape = -1, rate = 1))
    expect_identical(invalid, rep(NaN, 3))
    # A valid call afterwards is unaffected
    expect_identical(ll(shape = 2, rate = 1), valid)
  }
})

test_that("pcens_loglik_fn gives NaN when the first call is invalid", {
  ll <- pcens_loglik_fn(0:2, pgamma)
  expect_identical(
    suppressWarnings(ll(shape = -1, rate = 1)),
    rep(NaN, 3)
  )
  # pdist is still checked once a call gives a valid result
  expect_true(all(is.finite(ll(shape = 2, rate = 1))))
  decreasing <- function(q, shape) rep(2 * shape, length(q))
  bad <- pcens_loglik_fn(0:2, decreasing)
  expect_error(bad(shape = 1), "not a valid cumulative")
})

test_that("pcens_loglik_fn matches pcens_pmf with all truncation forms", {
  case <- loglik_cases$lnorm_unif
  x <- c(1, 2, 3, 4)
  settings <- list(
    list(L = -Inf, D = Inf), list(L = 0, D = Inf), list(L = -Inf, D = 6),
    list(L = 0.5, D = 6), list(L = 0, D = 4.5)
  )
  for (st in settings) {
    for (sw in c(0, 1, 2.5)) {
      ll <- suppressMessages(make_loglik(
        case, x,
        pwindow = 1.5, swindow = sw, L = st$L, D = st$D
      ))
      expect_equal(
        do.call(ll, case$pars),
        reference_loglik(case, x, 1.5, sw, st$L, st$D),
        tolerance = tol
      )
    }
  }
})

test_that("pcens_loglik_fn clips secondary windows at D with one message", {
  x <- c(1, 4, 8, 9.5)
  case <- loglik_cases$gamma_unif
  expect_message(
    make_loglik(case, x, pwindow = 1, swindow = 1, D = 10),
    "clipping"
  )
  ll <- suppressMessages(
    make_loglik(case, x, pwindow = 1, swindow = 1, D = 10)
  )
  expect_no_message(do.call(ll, case$pars))
  out <- do.call(ll, case$pars)
  expect_equal(
    out, reference_loglik(case, x, 1, 1, -Inf, 10),
    tolerance = tol
  )
})

test_that("pcens_loglik_fn accepts a pdist name and a pprimary", {
  x <- 0:5
  expect_identical(
    pcens_loglik_fn(x, "gamma")(shape = 2, scale = 1),
    pcens_loglik_fn(x, pgamma)(shape = 2, scale = 1)
  )
  expect_identical(
    pcens_loglik_fn(x, pgamma, pprimary = "uniform")(shape = 2, scale = 1),
    pcens_loglik_fn(x, pgamma)(shape = 2, scale = 1)
  )
})

test_that("pcens_loglik_fn handles empty and single-row input", {
  expect_identical(
    pcens_loglik_fn(numeric(0), pgamma)(shape = 2, scale = 1),
    numeric(0)
  )
  expect_equal(
    pcens_loglik_fn(2, pgamma)(shape = 2, scale = 1),
    log(dprimarycensored(2, pgamma, shape = 2, scale = 1)),
    tolerance = tol
  )
})

test_that("pcens_loglik_fn is the likelihood of fitdistdoublecens", {
  skip_if_not_installed("fitdistrplus")
  skip_if_not_installed("withr")
  set.seed(42)
  n <- 300
  pw <- sample.int(2, n, replace = TRUE)
  samples <- rprimarycensored(
    n, rgamma,
    shape = 3, rate = 1.2, pwindow = pw, swindow = 1, D = 15
  )
  dat <- data.frame(
    left = samples, right = samples + 1, pwindow = pw, D = 15
  )
  fit <- suppressMessages(fitdistdoublecens(
    dat, "gamma",
    start = list(shape = 2, rate = 1), truncation_check_multiplier = NULL
  ))
  ll <- pcens_loglik_fn(
    dat$left, pgamma,
    pwindow = dat$pwindow, swindow = 1, D = dat$D
  )
  expect_equal(
    sum(do.call(ll, as.list(fit$estimate))), fit$loglik,
    tolerance = 1e-8
  )
})

test_that("pcens_loglik_fn can be optimised directly", {
  set.seed(7)
  n <- 400
  samples <- rprimarycensored(
    n, rlnorm,
    meanlog = 1.3, sdlog = 0.5,
    pwindow = 1, swindow = 1, D = 20
  )
  ll <- pcens_loglik_fn(samples, plnorm, pwindow = 1, swindow = 1, D = 20)
  fit <- optim(
    c(1, 0),
    function(par) -sum(ll(meanlog = par[1], sdlog = exp(par[2]))),
    method = "BFGS"
  )
  expect_equal(fit$par[1], 1.3, tolerance = 0.1)
  expect_equal(exp(fit$par[2]), 0.5, tolerance = 0.15)
})

test_that(".pcens_row_groups collapses equal settings and repeated x", {
  rg <- .pcens_row_groups(
    c(1, 1, 2, 2, 3),
    pwindow = 1, swindow = 1, L = -Inf, D = Inf
  )
  groups <- rg$groups
  expect_length(groups, 1)
  expect_identical(groups[[1]]$x, c(1, 2, 3))
  expect_identical(groups[[1]]$map, c(1L, 1L, 2L, 2L, 3L))
  expect_null(groups[[1]]$idx)
  rg <- .pcens_row_groups(
    c(1, 2, 3, 4),
    pwindow = c(1, 2, 1, 2), swindow = 1, L = -Inf, D = Inf
  )
  groups <- rg$groups
  expect_length(groups, 2)
  expect_identical(groups[[1]]$idx, c(1L, 3L))
  expect_identical(groups[[2]]$idx, c(2L, 4L))
  expect_null(groups[[1]]$map)
  expect_identical(groups[[2]]$pwindow, 2)
  expect_length(.pcens_row_groups(numeric(0), 1, 1, -Inf, Inf)$groups, 0)
})

test_that(".pcens_row_groups shares one sorted set of points per pwindow", {
  rg <- .pcens_row_groups(
    c(3, 1, 2),
    pwindow = 1, swindow = 1, L = 0, D = 3.5
  )
  expect_length(rg$sets, 1)
  pts <- rg$sets[[1]]$points
  g <- rg$groups[[1]]
  # Upper ends are clipped at D, and L and D are in the same set
  expect_identical(pts, c(0, 1, 2, 3, 3.5))
  expect_identical(pts[g$lower], c(3, 1, 2))
  expect_identical(pts[g$upper], c(3.5, 2, 3))
  expect_identical(g$pos_L, 1L)
  expect_identical(g$pos_D, 5L)
  expect_true(g$truncated)
  expect_false(g$exact)
  expect_identical(g$set, 1L)
  rg <- .pcens_row_groups(c(1, 2), 1, swindow = 0, L = -Inf, D = Inf)
  g <- rg$groups[[1]]
  expect_true(g$exact)
  expect_length(rg$sets[[1]]$points, 0)
  expect_null(g$lower)
  expect_true(is.na(g$pos_D))
  expect_false(g$truncated)
  rg <- .pcens_row_groups(1, 1, swindow = Inf, L = -Inf, D = Inf)
  expect_identical(rg$sets[[1]]$at_inf, 2L)
  # Groups with the same pwindow share a set, others do not
  rg <- .pcens_row_groups(
    c(1, 2, 3, 4),
    pwindow = c(1, 1, 2, 2), swindow = c(1, 2, 1, 1), L = -Inf, D = Inf
  )
  expect_length(rg$groups, 3)
  expect_length(rg$sets, 2)
  expect_identical(vapply(rg$groups, `[[`, integer(1), "set"), c(1L, 1L, 2L))
  expect_identical(rg$sets[[1]]$points, c(1, 2, 4))
  expect_identical(rg$sets[[1]]$pwindow, 1)
  expect_identical(rg$sets[[2]]$points, c(3, 4, 5))
})

test_that("fitdistdoublecens is unchanged with repeated delays", {
  skip_if_not_installed("fitdistrplus")
  skip_if_not_installed("withr")
  set.seed(3)
  samples <- rprimarycensored(
    500, rgamma,
    shape = 2.5, scale = 2, pwindow = 1, swindow = 1, D = 25
  )
  dat <- data.frame(
    left = floor(samples), right = floor(samples) + 1,
    pwindow = 1, D = 25
  )
  fit <- suppressMessages(fitdistdoublecens(
    dat, "gamma",
    start = list(shape = 2, scale = 2), truncation_check_multiplier = NULL
  ))
  # Reference log-likelihood computed one row at a time
  est <- as.list(fit$estimate)
  expected <- sum(vapply(seq_len(nrow(dat)), function(i) {
    log(do.call(
      dprimarycensored,
      c(list(dat$left[i], pgamma, pwindow = 1, swindow = 1, D = 25), est)
    ))
  }, numeric(1)))
  expect_equal(fit$loglik, expected, tolerance = 1e-8)
})

test_that(".dpcens gives NaN for delays outside the truncation limits", {
  one <- data.frame(swindow = 1, pwindow = 1, L = 2, D = 10)
  dens <- function(x, params) {
    .dpcens(x, params, pgamma, dunif, list(), shape = 2, rate = 1)
  }
  # Short vectors, as probed by fitdistrplus
  expect_identical(dens(c(0, 1, 5, 12), one), rep(NaN, 4))
  expect_true(all(is.finite(dens(c(2, 5, 9, 9.5), one))))
  # Full-length vectors, below L and at or above D
  three <- one[rep(1, 3), ]
  expect_identical(dens(c(1, 5, 6), three), rep(NaN, 3))
  expect_identical(dens(c(2, 5, 10), three), rep(NaN, 3))
  expect_true(all(is.finite(dens(c(2, 5, 9), three))))
})
