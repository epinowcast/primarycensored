test_that("pprimarycensored with pwindow = 0 is the delay CDF", {
  q <- c(0, 0.5, 1, 2.5, 4, 8)
  expect_equal(
    pprimarycensored(q, pgamma, pwindow = 0, shape = 2, rate = 0.5),
    pgamma(q, shape = 2, rate = 0.5),
    tolerance = 1e-12
  )
  expect_equal(
    pprimarycensored(q, plnorm, pwindow = 0, meanlog = 1, sdlog = 0.5),
    plnorm(q, meanlog = 1, sdlog = 0.5),
    tolerance = 1e-12
  )
  expect_equal(
    pprimarycensored(q, pweibull, pwindow = 0, shape = 1.5, scale = 3),
    pweibull(q, shape = 1.5, scale = 3),
    tolerance = 1e-12
  )
  # Numeric default and a non-uniform primary distribution
  expect_equal(
    pprimarycensored(q, pnorm, pwindow = 0, mean = 3, sd = 1),
    pnorm(q, mean = 3, sd = 1),
    tolerance = 1e-12
  )
  expect_equal(
    pprimarycensored(
      q, pgamma,
      pwindow = 0, dprimary = dexpgrowth, primary_args = list(r = 0.2),
      shape = 2, rate = 0.5
    ),
    pgamma(q, shape = 2, rate = 0.5),
    tolerance = 1e-12
  )
})

test_that("pprimarycensored with pwindow = 0 handles truncation", {
  q <- c(1, 2, 3, 5)
  cdf <- function(x) pgamma(x, shape = 2, rate = 0.5)
  expect_equal(
    pprimarycensored(q, pgamma, pwindow = 0, D = 6, shape = 2, rate = 0.5),
    cdf(q) / cdf(6),
    tolerance = 1e-12
  )
  expect_equal(
    pprimarycensored(
      q, pgamma,
      pwindow = 0, L = 1, D = 6, shape = 2, rate = 0.5
    ),
    (cdf(q) - cdf(1)) / (cdf(6) - cdf(1)),
    tolerance = 1e-12
  )
})

test_that("dprimarycensored with pwindow = 0 is the delay PMF (#345)", {
  x <- 0:7
  cdf <- function(q) pgamma(q, shape = 2, rate = 0.5)
  expect_equal(
    dprimarycensored(3, pgamma, shape = 2, rate = 0.5, pwindow = 0),
    cdf(4) - cdf(3),
    tolerance = 1e-12
  )
  expect_equal(
    dprimarycensored(x, pgamma, pwindow = 0, shape = 2, rate = 0.5),
    cdf(x + 1) - cdf(x),
    tolerance = 1e-12
  )
  expect_equal(
    dprimarycensored(x, pgamma, pwindow = 0, D = 8, shape = 2, rate = 0.5),
    (cdf(x + 1) - cdf(x)) / cdf(8),
    tolerance = 1e-12
  )
  expect_equal(
    dprimarycensored(
      1:7, pgamma,
      pwindow = 0, L = 1, D = 8, shape = 2, rate = 0.5
    ),
    (cdf(2:8) - cdf(1:7)) / (cdf(8) - cdf(1)),
    tolerance = 1e-12
  )
  expect_equal(
    dprimarycensored(x, pnorm, pwindow = 0, swindow = 0.5, mean = 3, sd = 1),
    pnorm(x + 0.5, 3, 1) - pnorm(x, 3, 1),
    tolerance = 1e-12
  )
})

test_that("pcens objects support pwindow = 0", {
  obj <- new_pcens(pgamma, dunif, list(), shape = 2, rate = 0.5)
  q <- c(0.5, 2, 4)
  expect_equal(
    pcens_cdf(obj, q, pwindow = 0), pgamma(q, 2, 0.5),
    tolerance = 1e-12
  )
  expect_equal(
    pcens_cdf(obj, q, pwindow = 0, use_numeric = TRUE), pgamma(q, 2, 0.5),
    tolerance = 1e-12
  )
  expect_equal(
    pcens_pmf(obj, 0:4, pwindow = 0, D = 10),
    (pgamma(1:5, 2, 0.5) - pgamma(0:4, 2, 0.5)) / pgamma(10, 2, 0.5),
    tolerance = 1e-12
  )
  obj <- update(obj, shape = 3)
  expect_equal(
    pcens_cdf(obj, q, pwindow = 0), pgamma(q, 3, 0.5),
    tolerance = 1e-12
  )
})

test_that("qprimarycensored works with pwindow = 0", {
  expect_equal(
    qprimarycensored(c(0.25, 0.5), pgamma, pwindow = 0, shape = 2, rate = 1),
    qgamma(c(0.25, 0.5), shape = 2, rate = 1),
    tolerance = 1e-4
  )
})

test_that("check_dprimary accepts pwindow = 0", {
  expect_null(check_dprimary(dunif, pwindow = 0))
  expect_null(check_dprimary(dexpgrowth, pwindow = 0, list(r = 0.2)))
})

test_that("dprimarycensored with swindow = 0 returns a density", {
  x <- c(0.5, 1, 2.5, 4)
  # Exact primary and secondary events: the delay density
  expect_equal(
    dprimarycensored(
      x, pgamma,
      pwindow = 0, swindow = 0, shape = 2, rate = 0.5
    ),
    dgamma(x, shape = 2, rate = 0.5),
    tolerance = 1e-12
  )
  # Uniform primary: the primary censored density
  expect_equal(
    dprimarycensored(
      x, pgamma,
      pwindow = 1, swindow = 0, shape = 2, rate = 0.5
    ),
    pgamma(x, 2, 0.5) - pgamma(x - 1, 2, 0.5),
    tolerance = 1e-10
  )
  expect_equal(
    dprimarycensored(
      x, plnorm,
      pwindow = 2, swindow = 0, meanlog = 1, sdlog = 1
    ),
    (plnorm(x, 1, 1) - plnorm(x - 2, 1, 1)) / 2,
    tolerance = 1e-10
  )
  # Non-uniform primary: integrate the delay density over the primary
  expected <- vapply(x, function(xi) {
    integrate(
      function(p) dgamma(xi - p, 2, 0.5) * dexpgrowth(p, 0, 1, r = 0.5),
      0, 1
    )$value
  }, numeric(1))
  expect_equal(
    dprimarycensored(
      x, pgamma,
      pwindow = 1, swindow = 0, dprimary = dexpgrowth,
      primary_args = list(r = 0.5), shape = 2, rate = 0.5
    ),
    expected,
    tolerance = 1e-6
  )
})

test_that("swindow = 0 gives the limit of the pmf over swindow", {
  h <- 1e-6
  x <- c(0.5, 2, 3.5)
  dens <- dprimarycensored(
    x, pweibull,
    pwindow = 1, swindow = 0, D = 10, shape = 1.5, scale = 2
  )
  approx <- dprimarycensored(
    x, pweibull,
    pwindow = 1, swindow = h, D = 10, shape = 1.5, scale = 2
  ) / h
  expect_equal(dens, approx, tolerance = 1e-5)
})

test_that("dprimarycensored with swindow = 0 handles truncation and log", {
  x <- c(1.5, 2, 3.5)
  cdf <- function(q) pgamma(q, 2, 0.5)
  expect_equal(
    dprimarycensored(
      x, pgamma,
      pwindow = 0, swindow = 0, L = 1, D = 6, shape = 2, rate = 0.5
    ),
    dgamma(x, 2, 0.5) / (cdf(6) - cdf(1)),
    tolerance = 1e-12
  )
  expect_equal(
    dprimarycensored(
      x, pgamma,
      pwindow = 0, swindow = 0, D = 6, log = TRUE, shape = 2, rate = 0.5
    ),
    dgamma(x, 2, 0.5, log = TRUE) - log(cdf(6)),
    tolerance = 1e-12
  )
})

test_that("pcens_pmf mixes densities and probabilities by swindow", {
  obj <- new_pcens(pgamma, dunif, list(), shape = 2, rate = 0.5)
  x <- c(1, 2, 3)
  sw <- c(0, 1, 0)
  expect_equal(
    pcens_pmf(obj, x, pwindow = 0, swindow = sw),
    c(dgamma(1, 2, 0.5), pgamma(3, 2, 0.5) - pgamma(2, 2, 0.5),
      dgamma(3, 2, 0.5)),
    tolerance = 1e-12
  )
})

test_that("swindow = 0 errors when no delay density can be found", {
  pmystery <- add_name_attribute(
    function(q, rate) pexp(q, rate), "pmystery"
  )
  expect_error(
    dprimarycensored(1, pmystery, pwindow = 0, swindow = 0, rate = 1),
    "density"
  )
  # A uniform primary with a positive window only needs the CDF
  expect_equal(
    dprimarycensored(1, pmystery, pwindow = 1, swindow = 0, rate = 1),
    pexp(1) - pexp(0),
    tolerance = 1e-8
  )
})

# ---- fitdistdoublecens ----------------------------------------------------

fit_base_data <- function() {
  data.frame(
    left = c(2, 3, 4, 5, 3, 6, 2), right = c(3, 4, 5, 6, 4, 7, 3),
    pwindow = 1, D = Inf
  )
}

test_that("fitdistdoublecens fits rows with pwindow = 0 (#345)", {
  skip_if_not_installed("fitdistrplus")
  p0 <- fit_base_data()
  p0$pwindow[c(2, 5)] <- 0
  fit <- fitdistdoublecens(
    p0, "gamma",
    start = list(shape = 2, rate = 0.5)
  )
  expect_true(all(is.finite(fit$estimate)))
  expect_false(isTRUE(all.equal(unname(fit$estimate), c(2, 0.5))))
})

test_that("fitdistdoublecens fits rows with left == right (#345)", {
  skip_if_not_installed("fitdistrplus")
  s0 <- fit_base_data()
  s0$right[c(2, 5)] <- s0$left[c(2, 5)]
  fit <- fitdistdoublecens(
    s0, "gamma",
    start = list(shape = 2, rate = 0.5)
  )
  expect_true(all(is.finite(fit$estimate)))
  expect_false(isTRUE(all.equal(unname(fit$estimate), c(2, 0.5))))
})

# Simulate delays of four observation types: interval or exact primary
# event, crossed with interval or exact secondary event.
simulate_mixed <- function(n, rdist, ...) {
  delay <- rdist(n, ...)
  pexact <- rep_len(c(TRUE, FALSE), n)
  sexact <- rep_len(c(TRUE, TRUE, FALSE, FALSE), n)
  primary <- ifelse(pexact, 0, runif(n))
  secondary <- primary + delay
  left <- ifelse(sexact, secondary, floor(secondary))
  data.frame(
    left = left,
    right = ifelse(sexact, left, left + 1),
    pwindow = as.numeric(!pexact),
    D = Inf
  )
}

# Hand-written likelihood for mixed data with a uniform primary event.
mixed_nll <- function(par, data, pdist, ddist) {
  a <- exp(par[1])
  b <- exp(par[2])
  cdf <- function(q, pw) {
    if (pw == 0) {
      return(pdist(q, a, b))
    }
    integrate(function(p) pdist(q - p, a, b), 0, pw)$value / pw
  }
  ll <- vapply(seq_len(nrow(data)), function(i) {
    l <- data$left[i]
    r <- data$right[i]
    pw <- data$pwindow[i]
    if (l == r) {
      dens <- if (pw == 0) {
        ddist(l, a, b)
      } else {
        (pdist(l, a, b) - pdist(l - pw, a, b)) / pw
      }
      return(log(dens))
    }
    log(cdf(r, pw) - cdf(l, pw))
  }, numeric(1))
  -sum(ll)
}

test_that("fitdistdoublecens recovers parameters from mixed-type data", {
  skip_if_not_installed("fitdistrplus")
  set.seed(345)
  data <- simulate_mixed(400, rgamma, shape = 4, rate = 1)
  fit <- fitdistdoublecens(
    data, "gamma",
    start = list(shape = 2, rate = 0.5)
  )
  expect_equal(
    unname(fit$estimate), c(4, 1),
    tolerance = 0.2
  )
  hand <- optim(
    log(c(2, 0.5)), mixed_nll,
    data = data, pdist = pgamma, ddist = dgamma
  )
  expect_equal(
    unname(fit$estimate), exp(hand$par),
    tolerance = 1e-2
  )
})

test_that("fitdistdoublecens matches coarseDataTools::dic.fit", {
  skip_if_not_installed("fitdistrplus")
  skip_if_not_installed("coarseDataTools")
  set.seed(262)
  data <- simulate_mixed(200, rlnorm, meanlog = 1.2, sdlog = 0.5)
  # dic.fit() approximates doubly interval censored rows whose windows
  # overlap by a single interval, so drop them.
  data <- data[!(data$pwindow > 0 & data$right > data$left & data$left < 1), ]
  fit <- fitdistdoublecens(
    data, "lnorm",
    start = list(meanlog = 1, sdlog = 1)
  )
  # coarseDataTools columns are absolute times with the primary window
  # starting at 0: EL = 0, ER = pwindow, SL = left, SR = right.
  # Its types are 0 (doubly interval censored), 1 (single interval
  # censored) and 2 (exact).
  pexact <- data$pwindow == 0
  sexact <- data$left == data$right
  cdt_data <- data.frame(
    EL = 0, ER = data$pwindow, SL = data$left, SR = data$right,
    type = 2 * (pexact & sexact) + (xor(pexact, sexact))
  )
  cdt_fit <- suppressMessages(suppressWarnings(
    coarseDataTools::dic.fit(cdt_data, dist = "L")
  ))
  expect_equal(
    unname(fit$estimate), unname(cdt_fit@ests[1:2, 1]),
    tolerance = 1e-2
  )
})
