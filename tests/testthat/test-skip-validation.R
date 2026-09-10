test_that("pprimarycensored validates pdist by default and skips it when
   check = FALSE", {
  # dnorm is a density, not a CDF, so validation must reject it.
  expect_error(
    pprimarycensored(c(1, 2), stats::dnorm, pwindow = 1),
    "not a valid cumulative distribution function"
  )
  expect_no_error(
    pprimarycensored(c(1, 2), stats::dnorm, pwindow = 1, check = FALSE)
  )
})

test_that("dprimarycensored validates pdist by default and skips it when
   check = FALSE", {
  expect_error(
    dprimarycensored(c(1, 2), stats::dnorm, pwindow = 1, D = 10),
    "not a valid cumulative distribution function"
  )
  expect_no_error(
    dprimarycensored(c(1, 2), stats::dnorm, pwindow = 1, D = 10, check = FALSE)
  )
})

test_that("qprimarycensored validates pdist by default and skips it when
   check = FALSE", {
  expect_error(
    qprimarycensored(0.5, stats::dnorm, pwindow = 1),
    "not a valid cumulative distribution function"
  )
  expect_no_error(
    qprimarycensored(0.5, stats::dnorm, pwindow = 1, check = FALSE)
  )
})

test_that("check = FALSE leaves the RNG stream untouched", {
  set.seed(42)
  pprimarycensored(c(1, 2), stats::plnorm, pwindow = 1, D = 10, check = FALSE)
  after_p <- runif(3)

  set.seed(42)
  dprimarycensored(c(1, 2), stats::plnorm, pwindow = 1, D = 10, check = FALSE)
  after_d <- runif(3)

  set.seed(42)
  expected <- runif(3)

  expect_identical(after_p, expected)
  expect_identical(after_d, expected)
})

test_that("check = TRUE advances the RNG stream", {
  set.seed(42)
  pprimarycensored(c(1, 2), stats::plnorm, pwindow = 1, D = 10)
  after <- runif(3)

  set.seed(42)
  expect_false(identical(runif(3), after))
})

test_that("check = FALSE gives the same result as check = TRUE", {
  args <- list(
    c(1, 2, 3), stats::plnorm,
    pwindow = 1, D = 10, meanlog = 1, sdlog = 0.5
  )
  expect_identical(
    do.call(pprimarycensored, c(args, list(check = FALSE))),
    do.call(pprimarycensored, args)
  )
  expect_identical(
    do.call(dprimarycensored, c(args, list(check = FALSE))),
    do.call(dprimarycensored, args)
  )
})

test_that("dprimarycensored validates once, not once per internal
   pprimarycensored call", {
  n_pdist <- 0L
  n_dprimary <- 0L
  testthat::local_mocked_bindings(
    check_pdist = function(...) {
      n_pdist <<- n_pdist + 1L
      invisible(NULL)
    },
    check_dprimary = function(...) {
      n_dprimary <<- n_dprimary + 1L
      invisible(NULL)
    }
  )

  # L/D finite forces the cdf_D and cdf_L branches, so `dprimarycensored`
  # makes three internal `pprimarycensored` calls here.
  dprimarycensored(
    c(1, 2), stats::plnorm,
    pwindow = 1, L = 0.5, D = 10, meanlog = 1, sdlog = 0.5
  )

  expect_identical(n_pdist, 1L)
  expect_identical(n_dprimary, 1L)
})

test_that("fitdistdoublecens validates once, not once per likelihood
   evaluation", {
  skip_if_not_installed("fitdistrplus")
  withr::local_seed(123)

  n <- 100
  data <- data.frame(
    left = rlnorm(n, 1, 0.5),
    pwindow = 1,
    swindow = 1,
    D = 100
  )
  data$right <- data$left + data$swindow

  n_pdist <- 0L
  testthat::local_mocked_bindings(
    check_pdist = function(...) {
      n_pdist <<- n_pdist + 1L
      invisible(NULL)
    }
  )

  suppressWarnings(
    fitdistdoublecens(
      data,
      distr = "lnorm",
      start = list(meanlog = 1, sdlog = 0.5)
    )
  )

  # The optimiser calls the likelihood many times; validation must not
  # scale with it.
  expect_lte(n_pdist, 1L)
})
