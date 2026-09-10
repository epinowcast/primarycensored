# A function with a CDF-compatible signature that is not a CDF: it returns
# values outside [0, 1], so `check_pdist()` must reject it.
bad_pdist <- function(q, ...) rep(2, length(q))

# Counts for the mocked validation helpers. These live in an environment so
# the stubs mutate it rather than rebinding with `<<-`.
new_counter <- function() {
  counts <- new.env(parent = emptyenv())
  counts$pdist <- 0L
  counts$dprimary <- 0L
  counts
}

test_that("pprimarycensored validates pdist by default and skips it when
   check = FALSE", {
  expect_error(
    pprimarycensored(c(1, 2), bad_pdist, pwindow = 1, D = 10),
    "not a valid cumulative distribution function"
  )
  expect_no_error(
    pprimarycensored(c(1, 2), bad_pdist, pwindow = 1, D = 10, check = FALSE)
  )
})

test_that("dprimarycensored validates pdist by default and skips it when
   check = FALSE", {
  expect_error(
    dprimarycensored(c(1, 2), bad_pdist, pwindow = 1, D = 10),
    "not a valid cumulative distribution function"
  )
  expect_no_error(
    dprimarycensored(c(1, 2), bad_pdist, pwindow = 1, D = 10, check = FALSE)
  )
})

test_that("qprimarycensored validates pdist by default and skips it when
   check = FALSE", {
  expect_error(
    qprimarycensored(0.5, bad_pdist, pwindow = 1, D = 10),
    "not a valid cumulative distribution function"
  )
  # With `check = FALSE` the validation error is skipped. An invalid CDF then
  # fails further downstream instead, which is the caller's responsibility.
  err <- tryCatch(
    qprimarycensored(0.5, bad_pdist, pwindow = 1, D = 10, check = FALSE),
    error = conditionMessage
  )
  expect_no_match(err, "not a valid cumulative distribution function")
})

test_that("check = FALSE leaves the RNG stream untouched", {
  set.seed(42)
  pprimarycensored(
    c(1, 2), stats::plnorm,
    pwindow = 1, D = 10, meanlog = 1, sdlog = 0.5, check = FALSE
  )
  after_p <- runif(3)

  set.seed(42)
  dprimarycensored(
    c(1, 2), stats::plnorm,
    pwindow = 1, D = 10, meanlog = 1, sdlog = 0.5, check = FALSE
  )
  after_d <- runif(3)

  set.seed(42)
  expected <- runif(3)

  expect_identical(after_p, expected)
  expect_identical(after_d, expected)
})

test_that("check = TRUE advances the RNG stream", {
  set.seed(42)
  pprimarycensored(
    c(1, 2), stats::plnorm,
    pwindow = 1, D = 10, meanlog = 1, sdlog = 0.5
  )
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
  counts <- new_counter()
  testthat::local_mocked_bindings(
    check_pdist = function(...) {
      counts$pdist <- counts$pdist + 1L
      invisible(NULL)
    },
    check_dprimary = function(...) {
      counts$dprimary <- counts$dprimary + 1L
      invisible(NULL)
    }
  )

  # Finite L and D force the cdf_D and cdf_L branches, so `dprimarycensored`
  # makes three internal `pprimarycensored` calls here.
  dprimarycensored(
    c(1, 2), stats::plnorm,
    pwindow = 1, L = 0.5, D = 10, meanlog = 1, sdlog = 0.5
  )

  expect_identical(counts$pdist, 1L)
  expect_identical(counts$dprimary, 1L)
})

test_that(".ppcens and .dpcens validate once per call, not once per
   observation", {
  counts <- new_counter()
  testthat::local_mocked_bindings(
    check_pdist = function(...) {
      counts$pdist <- counts$pdist + 1L
      invisible(NULL)
    },
    check_dprimary = function(...) {
      counts$dprimary <- counts$dprimary + 1L
      invisible(NULL)
    }
  )

  # Varying pwindow and D exercise the grouped `.dpcens` branch and the
  # per-observation `mapply` in `.ppcens`.
  params <- data.frame(
    swindow = 1,
    pwindow = c(1, 1, 2, 2),
    L = -Inf,
    D = c(10, 10, 20, 20)
  )

  .ppcens(
    q = c(1, 2, 3, 4), params = params, pdist = stats::plnorm,
    dprimary = stats::dunif, dprimary_args = list(),
    meanlog = 1, sdlog = 0.5
  )
  expect_identical(counts$pdist, 1L)
  # One `check_dprimary` per unique pwindow, not one per observation.
  expect_identical(counts$dprimary, 2L)

  counts$pdist <- 0L
  counts$dprimary <- 0L
  .dpcens(
    x = c(1, 2, 3, 4), params = params, pdist = stats::plnorm,
    dprimary = stats::dunif, dprimary_args = list(),
    meanlog = 1, sdlog = 0.5
  )
  expect_identical(counts$pdist, 1L)
  expect_identical(counts$dprimary, 2L)
})

test_that(".ppcens and .dpcens skip validation entirely with check = FALSE", {
  counts <- new_counter()
  testthat::local_mocked_bindings(
    check_pdist = function(...) {
      counts$pdist <- counts$pdist + 1L
      invisible(NULL)
    },
    check_dprimary = function(...) {
      counts$dprimary <- counts$dprimary + 1L
      invisible(NULL)
    }
  )
  params <- data.frame(swindow = 1, pwindow = 1, L = -Inf, D = 10)

  .ppcens(
    q = 1, params = params, pdist = stats::plnorm,
    dprimary = stats::dunif, dprimary_args = list(),
    check = FALSE, meanlog = 1, sdlog = 0.5
  )
  .dpcens(
    x = 1, params = params, pdist = stats::plnorm,
    dprimary = stats::dunif, dprimary_args = list(),
    check = FALSE, meanlog = 1, sdlog = 0.5
  )
  expect_identical(counts$pdist + counts$dprimary, 0L)
})

test_that("fitdistdoublecens validates once per fit, not once per likelihood
   evaluation", {
  skip_if_not_installed("fitdistrplus")
  withr::local_seed(123)

  n <- 100
  data <- data.frame(
    left = rlnorm(n, 1, 0.5),
    # Varying pwindow and D so the grouped and mapply paths are both used.
    pwindow = rep_len(c(1, 2), n),
    swindow = 1,
    D = rep_len(c(50, 100), n)
  )
  data$right <- data$left + data$swindow

  counts <- new_counter()
  testthat::local_mocked_bindings(
    check_pdist = function(...) {
      counts$pdist <- counts$pdist + 1L
      invisible(NULL)
    }
  )

  suppressMessages(suppressWarnings(
    fitdistdoublecens(
      data,
      distr = "lnorm",
      start = list(meanlog = 1, sdlog = 0.5)
    )
  ))

  # The optimiser calls the likelihood many times; validation must not
  # scale with it.
  expect_identical(counts$pdist, 1L)
})

test_that("fitdistdoublecens skips validation with check = FALSE", {
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

  counts <- new_counter()
  testthat::local_mocked_bindings(
    check_pdist = function(...) {
      counts$pdist <- counts$pdist + 1L
      invisible(NULL)
    }
  )

  suppressMessages(suppressWarnings(
    fitdistdoublecens(
      data,
      distr = "lnorm",
      start = list(meanlog = 1, sdlog = 0.5),
      check = FALSE
    )
  ))

  expect_identical(counts$pdist, 0L)
})
