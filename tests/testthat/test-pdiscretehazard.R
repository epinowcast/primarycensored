test_that("pdiscretehazard carries name attribute", {
  expect_identical(attr(pdiscretehazard, "name"), "pdiscretehazard")
})

# `pdiscretehazard` is a thin wrapper around `pdiscretestep` (see roxygen).
# These identity tests assert the wrapper relationship across the family;
# the underlying behaviour and edge cases are exercised in
# test-pdiscretestep.R rather than duplicated here.
test_that("pdiscretehazard equals pdiscretestep with hazards_to_pmf", {
  hazards <- c(0.3, 0.5, 1)
  boundaries <- 0:3
  pmf <- hazards_to_pmf(hazards)
  q <- c(-0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5)
  expect_identical(
    pdiscretehazard(q, boundaries, hazards),
    pdiscretestep(q, boundaries, pmf)
  )
})

test_that("ddiscretehazard equals ddiscretestep with hazards_to_pmf", {
  hazards <- c(0.3, 0.5, 1)
  boundaries <- 0:3
  pmf <- hazards_to_pmf(hazards)
  x <- c(0, 1, 1.5, 2, 3)
  expect_identical(
    ddiscretehazard(x, boundaries, hazards),
    ddiscretestep(x, boundaries, pmf)
  )
})

test_that("rdiscretehazard equals rdiscretestep with hazards_to_pmf", {
  hazards <- c(0.3, 0.5, 1)
  boundaries <- 0:3
  pmf <- hazards_to_pmf(hazards)
  set.seed(1)
  s_h <- rdiscretehazard(20, boundaries, hazards)
  set.seed(1)
  s_p <- rdiscretestep(20, boundaries, pmf)
  expect_identical(s_h, s_p)
})

test_that(
  "pdiscretehazard appends trailing 1 to hazards of length K-1",
  {
    boundaries <- 0:3
    q <- seq(0, 3, by = 0.5)
    expect_identical(
      pdiscretehazard(q, boundaries, hazards = c(0.3, 0.5)),
      pdiscretehazard(q, boundaries, hazards = c(0.3, 0.5, 1))
    )
  }
)

test_that("pdiscretehazard propagates validation errors from hazards_to_pmf", {
  expect_error(
    pdiscretehazard(1, boundaries = 0:2, hazards = c(-0.1, 1)),
    "\\[0, 1\\]"
  )
})

test_that("default boundaries are 0:length(hazards) when omitted", {
  h <- c(0.3, 0.5, 1)
  expect_identical(
    pdiscretehazard(c(0.5, 1, 2, 3), hazards = h),
    pdiscretehazard(c(0.5, 1, 2, 3), boundaries = 0:3, hazards = h)
  )
  expect_identical(
    ddiscretehazard(c(0.5, 1, 2, 3), hazards = h),
    ddiscretehazard(c(0.5, 1, 2, 3), boundaries = 0:3, hazards = h)
  )
  set.seed(1)
  s1 <- rdiscretehazard(20, hazards = h)
  set.seed(1)
  s2 <- rdiscretehazard(20, boundaries = 0:3, hazards = h)
  expect_identical(s1, s2)
})

test_that(".resolve_hazard_prior ignores unknown components", {
  out <- primarycensored:::.resolve_hazard_prior(
    list(unknown = list(mean = 99))
  )
  defaults <- primarycensored:::.resolve_hazard_prior(NULL)
  expect_identical(out, defaults)
})

test_that("user prior partially overrides defaults in discretehazard", {
  set.seed(404)
  K <- 4L
  true_pmf <- c(0.4, 0.3, 0.2, 0.1)
  D_val <- K + 1L
  n <- 200
  samples <- rprimarycensored(
    n, rdiscretestep,
    boundaries = 0:K, pmf = true_pmf,
    pwindow = 1, swindow = 1, D = D_val
  )
  haz_data <- data.frame(
    left    = samples,
    right   = samples + 1,
    pwindow = rep(1, n),
    D       = rep(D_val, n)
  )
  start <- c(
    list(alpha = -1, log_sigma = log(0.5)),
    as.list(stats::setNames(
      rep(0, K - 1L),
      paste0("eps_", seq_len(K - 1L))
    ))
  )
  fit <- fitdistdoublecens(
    haz_data,
    distr      = "discretehazard",
    start      = start,
    boundaries = 0:K,
    prior      = list(alpha = list(mean = -2))
  )
  expect_s3_class(fit, "fitdist")
  expect_false(is.na(fit$loglik))
})
