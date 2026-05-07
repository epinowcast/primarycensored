test_that("pdiscretehazard carries name attribute", {
  expect_identical(attr(pdiscretehazard, "name"), "pdiscretehazard")
})

test_that("pdiscretehazard round-trips with pdiscretestep via hazards_to_pmf", {
  hazards <- c(0.3, 0.5, 1)
  boundaries <- 0:3
  pmf <- hazards_to_pmf(hazards)
  q <- c(-0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5)
  expect_identical(
    pdiscretehazard(q, boundaries, hazards),
    pdiscretestep(q, boundaries, pmf)
  )
})

test_that("pdiscretehazard returns values in [0, 1]", {
  hazards <- c(0.2, 0.4, 0.6, 1)
  boundaries <- 0:4
  q <- seq(-1, 5, by = 0.25)
  result <- pdiscretehazard(q, boundaries, hazards)
  expect_true(all(result >= 0 & result <= 1))
})

test_that("pdiscretehazard is non-decreasing", {
  hazards <- c(0.2, 0.4, 0.6, 1)
  boundaries <- 0:4
  q <- seq(-1, 5, by = 0.1)
  result <- pdiscretehazard(q, boundaries, hazards)
  expect_true(all(diff(result) >= 0))
})

test_that("pdiscretehazard equals 0 below first right edge", {
  hazards <- c(0.5, 0.5, 1)
  boundaries <- 1:4
  expect_identical(
    pdiscretehazard(c(0, 1, 1.999), boundaries, hazards), c(0, 0, 0)
  )
})

test_that("pdiscretehazard equals 1 at and above last right edge", {
  hazards <- c(0.3, 0.5, 1)
  boundaries <- 0:3
  expect_equal(
    pdiscretehazard(c(3, 4, 100), boundaries, hazards), c(1, 1, 1),
    tolerance = 1e-12
  )
})

test_that("ddiscretehazard round-trips with ddiscretestep", {
  hazards <- c(0.3, 0.5, 1)
  boundaries <- 0:3
  pmf <- hazards_to_pmf(hazards)
  x <- c(0, 1, 2, 3, 1.5)
  expect_identical(
    ddiscretehazard(x, boundaries, hazards),
    ddiscretestep(x, boundaries, pmf)
  )
})

test_that("ddiscretehazard returns pmf at right boundaries, 0 elsewhere", {
  hazards <- c(0.3, 0.5, 1)
  boundaries <- 0:3
  pmf <- hazards_to_pmf(hazards)
  # Mass only at right edges 1, 2, 3
  expect_identical(
    ddiscretehazard(c(1, 2, 3), boundaries, hazards), pmf
  )
  expect_identical(
    ddiscretehazard(c(0, 0.5, 1.5), boundaries, hazards), c(0, 0, 0)
  )
})

test_that("rdiscretehazard samples converge to expected pmf", {
  set.seed(99)
  hazards <- c(0.3, 0.5, 1)
  boundaries <- 0:3
  pmf <- hazards_to_pmf(hazards)
  n <- 50000
  samp <- rdiscretehazard(n, boundaries, hazards)
  observed <- table(samp) / n
  expect_equal(as.numeric(observed), pmf, tolerance = 0.03)
})

test_that(
  "pdiscretehazard appends trailing 1 to hazards of length K-1",
  {
    hazards_short <- c(0.3, 0.5) # trailing 1 omitted
    hazards_full <- c(0.3, 0.5, 1)
    boundaries <- 0:3
    q <- seq(0, 3, by = 0.5)
    expect_identical(
      pdiscretehazard(q, boundaries, hazards_short),
      pdiscretehazard(q, boundaries, hazards_full)
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
    ddiscretehazard(c(0, 1, 2, 3), hazards = h),
    ddiscretehazard(c(0, 1, 2, 3), boundaries = 0:3, hazards = h)
  )
  set.seed(1)
  s1 <- rdiscretehazard(20, hazards = h)
  set.seed(1)
  s2 <- rdiscretehazard(20, boundaries = 0:3, hazards = h)
  expect_identical(s1, s2)
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
