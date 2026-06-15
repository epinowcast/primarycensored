test_that("pdiscretehazard carries name attribute", {
  expect_identical(attr(pdiscretehazard, "name"), "pdiscretehazard")
})

# `pdiscretehazard` is a thin wrapper around `pdiscretestep` (see roxygen).
# These identity tests assert the wrapper relationship across the p/d/r
# family in one loop; the underlying behaviour and edge cases are
# exercised in test-pdiscretestep.R rather than duplicated here.
test_that("discretehazard p/d/r equal the step equivalents under hazards_to_pmf", {
  hazards <- c(0.3, 0.5, 1)
  boundaries <- 0:3
  pmf <- hazards_to_pmf(hazards)
  q <- c(-0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5)
  identities <- list(
    p = function() list(
      pdiscretehazard(q, boundaries, hazards),
      pdiscretestep(q, boundaries, pmf)
    ),
    d = function() list(
      ddiscretehazard(c(0, 1, 1.5, 2, 3), boundaries, hazards),
      ddiscretestep(c(0, 1, 1.5, 2, 3), boundaries, pmf)
    ),
    r = function() {
      set.seed(1); s_h <- rdiscretehazard(20, boundaries, hazards)
      set.seed(1); s_p <- rdiscretestep(20, boundaries, pmf)
      list(s_h, s_p)
    }
  )
  for (nm in names(identities)) {
    pair <- identities[[nm]]()
    expect_identical(pair[[1]], pair[[2]], info = sprintf("variant: %s", nm))
  }
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

test_that("discretehazard_start validates K and the eps length", {
  expect_error(discretehazard_start(K = 0), "positive integer")
  expect_error(discretehazard_start(K = NA), "positive integer")
  expect_error(
    discretehazard_start(K = 4, eps = c(0.1, 0.2)),
    "must be a scalar or have length K - 1"
  )
  # K = 1: no innovations, just alpha and log_sigma.
  out <- discretehazard_start(K = 1)
  expect_named(out, c("alpha", "log_sigma"))
  # Vector eps of the right length flows through unchanged.
  out2 <- discretehazard_start(K = 4, eps = c(-0.1, 0, 0.1))
  expect_equal(
    unname(unlist(out2[paste0("eps_", 1:3)])),
    c(-0.1, 0, 0.1)
  )
})
