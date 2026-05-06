test_that("pdiscretestep returns values in [0, 1]", {
  boundaries <- 0:4
  pmf <- c(0.1, 0.3, 0.4, 0.2)
  q <- seq(-1, 5, by = 0.25)
  result <- pdiscretestep(q, boundaries, pmf)
  expect_true(all(result >= 0 & result <= 1))
})

test_that("pdiscretestep is non-decreasing", {
  boundaries <- 0:4
  pmf <- c(0.1, 0.3, 0.4, 0.2)
  q <- seq(-1, 5, by = 0.1)
  result <- pdiscretestep(q, boundaries, pmf)
  expect_true(all(diff(result) >= 0))
})

test_that("pdiscretestep equals 0 below first right edge", {
  # boundaries = 1:4, right edges at 2, 3, 4
  boundaries <- 1:4
  pmf <- c(0.2, 0.5, 0.3)
  # q below right_edges[1] = 2 should return 0
  expect_identical(
    pdiscretestep(c(0, 1, 1.999), boundaries, pmf), c(0, 0, 0)
  )
})

test_that("pdiscretestep equals 1 at and above last right edge", {
  boundaries <- 0:3
  pmf <- c(0.2, 0.5, 0.3)
  expect_equal(
    pdiscretestep(c(3, 4, 100), boundaries, pmf), c(1, 1, 1),
    tolerance = 1e-12
  )
})

test_that(
  "pdiscretestep jumps at right edges (right-continuous, mass at right edge)",
  {
    boundaries <- 0:3
    pmf <- c(0.2, 0.5, 0.3)
    # right edges: 1, 2, 3
    # Below first right edge (1): F = 0
    expect_identical(pdiscretestep(0.999, boundaries, pmf), 0)
    # At first right edge (1): F = 0.2
    expect_equal(pdiscretestep(1, boundaries, pmf), 0.2, tolerance = 1e-12)
    # Between first and second right edge: F = 0.2
    expect_equal(pdiscretestep(1.5, boundaries, pmf), 0.2, tolerance = 1e-12)
    # At second right edge (2): F = 0.7
    expect_equal(pdiscretestep(2, boundaries, pmf), 0.7, tolerance = 1e-12)
    # At last right edge (3): F = 1
    expect_equal(pdiscretestep(3, boundaries, pmf), 1, tolerance = 1e-12)
  }
)

test_that("pdiscretestep carries name attribute", {
  expect_identical(attr(pdiscretestep, "name"), "pdiscretestep")
})

test_that("pdiscretestep validates structural arguments", {
  # Length mismatch and non-monotone boundaries are structural errors.
  expect_error(
    pdiscretestep(1, boundaries = 0:3, pmf = c(0.2, 0.5)),
    "length\\(boundaries\\) must equal length\\(pmf\\) \\+ 1" # nolint: nonportable_path_linter
  )
  expect_error(
    pdiscretestep(1, boundaries = c(0, 2, 1), pmf = c(0.5, 0.5)),
    "strictly increasing"
  )
})

test_that(
  "pdiscretestep returns 0 (soft simplex penalty) for invalid pmf",
  {
    # Negative pmf -> 0 (soft penalty for use inside fitting closures)
    expect_identical(
      pdiscretestep(1, boundaries = 0:2, pmf = c(-0.1, 1.1)), 0
    )
    # Pmf not summing to 1 -> 0
    expect_identical(
      pdiscretestep(1, boundaries = 0:2, pmf = c(0.3, 0.5)), 0
    )
  }
)

test_that("ddiscretestep returns near-zero density for invalid pmf", {
  # Negative pmf -> .Machine$double.eps
  expect_identical(
    ddiscretestep(1, boundaries = 0:2, pmf = c(-0.1, 1.1)),
    .Machine$double.eps
  )
  # Pmf not summing to 1 -> .Machine$double.eps
  expect_identical(
    ddiscretestep(1, boundaries = 0:2, pmf = c(0.3, 0.5)),
    .Machine$double.eps
  )
})

test_that("ddiscretestep and pdiscretestep carry vector_param attribute", {
  expect_identical(attr(pdiscretestep, "vector_param"), "pmf")
  expect_identical(attr(ddiscretestep, "vector_param"), "pmf")
  expect_identical(attr(rdiscretestep, "vector_param"), "pmf")
})

test_that("ddiscretestep returns pmf at right boundaries, 0 elsewhere", {
  boundaries <- 0:3
  pmf <- c(0.2, 0.5, 0.3)
  # At right edges
  expect_identical(ddiscretestep(c(1, 2, 3), boundaries, pmf), pmf)
  # Not at edges
  expect_identical(
    ddiscretestep(c(0, 0.5, 1.5), boundaries, pmf), c(0, 0, 0)
  )
})

test_that("rdiscretestep samples converge to pmf", {
  set.seed(123)
  boundaries <- 0:3
  pmf <- c(0.2, 0.5, 0.3)
  n <- 50000
  samp <- rdiscretestep(n, boundaries, pmf)
  observed <- table(samp) / n
  expect_equal(as.numeric(observed), pmf, tolerance = 0.03)
})

test_that("hazards_to_pmf produces valid PMF", {
  h <- c(0.2, 0.3, 1)
  p <- hazards_to_pmf(h)
  expect_true(all(p >= 0))
  expect_equal(sum(p), 1, tolerance = 1e-10)
  expect_length(p, 3)
})

test_that("hazards_to_pmf appends trailing 1 if missing", {
  h <- c(0.2, 0.3)
  p <- hazards_to_pmf(h)
  expect_length(p, 3)
  expect_equal(sum(p), 1, tolerance = 1e-10)
})

test_that("pmf_to_hazards produces valid hazards", {
  pmf <- c(0.2, 0.3, 0.5)
  h <- pmf_to_hazards(pmf)
  expect_true(all(h >= 0 & h <= 1))
  expect_equal(h[length(h)], 1, tolerance = 1e-12)
})

test_that("hazards_to_pmf and pmf_to_hazards are mutual inverses", {
  pmf_orig <- c(0.1, 0.4, 0.3, 0.2)
  h <- pmf_to_hazards(pmf_orig)
  pmf_recovered <- hazards_to_pmf(h)
  expect_equal(pmf_orig, pmf_recovered, tolerance = 1e-10)
})

test_that("hazards_to_pmf validates inputs", {
  expect_error(
    hazards_to_pmf(c(-0.1, 1)),
    "\\[0, 1\\]"
  )
  expect_error(
    hazards_to_pmf(c(0.5, 1.2)),
    "\\[0, 1\\]"
  )
})

test_that("pmf_to_hazards validates inputs", {
  expect_error(
    pmf_to_hazards(c(-0.1, 1.1)),
    "non-negative"
  )
  expect_error(
    pmf_to_hazards(c(0.3, 0.3)),
    "sum to 1"
  )
})

test_that("pdiscretestep handles structural zero pmf entries", {
  boundaries <- 0:4
  pmf <- c(0.0, 0.5, 0.0, 0.5)
  expect_equal(sum(pmf), 1, tolerance = 1e-12)
  # Right edges: 1, 2, 3, 4
  result <- pdiscretestep(c(1, 2, 3, 4), boundaries, pmf)
  # CDF at right edges should be cumulative sums
  expect_equal(result, cumsum(pmf), tolerance = 1e-12)
})
