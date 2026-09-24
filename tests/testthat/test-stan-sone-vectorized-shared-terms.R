skip_on_cran()

# With a uniform primary and an integer pwindow,
# primarycensored_sone_lpmf_vectorized() computes the analytical uniform
# primary terms once per integer node and shares them. These tests check the
# shared path against the per-node path.

shared_terms_dists <- list(
  list(dist_id = 1L, params = list(
    c(1.5, 0.5), c(0.2, 1.2), c(-1, 0.3), c(3.5, 0.05), c(4, 2)
  )),
  list(dist_id = 2L, params = list(
    c(2, 0.5), c(0.5, 0.2), c(20, 4), c(1.2, 0.02), c(50, 0.5)
  )),
  list(dist_id = 3L, params = list(
    c(1.5, 5), c(0.7, 2), c(3, 30), c(5, 0.8)
  )),
  list(dist_id = 5L, params = list(c(1.5, 3, 2), c(0.8, 1, 0.5)))
)

test_that("check_for_uniform_terms selects the analytical uniform delays", {
  for (dist_id in c(1L, 2L, 3L, 5L)) {
    expect_identical(check_for_uniform_terms(dist_id, 1L), 1L)
    expect_identical(check_for_uniform_terms(dist_id, 2L), 0L)
  }
  for (dist_id in c(4L, 18L, 26L, 27L, 28L)) {
    expect_identical(check_for_uniform_terms(dist_id, 1L), 0L)
  }
})

test_that("uniform primary terms are -Inf for t <= 0", {
  for (dist in shared_terms_dists) {
    for (t in c(0, -0.5, -3)) {
      expect_identical(
        primarycensored_uniform_terms(t, dist$dist_id, dist$params[[1]]),
        c(-Inf, -Inf)
      )
    }
  }
})

test_that("primarycensored_uniform_lcdf_from_terms is -Inf when all terms
  underflow", {
  expect_identical(
    primarycensored_uniform_lcdf_from_terms(c(-Inf, -Inf), c(-Inf, -Inf), 2),
    -Inf
  )
})

test_that("shared node log CDFs match per-node log CDFs", {
  n <- 41L
  for (dist in shared_terms_dists) {
    for (params in dist$params) {
      for (pwindow in c(1L, 2L, 3L, 7L)) {
        for (start in c(1L, 4L, 10L)) {
          info <- paste(
            "dist", dist$dist_id, "params", toString(params),
            "pwindow", pwindow, "start", start
          )
          shared <- primarycensored_uniform_node_log_cdfs(
            start, n, dist$dist_id, params, pwindow
          )
          node <- primarycensored_node_log_cdfs(
            start, n, dist$dist_id, params, pwindow, 1L, numeric(0)
          )
          expect_length(shared, n)
          expect_identical(shared[start:n], node[start:n], info = info)
        }
      }
    }
  }
})

test_that("primarycensored_sone_lpmf_vectorized with shared terms matches
  primarycensored_lpmf", {
  dists <- list(
    list(dist_id = 1L, params = c(1.5, 0.5)),
    list(dist_id = 1L, params = c(0.2, 1.2)),
    list(dist_id = 2L, params = c(2, 0.5)),
    list(dist_id = 2L, params = c(20, 4)),
    list(dist_id = 3L, params = c(1.5, 5)),
    list(dist_id = 3L, params = c(0.7, 2)),
    list(dist_id = 5L, params = c(1.5, 3, 2))
  )
  settings <- list(
    list(max_delay = 0, L = 0, D = 1),
    list(max_delay = 1, L = 0, D = Inf),
    list(max_delay = 20, L = 0, D = 21),
    list(max_delay = 20, L = -Inf, D = 21),
    list(max_delay = 20, L = 3, D = 21),
    list(max_delay = 20, L = 2.5, D = 25.5),
    list(max_delay = 15, L = 0, D = 30),
    list(max_delay = 20, L = 4, D = Inf)
  )
  for (dist in dists) {
    for (s in settings) {
      for (pwindow in c(1, 2, 3, 7)) {
        info <- paste(
          "dist", dist$dist_id, "params", toString(dist$params),
          "max_delay", s$max_delay, "L", s$L, "D", s$D, "pwindow", pwindow
        )
        vectorised <- primarycensored_sone_lpmf_vectorized(
          s$max_delay, s$L, s$D, dist$dist_id, dist$params, pwindow,
          1L, numeric(0)
        )
        # A bin that L cuts is left out: primarycensored_lpmf() only takes
        # integer lower bounds
        delays <- 0:s$max_delay
        full <- delays >= s$L
        below <- delays + 1 <= s$L
        per_delay <- vapply(delays[full], function(d) {
          primarycensored_lpmf(
            d, dist$dist_id, dist$params, pwindow, d + 1, s$L, s$D,
            1L, numeric(0)
          )
        }, numeric(1))
        expect_equal(
          vectorised[full], per_delay, tolerance = 1e-10, info = info
        )
        expect_identical(vectorised[below], rep(-Inf, sum(below)))
      }
    }
  }
})
