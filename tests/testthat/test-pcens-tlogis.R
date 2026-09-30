families <- exptilt_families()
normals <- families[vapply(families, function(f) !f$positive, logical(1))]
positives <- families[vapply(families, function(f) f$positive, logical(1))]

# The grid of the acceptance criteria plus a midpoint far outside the window
# and scales from sharp to nearly uniform
locations <- c(-0.5, 0, 0.5, 1, 1.5)
scales <- c(0.1, 0.2, 1)
windows <- c(1, 2)

# Returns the q used against a reference, including the tails, q below the
# window and q close to zero
tlogis_quantiles <- function(family, pwindow) {
  sort(c(
    1e-6, 1e-3, 0.3 * pwindow, pwindow - 1e-3, pwindow, pwindow + 1e-3,
    2, 3, 6, 12, 25,
    if (!family$positive) c(-30, -10, -3, -0.5)
  ))
}

test_that("truncated logistic primaries dispatch to analytic methods", {
  expect_s3_class(
    tlogis_object(families[[1]], 0.5, 0.2), "pcens_pexp_dtlogis"
  )
  expect_s3_class(
    tlogis_object(families[[3]], 0.5, 0.2), "pcens_pgamma_dtlogis"
  )
  expect_s3_class(
    tlogis_object(families[[6]], 0.5, 0.2), "pcens_pnorm_dtlogis"
  )
  for (cls in c(
    "pcens_pexp_dtlogis", "pcens_pgamma_dtlogis", "pcens_pnorm_dtlogis"
  )) {
    expect_false(
      is.null(utils::getS3method("pcens_cdf", cls, optional = TRUE)),
      info = cls
    )
  }
})

test_that("the truncation rule bounds the error of the averaged series", {
  # The series 1 / (1 + x) = sum (-x)^n is summed with weights W_n that
  # are 1 up to n0 and taper to 0 at n0 + M. The error must be below the
  # tolerance for every x in the range the rule was given.
  ranges <- list(
    c(0, 1), c(0.2, 1), c(0.9, 1), c(0, 0.6), c(0.1, 0.3), c(1e-3, 1e-2),
    c(0, 1e-40)
  )
  for (tol in c(1e-8, 1e-12)) {
    for (range in ranges) {
      terms <- .tlogis_series_terms(log(range[1]), log(range[2]), log(tol))
      expect_false(is.null(terms), info = toString(range))
      weights <- .tlogis_weights(terms[["n0"]], terms[["M"]])
      expect_length(weights, terms[["n0"]] + terms[["M"]])
      x <- seq(range[1], range[2], length.out = 200)
      n <- seq_along(weights) - 1
      approx <- vapply(
        x, function(xx) sum((-xx)^n * weights), numeric(1)
      )
      expect_lt(
        max(abs(approx - 1 / (1 + x))), tol,
        label = sprintf("tol %g, range %s", tol, toString(range))
      )
    }
  }
})

test_that("the truncation rule gives no rule where too many terms are needed", {
  expect_null(.tlogis_series_terms(log(1e-300), 0, log(1e-300), max_terms = 20))
  expect_null(.tlogis_series_terms(-Inf, 0, -Inf))
})

test_that("the weights taper from 1 to 0", {
  weights <- .tlogis_weights(3L, 5L)
  expect_identical(weights[1:3], rep(1, 3))
  expect_identical(weights[4:8], pbinom(0:4, 5, 0.5, lower.tail = FALSE))
  expect_true(all(diff(weights) <= 0))
  expect_identical(.tlogis_weights(4L, 0L), rep(1, 4))
})

test_that("the analytic CDF of a normal delay matches a reference integral", {
  for (family in normals) {
    for (pwindow in windows) {
      q <- tlogis_quantiles(family, pwindow)
      for (location in locations) {
        for (scale in scales) {
          obj <- tlogis_object(family, location, scale)
          expected <- tlogis_reference(
            q, pwindow, location, scale, exptilt_cdf(family), FALSE
          )
          actual <- pcens_cdf(obj, q, pwindow)
          expect_lt(
            max_rel_diff(actual, expected), 1e-8,
            label = tlogis_label(family, pwindow, location, scale)
          )
        }
      }
    }
  }
})

test_that("the analytic CDF covers midpoints far outside the window and
  extreme scales", {
  family <- normals[[1]]
  cases <- list(
    c(pwindow = 0.5, location = -6, scale = 0.02),
    c(pwindow = 3, location = -2, scale = 5),
    c(pwindow = 3, location = 12, scale = 0.05),
    c(pwindow = 3, location = 8, scale = 2),
    c(pwindow = 0.5, location = 0.25, scale = 50),
    c(pwindow = 3, location = 1.5, scale = 0.02),
    c(pwindow = 1, location = 1e-9, scale = 0.3),
    c(pwindow = 1, location = 1 - 1e-9, scale = 0.3)
  )
  for (cs in cases) {
    q <- tlogis_quantiles(family, cs[["pwindow"]])
    obj <- tlogis_object(family, cs[["location"]], cs[["scale"]])
    expected <- tlogis_reference(
      q, cs[["pwindow"]], cs[["location"]], cs[["scale"]],
      exptilt_cdf(family), FALSE
    )
    expect_lt(
      max_rel_diff(pcens_cdf(obj, q, cs[["pwindow"]]), expected), 1e-8,
      label = toString(cs)
    )
  }
})

test_that("gamma and exponential delays use the analytic CDF when the midpoint
  is after the window", {
  # Only negative tilts are needed, which always exist
  for (family in positives) {
    for (pwindow in windows) {
      q <- tlogis_quantiles(family, pwindow)
      for (location in c(pwindow + 1e-3, pwindow + 0.5, pwindow + 5)) {
        for (scale in c(0.1, 0.5, 2)) {
          obj <- tlogis_object(family, location, scale)
          plan <- .tlogis_plan(obj, pwindow)
          expect_false(is.null(plan), label = family$label)
          expected <- tlogis_reference(
            q, pwindow, location, scale, exptilt_cdf(family), TRUE
          )
          expect_lt(
            max_rel_diff(pcens_cdf(obj, q, pwindow), expected), 1e-8,
            label = tlogis_label(family, pwindow, location, scale)
          )
        }
      }
    }
  }
})

test_that("gamma delays use the analytic CDF with positive tilts when the
  rate allows them", {
  # The positive tilts reach (n0 + M) / scale, which must be below the rate
  family <- list(
    label = "gamma shape 2.5, rate 1000", pdist = pgamma,
    args = list(shape = 2.5, rate = 1000), rate = 1000, positive = TRUE
  )
  for (pwindow in windows) {
    q <- c(1e-4, 1e-3, 3e-3, 0.01, 0.1, 0.5, 1, 2)
    for (location in c(-3, -0.5, 0, 0.5, 1)) {
      for (scale in c(0.5, 1, 5)) {
        obj <- tlogis_object(family, location, scale)
        expect_false(is.null(.tlogis_plan(obj, pwindow)))
        expected <- tlogis_reference(
          q, pwindow, location, scale, exptilt_cdf(family), TRUE
        )
        expect_lt(
          max_rel_diff(pcens_cdf(obj, q, pwindow), expected), 1e-8,
          label = tlogis_label(family, pwindow, location, scale)
        )
      }
    }
  }
})

test_that("inadmissible positive tilts use the numerical method", {
  # rate 0.4 is far below the tilts of the series for these scales
  for (family in positives[c(2, 4)]) {
    obj <- tlogis_object(family, 0.5, 0.2)
    expect_null(.tlogis_plan(obj, 2))
    q <- c(0.5, 2, 5, 10)
    expect_identical(
      pcens_cdf(obj, q, 2),
      pcens_cdf.default(obj, q, 2),
      info = family$label
    )
    # The fallback is the correct CDF
    expected <- tlogis_reference(
      q, 2, 0.5, 0.2, exptilt_cdf(family), TRUE
    )
    expect_lt(
      max_rel_diff(pcens_cdf(obj, q, 2), expected), 1e-4,
      label = family$label
    )
  }
})

test_that("the analytic CDF agrees with use_numeric = TRUE", {
  q <- c(0.05, 0.5, 1.5, 3, 6, 12, 20)
  for (family in families) {
    for (pwindow in windows) {
      for (location in locations) {
        for (scale in scales) {
          obj <- tlogis_object(family, location, scale)
          analytic <- pcens_cdf(obj, q, pwindow)
          numeric <- pcens_cdf(obj, q, pwindow, use_numeric = TRUE)
          expect_equal(
            analytic, numeric,
            tolerance = 1e-6,
            info = tlogis_label(family, pwindow, location, scale)
          )
        }
      }
    }
  }
})

test_that("use_numeric = TRUE uses the default method", {
  obj <- tlogis_object(families[[6]], 0.5, 0.2)
  expect_identical(
    pcens_cdf(obj, c(1, 4), 2, use_numeric = TRUE),
    pcens_cdf.default(obj, c(1, 4), 2)
  )
})

test_that("the analytic CDF is continuous in the location across 0 and the
  window", {
  family <- normals[[1]]
  q <- c(-2, 0.3, 1.5, 4, 9)
  for (pwindow in windows) {
    for (scale in scales) {
      for (edge in c(0, pwindow)) {
        below <- pcens_cdf(
          tlogis_object(family, edge - 1e-7, scale), q, pwindow
        )
        at <- pcens_cdf(tlogis_object(family, edge, scale), q, pwindow)
        above <- pcens_cdf(
          tlogis_object(family, edge + 1e-7, scale), q, pwindow
        )
        expect_equal(below, at, tolerance = 1e-6)
        expect_equal(above, at, tolerance = 1e-6)
      }
    }
  }
})

test_that("the analytic CDF is monotone, bounded and handles special q", {
  family <- normals[[1]]
  obj <- tlogis_object(family, 0.5, 0.2)
  q <- seq(-40, 40, by = 0.25)
  cdf <- pcens_cdf(obj, q, 2)
  expect_true(all(cdf >= 0 & cdf <= 1))
  expect_true(all(diff(cdf) >= -1e-12))
  expect_identical(pcens_cdf(obj, c(-Inf, Inf), 2), c(0, 1))
  expect_identical(pcens_cdf(obj, numeric(0), 2), numeric(0))
  expect_true(is.na(pcens_cdf(obj, NA_real_, 2)))
  expect_identical(pcens_cdf(obj, c(3, 3, 1, 3), 2)[c(1, 2, 4)], rep(
    pcens_cdf(obj, 3, 2), 3
  ))
  # A positive delay has no mass at or below 0
  gamma_obj <- tlogis_object(families[[3]], 3, 0.5)
  expect_identical(pcens_cdf(gamma_obj, c(-3, 0), 2), c(0, 0))
  # pwindow = 0 is the delay CDF
  expect_identical(pcens_cdf(obj, 3, 0), pnorm(3, 3, 2))
})

test_that("a window that is not finite or positive uses the default method", {
  obj <- tlogis_object(normals[[1]], 0.5, 0.2)
  expect_null(.tlogis_plan(obj, Inf))
  expect_null(.tlogis_plan(obj, -1))
  expect_null(.tlogis_plan(obj, NA_real_))
})

test_that("bad location and scale are errors", {
  obj <- tlogis_object(normals[[1]], 0.5, 0.2)
  bad <- update(obj, primary_args = list(scale = -1))
  expect_error(pcens_cdf(bad, 1, 2), "scale")
  bad <- update(obj, primary_args = list(location = Inf))
  expect_error(pcens_cdf(bad, 1, 2), "location")
  bad <- update(obj, primary_args = list(scale = c(1, 2)))
  expect_error(pcens_cdf(bad, 1, 2), "scale")
})

test_that("pcens_pmf, pprimarycensored and truncation use the analytic CDF", {
  family <- normals[[1]]
  obj <- tlogis_object(family, 0.5, 0.2)
  pmf <- pcens_pmf(obj, 0:12, pwindow = 2, swindow = 1)
  cdf <- pcens_cdf(obj, 0:13, 2)
  expect_identical(pmf, diff(cdf))
  ref <- tlogis_reference(
    0:13, 2, 0.5, 0.2, exptilt_cdf(family), FALSE
  )
  expect_equal(pmf, diff(ref), tolerance = 1e-8)
  truncated <- pcens_pmf(obj, 1:12, pwindow = 2, L = 1, D = 13)
  expect_equal(
    truncated, diff(ref)[2:13] / (ref[14] - ref[2]),
    tolerance = 1e-8
  )
  expect_identical(
    pprimarycensored(
      c(1, 4), pnorm, dprimary = dtlogis, pwindow = 2,
      primary_args = list(location = 0.5, scale = 0.2),
      mean = 3, sd = 2
    ),
    pcens_cdf(obj, c(1, 4), 2)
  )
})

test_that("each endpoint is evaluated once for integer delays", {
  # The CDF at q needs the transforms at q and q - pwindow. Over neighbouring
  # integer q these are the same integer endpoints, and with an integer
  # location the split point q - location is one of them too. A transform
  # call evaluates each of its distinct points once.
  counts <- new.env()
  counts$points <- 0
  original <- .pcens_tilt_transform
  local_mocked_bindings(
    .pcens_tilt_transform = function(object, t, xi, upper = FALSE) {
      counts$points <- max(counts$points, length(unique(t)))
      original(object, t, xi, upper)
    }
  )
  family <- normals[[1]]
  for (location in c(-1, 0.5, 1, 2, 3.5)) {
    counts$points <- 0
    obj <- tlogis_object(family, location, 0.7)
    pcens_cdf(obj, 0:20, 2)
    # 21 delays and their neighbours at q - 2 give 23 endpoints. The points
    # q - location of a location inside the window that is not an integer
    # are all different.
    split_points <- if (location > 0 && location < 2 && location %% 1 != 0) {
      21
    } else {
      0
    }
    expect_lte(
      counts$points, 23 + split_points,
      label = paste("location", location)
    )
    expect_gte(counts$points, 23)
  }
})

test_that("non-parametric delays use the truncated logistic CDF", {
  boundaries <- c(0, 1, 3, 6, 10)
  pmf <- c(0.2, 0.3, 0.35, 0.15)
  obj <- new_pcens(
    pdist = pdiscretestep, dprimary = dtlogis,
    primary_args = list(location = 0.5, scale = 0.3),
    boundaries = boundaries, pmf = pmf
  )
  q <- c(0.5, 1.5, 3, 5, 8, 12)
  expect_equal(
    pcens_cdf(obj, q, 2),
    pcens_cdf(obj, q, 2, use_numeric = TRUE),
    tolerance = 1e-6
  )
})

# Narrow windows. The density of the primary has its mass in a few scales
# around the location (or around the edge nearest to it), so a quadrature
# that is not told where that is steps over it.
narrow_delays <- list(
  list(
    label = "exponential rate 1.5", pdist = pexp, args = list(rate = 1.5),
    positive = TRUE
  ),
  list(
    label = "gamma shape 2", pdist = pgamma,
    args = list(shape = 2, rate = 1), positive = TRUE
  ),
  list(
    label = "lognormal", pdist = plnorm,
    args = list(meanlog = 1, sdlog = 0.5), positive = TRUE
  ),
  list(
    label = "weibull", pdist = pweibull,
    args = list(shape = 2, scale = 2), positive = TRUE
  ),
  list(
    label = "normal", pdist = pnorm, args = list(mean = 3, sd = 2),
    positive = FALSE
  )
)

# Relative difference, absolute below 1e-10 where `stats::integrate()` has
# its absolute tolerance
narrow_diff <- function(actual, expected) {
  max(abs(actual - expected) / pmax(expected, 1e-10))
}

test_that("the numerical CDF resolves a narrow truncated logistic primary", {
  pwindow <- 2
  q <- c(0.5, 1, 3, 10)
  for (delay in narrow_delays) {
    cdf <- function(x) do.call(delay$pdist, c(list(x), delay$args))
    for (location in c(-0.5, 0.7, 2.5)) {
      for (scale in c(0.02, 0.005, 0.001)) {
        obj <- tlogis_object(delay, location, scale)
        expected <- tlogis_reference(
          q, pwindow, location, scale, cdf, delay$positive
        )
        label <- tlogis_label(delay, pwindow, location, scale)
        # The default dispatch, which falls back where the series does not
        # apply, and the forced numerical method
        expect_lt(
          narrow_diff(pcens_cdf(obj, q, pwindow), expected), 1e-6,
          label = label
        )
        expect_lt(
          narrow_diff(
            pcens_cdf(obj, q, pwindow, use_numeric = TRUE), expected
          ), 1e-6,
          label = paste("numeric:", label)
        )
      }
    }
  }
})

test_that("the numerical CDF resolves a named wrapper of dtlogis", {
  # A wrapper with the name attribute gets the tlogis classes, and so the
  # same dispatch as dtlogis, including the break points of a narrow primary
  wrapper <- function(x, min = 0, max = 1, location = 0, scale = 1,
                      log = FALSE) {
    dtlogis(x, min, max, location, scale, log)
  }
  named <- add_name_attribute(wrapper, "dtlogis")
  pwindow <- 2
  q <- c(0.72, 1, 1.5, 3)
  cdf <- function(x) pgamma(x, shape = 2, rate = 1)
  for (scale in c(1e-3, 1e-4)) {
    obj <- new_pcens(
      pdist = pgamma, dprimary = named,
      primary_args = list(location = 0.7, scale = scale),
      shape = 2, rate = 1
    )
    expect_s3_class(obj, "pcens_pgamma_dtlogis")
    expected <- tlogis_reference(q, pwindow, 0.7, scale, cdf)
    expect_lt(
      narrow_diff(pcens_cdf(obj, q, pwindow), expected), 1e-6,
      label = paste("scale", scale)
    )
    expect_lt(
      narrow_diff(pcens_cdf(obj, q, pwindow, use_numeric = TRUE), expected),
      1e-6, label = paste("numeric, scale", scale)
    )
  }
})

test_that("the normal tilt transform keeps the shift of a small sd", {
  # The CDF at q = mean - k sd is, with p = sd z,
  # int Phi(-k - z) f(sd z) sd dz, so the window is in units of sd
  normal_reference <- function(mean, sd, k, pwindow, location, scale) {
    integrand <- function(z) {
      stats::pnorm(-k - z) * dtlogis(sd * z, 0, pwindow, location, scale) * sd
    }
    upper <- min(40, pwindow / sd)
    stats::integrate(
      integrand, 0, upper, rel.tol = 1e-13, abs.tol = 0, subdivisions = 2000L
    )$value
  }
  cases <- list(
    list(mean = 1, sd = 1e-8, location = -0.5, scale = 1, pwindow = 1),
    list(mean = 30, sd = 1e-6, location = 0.5, scale = 20, pwindow = 1),
    list(mean = 30, sd = 1e-4, location = 0.5, scale = 20, pwindow = 1),
    list(mean = 100, sd = 1e-3, location = 0.5, scale = 20, pwindow = 1)
  )
  for (cs in cases) {
    for (k in c(0, 2)) {
      obj <- new_pcens(
        pdist = pnorm, dprimary = dtlogis,
        primary_args = list(location = cs$location, scale = cs$scale),
        mean = cs$mean, sd = cs$sd
      )
      expected <- normal_reference(
        cs$mean, cs$sd, k, cs$pwindow, cs$location, cs$scale
      )
      actual <- pcens_cdf(obj, cs$mean - k * cs$sd, cs$pwindow)
      expect_lt(
        abs(actual / expected - 1), 1e-6,
        label = paste(toString(unlist(cs)), "k =", k)
      )
    }
  }
})
