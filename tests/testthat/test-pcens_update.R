update_cases <- list(
  gamma_unif = list(
    pdist = pgamma, dprimary = dunif, primary_args = list(),
    old = list(shape = 1.5, scale = 2), new = list(shape = 3, scale = 0.7)
  ),
  lnorm_unif = list(
    pdist = plnorm, dprimary = dunif, primary_args = list(),
    old = list(meanlog = 0, sdlog = 1), new = list(meanlog = 1.2, sdlog = 0.4)
  ),
  weibull_unif = list(
    pdist = pweibull, dprimary = dunif, primary_args = list(),
    old = list(shape = 1.5, scale = 2), new = list(shape = 2.5, scale = 4)
  ),
  gamma_expgrowth = list(
    pdist = pgamma, dprimary = dexpgrowth, primary_args = list(r = 0.2),
    old = list(shape = 1.5, scale = 2), new = list(shape = 3, scale = 0.7)
  )
)

test_that("update.pcens matches a freshly built object", {
  for (case in update_cases) {
    obj <- do.call(
      new_pcens,
      c(list(case$pdist, case$dprimary, case$primary_args), case$old)
    )
    updated <- do.call(update, c(list(obj), case$new))
    fresh <- do.call(
      new_pcens,
      c(list(case$pdist, case$dprimary, case$primary_args), case$new)
    )
    expect_identical(updated, fresh)
    expect_s3_class(updated, class(obj), exact = TRUE)
  }
})

test_that("update.pcens + pcens_cdf matches pprimarycensored", {
  q <- seq(0, 10, by = 0.25)
  for (case in update_cases) {
    obj <- do.call(
      new_pcens,
      c(list(case$pdist, case$dprimary, case$primary_args), case$old)
    )
    updated <- do.call(update, c(list(obj), case$new))
    for (pwindow in c(1, 2)) {
      expected <- do.call(
        pprimarycensored,
        c(
          list(
            q, case$pdist,
            pwindow = pwindow, dprimary = case$dprimary,
            primary_args = case$primary_args
          ),
          case$new
        )
      )
      expect_identical(pcens_cdf(updated, q, pwindow), expected)
    }
  }
})

test_that("update.pcens + pcens_pmf matches dprimarycensored", {
  settings <- list(
    list(x = 0:20, pwindow = 1, swindow = 1, L = -Inf, D = Inf),
    list(x = 0:9, pwindow = 1, swindow = 1, L = -Inf, D = 10),
    list(x = 1:9, pwindow = 1, swindow = 1, L = 1, D = 10),
    list(x = seq(0, 8, by = 2), pwindow = 2, swindow = 2, L = 0, D = 10),
    list(x = c(0, 3, 9), pwindow = 1, swindow = 2, L = -Inf, D = 10)
  )
  for (case in update_cases) {
    obj <- do.call(
      new_pcens,
      c(list(case$pdist, case$dprimary, case$primary_args), case$old)
    )
    updated <- do.call(update, c(list(obj), case$new))
    for (s in settings) {
      for (log in c(FALSE, TRUE)) {
        expected <- suppressMessages(do.call(
          dprimarycensored,
          c(
            list(
              s$x, case$pdist,
              pwindow = s$pwindow, swindow = s$swindow,
              L = s$L, D = s$D, dprimary = case$dprimary,
              primary_args = case$primary_args, log = log
            ),
            case$new
          )
        ))
        actual <- suppressMessages(pcens_pmf(
          updated, s$x,
          pwindow = s$pwindow, swindow = s$swindow,
          L = s$L, D = s$D, log = log
        ))
        expect_identical(actual, expected)
      }
    }
  }
})

test_that("update.pcens merges delay parameters", {
  obj <- new_pcens(pgamma, dunif, list(), shape = 1.5, scale = 2)
  updated <- update(obj, scale = 3)
  expect_identical(updated$args, list(shape = 1.5, scale = 3))
  expect_identical(update(obj), obj)
})

test_that("update.pcens adds parameters in the formals of pdist", {
  obj <- new_pcens(pgamma, dunif, list(), shape = 1.5)
  updated <- update(obj, rate = 0.5)
  expect_identical(updated$args, list(shape = 1.5, rate = 0.5))
})

test_that("update.pcens accepts any name when pdist takes dots", {
  pdist <- add_name_attribute(function(q, ...) pgamma(q, ...), "pgamma")
  obj <- new_pcens(pdist, dunif, list(), shape = 1.5)
  updated <- update(obj, rate = 2)
  expect_identical(updated$args, list(shape = 1.5, rate = 2))
})

test_that("update.pcens merges primary_args and keeps the alias", {
  obj <- new_pcens(
    pgamma, dexpgrowth, list(r = 0.2), shape = 1.5, scale = 2
  )
  updated <- update(obj, primary_args = list(r = 0.5))
  fresh <- new_pcens(
    pgamma, dexpgrowth, list(r = 0.5), shape = 1.5, scale = 2
  )
  expect_identical(updated, fresh)
  expect_identical(updated$dprimary_args, list(r = 0.5))
  expect_identical(updated$pprimary, obj$pprimary)
  expect_s3_class(updated, "pcens_pgamma_dexpgrowth")
})

test_that("update.pcens errors on invalid parameters", {
  obj <- new_pcens(pgamma, dunif, list(), shape = 1.5, scale = 2)
  expect_error(update(obj, 2), "must be named")
  expect_error(update(obj, shap = 2), "shap")
  expect_error(
    update(obj, primary_args = 0.5),
    "primary_args must be a list"
  )
  expect_error(
    update(obj, primary_args = list(0.5)),
    "primary_args must be named"
  )
})
