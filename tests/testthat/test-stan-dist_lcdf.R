skip_on_cran()

test_that(
  "Stan dist_lcdf matches R CDF for each distribution at R-assigned ID",
  {
    delay <- 2.5

    # Lognormal: R dist_id = 1
    stan_result <- dist_lcdf(delay, c(0.5, 0.8), 1L)
    r_result <- plnorm(delay, 0.5, 0.8, log.p = TRUE)
    expect_equal(
      stan_result, r_result,
      tolerance = 1e-6,
      info = "dist_id 1 should be lognormal"
    )

    # Gamma: R dist_id = 2
    stan_result <- dist_lcdf(delay, c(2.0, 1.0), 2L)
    r_result <- pgamma(delay, 2.0, 1.0, log.p = TRUE)
    expect_equal(
      stan_result, r_result,
      tolerance = 1e-6,
      info = "dist_id 2 should be gamma"
    )

    # Weibull: R dist_id = 3
    stan_result <- dist_lcdf(delay, c(1.5, 2.0), 3L)
    r_result <- pweibull(delay, 1.5, 2.0, log.p = TRUE)
    expect_equal(
      stan_result, r_result,
      tolerance = 1e-6,
      info = "dist_id 3 should be weibull"
    )

    # Exponential: R dist_id = 4
    stan_result <- dist_lcdf(delay, c(0.5, 0.0), 4L)
    r_result <- pexp(delay, 0.5, log.p = TRUE)
    expect_equal(
      stan_result, r_result,
      tolerance = 1e-6,
      info = "dist_id 4 should be exponential"
    )

    # Beta: R dist_id = 9
    stan_result <- dist_lcdf(0.7, c(2.0, 3.0), 9L)
    r_result <- pbeta(0.7, 2.0, 3.0, log.p = TRUE)
    expect_equal(
      stan_result, r_result,
      tolerance = 1e-6,
      info = "dist_id 9 should be beta"
    )

    # Cauchy: R dist_id = 12
    stan_result <- dist_lcdf(delay, c(0.0, 1.0), 12L)
    r_result <- pcauchy(delay, 0.0, 1.0, log.p = TRUE)
    expect_equal(
      stan_result, r_result,
      tolerance = 1e-6,
      info = "dist_id 12 should be cauchy"
    )

    # Chi-square: R dist_id = 13
    stan_result <- dist_lcdf(delay, c(3.0, 0.0), 13L)
    r_result <- pchisq(delay, 3.0, log.p = TRUE)
    expect_equal(
      stan_result, r_result,
      tolerance = 1e-6,
      info = "dist_id 13 should be chi-square"
    )

    # Inverse gamma: R dist_id = 16
    # Stan inv_gamma(alpha, beta): CDF = upper_gamma(alpha, beta/x)
    # In R: pgamma(beta/x, alpha, 1, lower.tail = FALSE)
    ig_alpha <- 3.0
    ig_beta <- 2.0
    stan_result <- dist_lcdf(delay, c(ig_alpha, ig_beta), 16L)
    r_result <- pgamma(
      ig_beta / delay, ig_alpha, 1,
      lower.tail = FALSE, log.p = TRUE
    )
    expect_equal(
      stan_result, r_result,
      tolerance = 1e-6,
      info = "dist_id 16 should be inverse gamma"
    )

    # Logistic: R dist_id = 17
    stan_result <- dist_lcdf(delay, c(1.0, 0.5), 17L)
    r_result <- plogis(delay, 1.0, 0.5, log.p = TRUE)
    expect_equal(
      stan_result, r_result,
      tolerance = 1e-6,
      info = "dist_id 17 should be logistic"
    )
  }
)

test_that(
  "Stan dist_lcdf IDs are consistent with pcd_stan_dist_id()",
  {
    delay <- 3.0

    dists <- list(
      list(
        name = "lnorm",
        params = c(0.5, 0.8),
        r_cdf = function(d, p) plnorm(d, p[1], p[2], log.p = TRUE)
      ),
      list(
        name = "gamma",
        params = c(2.0, 1.0),
        r_cdf = function(d, p) pgamma(d, p[1], p[2], log.p = TRUE)
      ),
      list(
        name = "weibull",
        params = c(1.5, 2.0),
        r_cdf = function(d, p) pweibull(d, p[1], p[2], log.p = TRUE)
      ),
      list(
        name = "exp",
        params = c(0.5, 0.0),
        r_cdf = function(d, p) pexp(d, p[1], log.p = TRUE)
      ),
      list(
        name = "beta",
        delay = 0.7,
        params = c(2.0, 3.0),
        r_cdf = function(d, p) pbeta(d, p[1], p[2], log.p = TRUE)
      ),
      list(
        name = "cauchy",
        params = c(0.0, 1.0),
        r_cdf = function(d, p) pcauchy(d, p[1], p[2], log.p = TRUE)
      ),
      list(
        name = "chisq",
        params = c(3.0, 0.0),
        r_cdf = function(d, p) pchisq(d, p[1], log.p = TRUE)
      ),
      list(
        name = "logis",
        params = c(1.0, 0.5),
        r_cdf = function(d, p) plogis(d, p[1], p[2], log.p = TRUE)
      )
    )

    for (dist in dists) {
      d <- dist$delay %||% delay
      dist_id <- pcd_stan_dist_id(dist$name)
      stan_result <- dist_lcdf(d, dist$params, dist_id)
      r_result <- dist$r_cdf(d, dist$params)
      expect_equal(
        stan_result, r_result,
        tolerance = 1e-6,
        info = sprintf(
          "dist_lcdf with pcd_stan_dist_id('%s') = %d should match R",
          dist$name, dist_id
        )
      )
    }
  }
)
