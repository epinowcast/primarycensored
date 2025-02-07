pcd_distributions <- data.frame(
  name = c(
    "lnorm", "gamma", "weibull", "exp", "gengamma", "nbinom",
    "pois", "bern", "beta", "binom", "cat", "cauchy", "chisq",
    "dirich", "gumbel", "invgamma", "logis"
  ),
  pdist = c(
    "plnorm", "pgamma", "pweibull", "pexp", NA, "pnbinom",
    "ppois", NA, "pbeta", "pbinom", NA, "pcauchy", "pchisq",
    NA, "pgumbel", NA, "plogis"
  ),
  aliases = c(
    "lognormal", "gamma", "weibull", "exponential", "generalized gamma",
    "negative binomial", "poisson", "bernoulli", "beta", "binomial",
    "categorical", "cauchy", "chi-square", "dirichlet", "gumbel",
    "inverse gamma", "logistic"
  ),
  stan_id = 1:17,
  stringsAsFactors = FALSE
)

pcd_primary_distributions <- data.frame(
  name = c("unif", "expgrowth"),
  dprimary = c("dunif", "dexpgrowth"),
  aliases = c("uniform", "exponential growth"),
  stan_id = 1:2,
  stringsAsFactors = FALSE
)

usethis::use_data(
  pcd_distributions,
  pcd_primary_distributions,
  overwrite = TRUE
)
