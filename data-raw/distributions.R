pcd_distributions <- data.frame(
  name = c(
    "lnorm", "gamma", "weibull", "exp", "gengamma", "nbinom",
    "pois", "bern", "beta", "binom", "cat", "cauchy", "chisq",
    "dirich", "gumbel", "invgamma", "logis",
    "norm", "invchisq", "dblexp", "pareto",
    "scaleinvchisq", "student_t", "unif", "vonmises"
  ),
  pdist = c(
    "plnorm", "pgamma", "pweibull", "pexp", NA, "pnbinom",
    "ppois", NA, "pbeta", "pbinom", NA, "pcauchy", "pchisq",
    NA, "pgumbel", NA, "plogis",
    "pnorm", NA, NA, NA,
    NA, "pt", "punif", NA
  ),
  aliases = c(
    "lognormal", "gamma", "weibull", "exponential",
    "generalized gamma", "negative binomial", "poisson",
    "bernoulli", "beta", "binomial", "categorical", "cauchy",
    "chi-square", "dirichlet", "gumbel", "inverse gamma",
    "logistic",
    "normal", "inverse chi-square", "double exponential",
    "pareto", "scaled inverse chi-square", "student t",
    "uniform", "von mises"
  ),
  stan_id = 1:25,
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
