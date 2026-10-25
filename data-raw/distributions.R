pcd_distributions <- data.frame(
  name = c(
    "lnorm", "gamma", "weibull", "exp", "gengamma", "nbinom",
    "pois", "bern", "beta", "binom", "cat", "cauchy", "chisq",
    "dirich", "gumbel", "invgamma", "logis",
    "norm", "invchisq", "dblexp", "pareto",
    "scaleinvchisq", "student_t", "unif", "vonmises",
    "discretestep", "discretehazard_rw", "discretehazard_re",
    "loglogistic"
  ),
  pdist = c(
    "plnorm", "pgamma", "pweibull", "pexp", NA, "pnbinom",
    "ppois", NA, "pbeta", "pbinom", NA, "pcauchy", "pchisq",
    NA, "pgumbel", NA, "plogis",
    "pnorm", NA, NA, NA,
    NA, "pt", "punif", NA,
    "pdiscretestep", "pdiscretehazard", "pdiscretehazard",
    "pllogis"
  ),
  aliases = c(
    "lognormal", "gamma", "weibull", "exponential",
    "generalized gamma", "negative binomial", "poisson",
    "bernoulli", "beta", "binomial", "categorical", "cauchy",
    "chi-square", "dirichlet", "gumbel", "inverse gamma",
    "logistic",
    "normal", "inverse chi-square", "double exponential",
    "pareto", "scaled inverse chi-square", "student t",
    "uniform", "von mises",
    "nonparametric", "hazard random walk", "hazard random effect",
    "log-logistic"
  ),
  stan_id = c(1:25, 26L, 27L, 28L, 29L),
  stringsAsFactors = FALSE
)

pcd_primary_distributions <- data.frame(
  name = c("unif", "expgrowth", "tlogis", "tgumbel"),
  dprimary = c("dunif", "dexpgrowth", "dtlogis", "dtgumbel"),
  pprimary = c("punif", "pexpgrowth", "ptlogis", "ptgumbel"),
  aliases = c(
    "uniform", "exponential growth", "truncated logistic", "truncated gumbel"
  ),
  stan_id = 1:4,
  stringsAsFactors = FALSE
)

usethis::use_data( # nolint: namespace_linter
  pcd_distributions,
  pcd_primary_distributions,
  overwrite = TRUE
)
