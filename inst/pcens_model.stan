functions {
  #include primary_censored_dist.stan
  #include primary_censored_dist_cdf_analytical.stan
  #include expgrowth.stan

  real partial_sum(array[] int dummy, int start, int end,
                   array[] int d, array[] int d_upper, array[] int n,
                   array[] int pwindow, array[] int D,
                   int dist_id, array[] real params,
                   int primary_dist_id, array[] real primary_params) {
    real partial_target = 0;
    for (i in start:end) {
      partial_target += n[i] * primary_censored_dist_lpmf(
        d[i] | dist_id, params,
        pwindow[i], d_upper[i], D[i],
        primary_dist_id, primary_params
      );
    }
    return partial_target;
  }
}

data {
  int<lower=0> N;  // number of observations
  array[N] int<lower=0> d;     // observed delays
  array[N] int<lower=0> d_upper;     // observed delays upper bound
  array[N] int<lower=0> n;     // number of occurrences for each delay
  array[N] int<lower=0> pwindow; // primary censoring window
  array[N] int<lower=0> D; // maximum delay
  int<lower=1, upper=17> dist_id; // distribution identifier
  int<lower=1, upper=2> primary_dist_id; // primary distribution identifier
  int<lower=0> n_params; // number of distribution parameters
  int<lower=0> n_primary_params; // number of primary distribution parameters
  int<lower=0, upper=1> compute_log_lik; // whether to compute log likelihood
  int<lower=0, upper=1> use_reduce_sum; // whether to use reduce_sum
  vector[n_params] param_lower_bounds; // lower bounds for parameters
  vector[n_params] param_upper_bounds; // upper bounds for parameters
  vector[n_primary_params] primary_param_lower_bounds; // lower bounds for primary parameters
  vector[n_primary_params] primary_param_upper_bounds; // upper bounds for primary parameters
  array[n_params] real prior_location; // location parameters for priors
  array[n_params] real prior_scale; // scale parameters for priors
  array[n_primary_params] real primary_prior_location; // location parameters for primary priors
  array[n_primary_params] real primary_prior_scale; // scale parameters for primary priors
}

transformed data {
  array[N] int indexes = linspaced_int_array(N, 1, N);
}

parameters {
  vector<lower=param_lower_bounds, upper=param_upper_bounds>[n_params] params; // distribution parameters
  vector<lower=primary_param_lower_bounds, upper=primary_param_upper_bounds>[n_primary_params] primary_params; // primary distribution parameters
}


model {
  // Priors
  for (i in 1:n_params) {
    params[i] ~ normal(prior_location[i], prior_scale[i]);
  }
  if (n_primary_params) {
    for (i in 1:n_primary_params) {
      primary_params[i] ~ normal(
        primary_prior_location[i],
        primary_prior_scale[i]
      );
    }
  }

  // Likelihood
  if (use_reduce_sum) {
    target += reduce_sum(partial_sum, indexes, 1, d,
                         d_upper, n, pwindow, D, dist_id, to_array_1d(params),
                         primary_dist_id, to_array_1d(primary_params));
  } else {
    for (i in 1:N) {
      target += n[i] * primary_censored_dist_lpmf(
        d[i] | dist_id, to_array_1d(params),
        pwindow[i], d_upper[i], D[i],
        primary_dist_id, to_array_1d(primary_params)
      );
    }
  }
}

generated quantities {
  vector[compute_log_lik ? N : 0] log_lik;
  if (compute_log_lik) {
    for (i in 1:N) {
      log_lik[i] = primary_censored_dist_lpmf(
        d[i] | dist_id, to_array_1d(params),
        pwindow[i], d_upper[i], D[i],
        primary_dist_id, to_array_1d(primary_params)
      );
    }
  }
}
