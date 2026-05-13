functions {
  #include nonparametric.stan
  #include primarycensored.stan
  #include primarycensored_ode.stan
  #include primarycensored_analytical_cdf.stan
  #include expgrowth.stan

  real partial_sum(array[] int dummy, int start, int end,
                   array[] int d, array[] int d_upper, array[] int n,
                   array[] int pwindow, data array[] real L,
                   data array[] real D,
                   int dist_id, array[] real params,
                   int primary_id, array[] real primary_params) {
    real partial_target = 0;
    for (i in start:end) {
      partial_target += n[i] * primarycensored_lpmf(
        d[i] | dist_id, params,
        pwindow[i], d_upper[i], L[i], D[i],
        primary_id, primary_params
      );
    }
    return partial_target;
  }
}

data {
  int<lower=0> N;  // number of observations
  // observed delays (may be negative for signed-support delays)
  array[N] int d;
  array[N] int d_upper;     // observed delays upper bound
  array[N] int<lower=0> n;     // number of occurrences for each delay
  array[N] int<lower=0> pwindow; // primary censoring window
  // lower truncation; -Inf for no lower truncation
  array[N] real L;
  // upper truncation; +Inf for no upper truncation
  array[N] real D;
  int<lower=1, upper=27> dist_id; // distribution identifier
  int<lower=1, upper=2> primary_id; // primary distribution identifier
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
  // Non-parametric block (dist_id 26 = discretestep, 27 = discretehazard).
  // Sized 0 when nonparametric = 0 so the parametric path is unaffected.
  int<lower=0, upper=1> nonparametric;
  int<lower=0> K_np; // number of bins
  vector[K_np + 1] np_boundaries; // bin boundaries (length K_np + 1)
  int<lower=1, upper=3> np_paramtype; // 1 = simplex/Dirichlet, 2 = RW hazards, 3 = RE hazards
  vector<lower=0>[K_np] np_dirichlet_alpha; // Dirichlet concentration
  real np_alpha_mean; // hazard intercept prior mean (logit scale)
  real<lower=0> np_alpha_sd; // hazard intercept prior sd
  real np_log_sigma_mean; // log RW sd prior mean
  real<lower=0> np_log_sigma_sd; // log RW sd prior sd
}

transformed data {
  array[N] int indexes = linspaced_int_array(N, 1, N);
  // Sizes for the non-parametric parameter branches. Exactly one is K_np;
  // both are 0 when nonparametric = 0. Stan forbids `simplex[0]`, so we
  // pad to size 1 when the simplex branch is inactive and ignore it.
  int simplex_active = (nonparametric == 1 && np_paramtype == 1) ? 1 : 0;
  int K_simplex = simplex_active == 1 ? K_np : 1;
  // Both hazard parameterisations share the same set of unconstrained
  // parameters (alpha, sigma, eps_1..eps_{K-1}); the only difference is
  // how `np_weights` is assembled in `transformed parameters`.
  int K_hazard = (nonparametric == 1 && (np_paramtype == 2
                  || np_paramtype == 3)) ? K_np : 0;
  int n_eps = K_hazard > 0 ? K_hazard - 1 : 0;
  // Length of the params array consumed by primarycensored_lpmf:
  // - parametric path: n_params
  // - non-parametric: 2 * K_np + 1 = (K_np + 1) boundaries + K_np weights
  int n_lpmf_params = nonparametric == 1 ? 2 * K_np + 1 : n_params;
  for (i in 1:N) {
    if (d[i] < L[i]) {
      reject("d[", i, "] = ", d[i], " is below the lower truncation L[", i,
             "] = ", L[i]);
    }
    if (d_upper[i] > D[i]) {
      reject("d_upper[", i, "] = ", d_upper[i],
             " is above the upper truncation D[", i, "] = ", D[i]);
    }
  }
}

parameters {
  vector<lower=param_lower_bounds, upper=param_upper_bounds>[n_params] params; // distribution parameters
  vector<lower=primary_param_lower_bounds, upper=primary_param_upper_bounds>[n_primary_params] primary_params; // primary distribution parameters
  // Non-parametric parameters. At most one branch is active (size > 0).
  simplex[K_simplex] np_pmf;
  array[K_hazard > 0 ? 1 : 0] real np_alpha;
  array[K_hazard > 0 ? 1 : 0] real<lower=0> np_sigma;
  vector[n_eps] np_eps;
}

transformed parameters {
  // Assemble the array consumed by primarycensored_lpmf. For parametric
  // distributions this is just the parametric `params`; for dist_id 26/27
  // it is [boundaries (K_np + 1), pmf or hazards (K_np)].
  array[n_lpmf_params] real lpmf_params;
  vector[K_np] np_weights; // pmf for paramtype 1, hazards for paramtype 2
  if (nonparametric == 1) {
    if (np_paramtype == 1) {
      np_weights = np_pmf;
    } else {
      // Build hazards on the logit scale, pinning the final hazard to 1
      // so the implied PMF sums to 1. paramtype 2 uses a Gaussian random
      // walk (logit(h_i) = logit(h_{i-1}) + sigma * eps_i); paramtype 3
      // uses IID logit random effects around the intercept
      // (logit(h_i) = alpha + sigma * eps_i, eps_i ~ N(0, 1)).
      vector[K_np] logit_h;
      logit_h[1] = np_alpha[1];
      if (K_np > 1) {
        if (np_paramtype == 2) {
          for (i in 2:K_np) {
            logit_h[i] = logit_h[i - 1] + np_sigma[1] * np_eps[i - 1];
          }
        } else {
          for (i in 2:K_np) {
            logit_h[i] = np_alpha[1] + np_sigma[1] * np_eps[i - 1];
          }
        }
      }
      np_weights = inv_logit(logit_h);
      np_weights[K_np] = 1;
    }
    for (k in 1:(K_np + 1)) lpmf_params[k] = np_boundaries[k];
    for (k in 1:K_np) lpmf_params[K_np + 1 + k] = np_weights[k];
  } else {
    lpmf_params = to_array_1d(params);
  }
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
  if (nonparametric == 1) {
    if (np_paramtype == 1) {
      np_pmf ~ dirichlet(np_dirichlet_alpha);
    } else {
      // Shared priors for both hazard parameterisations (RW and RE).
      np_alpha[1] ~ normal(np_alpha_mean, np_alpha_sd);
      log(np_sigma[1]) ~ normal(np_log_sigma_mean, np_log_sigma_sd);
      target += -log(np_sigma[1]); // Jacobian for log transform
      np_eps ~ std_normal();
    }
  }

  // Likelihood
  if (use_reduce_sum) {
    target += reduce_sum(partial_sum, indexes, 1, d,
                         d_upper, n, pwindow, L, D, dist_id, lpmf_params,
                         primary_id, to_array_1d(primary_params));
  } else {
    for (i in 1:N) {
      target += n[i] * primarycensored_lpmf(
        d[i] | dist_id, lpmf_params,
        pwindow[i], d_upper[i], L[i], D[i],
        primary_id, to_array_1d(primary_params)
      );
    }
  }
}

generated quantities {
  vector[compute_log_lik ? N : 0] log_lik;
  if (compute_log_lik) {
    for (i in 1:N) {
      log_lik[i] = primarycensored_lpmf(
        d[i] | dist_id, lpmf_params,
        pwindow[i], d_upper[i], L[i], D[i],
        primary_id, to_array_1d(primary_params)
      );
    }
  }
}
