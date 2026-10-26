functions {
  #include nonparametric.stan
  #include primarycensored.stan
  #include primarycensored_ode.stan
  #include primarycensored_analytical_cdf.stan
  #include expgrowth.stan

  real partial_sum(array[] int rows, int start, int end,
                   array[] int d, array[] int d_upper, array[] int n,
                   array[] int pwindow, data array[] real L,
                   data array[] real D,
                   int dist_id, array[] real params,
                   int primary_id, array[] real primary_params) {
    real partial_target = 0;
    for (i in rows) {
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
  // Distribution identifier. 1..25 = parametric families. The
  // non-parametric step CDF has the same likelihood for all three of
  // the entries below; they only differ in the prior on the PMF /
  // hazards:
  //   26 = step CDF, Dirichlet prior on the PMF;
  //   27 = step CDF, Gaussian random walk on the logit hazards;
  //   28 = step CDF, IID logit random effects on the hazards.
  int<lower=1, upper=28> dist_id;
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
  // Non-parametric block. `dist_id` (26/27/28) selects the family;
  // priors flow through `np_*` data fields populated from R by
  // `pcd_as_stan_data()` from the user-facing `priors` argument.
  //   dist_id == 26: Dirichlet on the simplex (concentration in
  //                  `np_dirichlet_alpha`).
  //   dist_id == 27: Gaussian random walk on the logit hazards
  //                  (`alpha ~ N(mean, sd)`, `log_sigma ~ N(mean, sd)`).
  //   dist_id == 28: IID logit random effects on the hazards (same
  //                  priors as 27).
  int<lower=0, upper=1> nonparametric;
  int<lower=0> K_np; // number of bins
  vector[K_np + 1] np_boundaries; // bin boundaries (length K_np + 1)
  vector<lower=0>[K_np] np_dirichlet_alpha; // Dirichlet concentration
  real np_alpha_mean; // hazard intercept prior mean (logit scale)
  real<lower=0> np_alpha_sd; // hazard intercept prior sd
  real np_log_sigma_mean; // log_sigma prior mean
  real<lower=0> np_log_sigma_sd; // log_sigma prior sd
}

transformed data {
  // Sizes for the non-parametric parameter branches. Exactly one is K_np;
  // both are 0 when nonparametric = 0. Stan forbids `simplex[0]`, so we
  // pad to size 1 when the simplex branch is inactive and ignore it.
  int simplex_active = (nonparametric == 1 && dist_id == 26) ? 1 : 0;
  int K_simplex = simplex_active == 1 ? K_np : 1;
  // Both hazard parameterisations share the same set of unconstrained
  // parameters (alpha, sigma, eps_1..eps_{K-1}); they only differ in
  // how `np_weights` is assembled in `transformed parameters`.
  int K_hazard = (nonparametric == 1 && (dist_id == 27
                  || dist_id == 28)) ? K_np : 0;
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

  // Rows with swindow = 1 and d >= 0 that share pwindow, L and D take their
  // log PMFs from one primarycensored_sone_lpmf_vectorized() call, which
  // computes each CDF (and F(D)) once rather than once per row. Groups are
  // found within runs of equal D after sorting on D. A group is only used
  // when its rows cover at least half of the delays 0:max_delay, so sparse
  // groups do not compute many unused CDFs; their rows, and all other rows,
  // stay on the per-row primarycensored_lpmf() path.
  array[N] int row_group = rep_array(0, N);
  array[N] int group_first; // first row of each group
  array[N] int group_max_d;
  array[N] int group_size;
  int n_groups = 0;
  {
    array[N] int by_D = sort_indices_asc(D);
    int run_start = 1;
    for (k in 1:N) {
      int i = by_D[k];
      if (k > 1 && D[i] != D[by_D[k - 1]]) {
        run_start = n_groups + 1;
      }
      if (d_upper[i] - d[i] == 1 && d[i] >= 0) {
        int g = 0;
        for (h in run_start:n_groups) {
          if (pwindow[group_first[h]] == pwindow[i]
              && L[group_first[h]] == L[i]) {
            g = h;
            break;
          }
        }
        if (g == 0) {
          n_groups += 1;
          g = n_groups;
          group_first[g] = i;
          group_max_d[g] = d[i];
          group_size[g] = 0;
        }
        row_group[i] = g;
        group_max_d[g] = max(group_max_d[g], d[i]);
        group_size[g] += 1;
      }
    }
  }
  array[n_groups] int use_group;
  for (g in 1:n_groups) {
    use_group[g] = group_max_d[g] + 1 <= 2 * group_size[g];
  }
  int n_vec_groups = sum(use_group);
  int n_scalar = 0;
  for (i in 1:N) {
    n_scalar += row_group[i] == 0 || !use_group[row_group[i]];
  }
  // Per vectorised group: its first row (for pwindow, L and D), max delay
  // and the span of its rows in vec_rows.
  array[n_vec_groups] int vec_first;
  array[n_vec_groups] int vec_max_d;
  array[n_vec_groups] int vec_start;
  array[n_vec_groups] int vec_end;
  array[N - n_scalar] int vec_rows;
  array[n_scalar] int scalar_rows;
  {
    array[n_groups] int vec_id = rep_array(0, n_groups);
    int pos = 1;
    int v = 0;
    for (g in 1:n_groups) {
      if (use_group[g]) {
        v += 1;
        vec_id[g] = v;
        vec_first[v] = group_first[g];
        vec_max_d[v] = group_max_d[g];
        vec_start[v] = pos;
        pos += group_size[g];
        vec_end[v] = vec_start[v] - 1;
      }
    }
    int s = 0;
    for (i in 1:N) {
      if (row_group[i] > 0 && use_group[row_group[i]]) {
        v = vec_id[row_group[i]];
        vec_end[v] += 1;
        vec_rows[vec_end[v]] = i;
      } else {
        s += 1;
        scalar_rows[s] = i;
      }
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
  // distributions this is just the parametric `params`; for dist_id
  // 26/27/28 it is [boundaries (K_np + 1), pmf or hazards (K_np)].
  array[n_lpmf_params] real lpmf_params;
  vector[K_np] np_weights; // pmf for dist_id 26, hazards for 27/28
  if (nonparametric == 1) {
    if (dist_id == 26) {
      np_weights = np_pmf;
    } else {
      // Build hazards on the logit scale, pinning the final hazard to 1
      // so the implied PMF sums to 1. dist_id == 27 uses a Gaussian
      // random walk (logit(h_i) = logit(h_{i-1}) + sigma * eps_i);
      // dist_id == 28 uses IID logit random effects around the
      // intercept (logit(h_i) = alpha + sigma * eps_i, eps_i ~ N(0, 1)).
      vector[K_np] logit_h;
      logit_h[1] = np_alpha[1];
      if (K_np > 1) {
        if (dist_id == 27) {
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
    if (dist_id == 26) {
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
  for (v in 1:n_vec_groups) {
    int first = vec_first[v];
    vector[vec_max_d[v] + 1] group_lpmf = primarycensored_sone_lpmf_vectorized(
      vec_max_d[v], L[first], D[first], dist_id, lpmf_params, pwindow[first],
      primary_id, to_array_1d(primary_params)
    );
    for (j in vec_start[v]:vec_end[v]) {
      int i = vec_rows[j];
      target += n[i] * group_lpmf[d[i] + 1];
    }
  }
  if (use_reduce_sum) {
    target += reduce_sum(partial_sum, scalar_rows, 1, d,
                         d_upper, n, pwindow, L, D, dist_id, lpmf_params,
                         primary_id, to_array_1d(primary_params));
  } else {
    for (i in scalar_rows) {
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
