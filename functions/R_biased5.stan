data {
  int<lower=1> N;                      // Number of observed days
  int<lower=1> h;                      // Number of forecast days
  int<lower=1> L;                      // Length of generation time distribution
  int<lower=1> k;                      // Number of strata
  int<lower=0> incidence_strat[N, k];  // Observed incidence per group
  vector<lower=0>[L] gen_weights;      // Generation time distribution (normalized)
  int<lower=1> S0[k];                  // Initial susceptibles
}

transformed data {
  int T = N + h;
  int B = (T + 2) / 3;                 // ceil(T / 3)
}

parameters {
  vector<lower=0>[k] R0;
  real<lower=0> rw_sigma;
  real<lower=0> phi;
  vector<lower=0>[k] seed_infections;
  real<lower=0> drift;
  
  matrix[B, k] log_Rt_block;           // RW state every 3 days (log scale, unconstrained)
}

transformed parameters {
  matrix<lower=0>[B, k] Rt_block;      // exp(log_Rt_block)
  matrix<lower=0>[T, k] infections;
  matrix<lower=0>[T, k] cumulative_infections;
  
  Rt_block = exp(log_Rt_block);
  
  // Build infections forward in time, stratum by stratum.
  // NOTE: we do NOT build Rt[t,j]; we index Rt_block directly via b = 1 + (t-1)/3.
  for (j in 1:k) {
    infections[1, j] = seed_infections[j];
    cumulative_infections[1, j] = seed_infections[j];
    
    for (t in 2:T) {
      int b = 1 + (t - 1) / 3;
      real Rt_val = Rt_block[b, j];
      
      int S = min(L, t - 1);           // number of lags available
      real infectivity = 0;
      
      // infectivity = sum_{s=1..S} infections[t-s,j] * gen_weights[s]
      // (branchless; only sum feasible lags)
      for (s in 1:S) {
        infectivity += infections[t - s, j] * gen_weights[s];
      }
      
      // susceptible fraction
      real sus_frac = 1 - cumulative_infections[t - 1, j] / S0[j];
      
      // Avoid non-differentiable fmax; keep strictly positive mean
      infections[t, j] = Rt_val * infectivity * sus_frac + 1e-6;
      cumulative_infections[t, j] = cumulative_infections[t - 1, j] + infections[t, j];
    }
  }
}

model {
  // Priors
  R0 ~ lognormal(log(10), 0.4);
  rw_sigma ~ lognormal(log(0.1), 0.2);
  phi ~ lognormal(log(30), 0.5);
  drift ~ lognormal(-2.5, 0.3);
  seed_infections ~ lognormal(log(1), 0.5);
  
  // Random walk prior on log_Rt_block
  for (j in 1:k) {
    log_Rt_block[1, j] ~ normal(log(R0[j]), 0.1);
    for (b in 2:B) {
      log_Rt_block[b, j] ~ normal(log_Rt_block[b - 1, j] - drift, rw_sigma);
    }
  }
  
  // Likelihood (vectorised over time for each stratum)
  for (j in 1:k) {
    incidence_strat[1:N, j] ~ neg_binomial_2(infections[1:N, j], phi);
  }
}

generated quantities {
  matrix[h, k] forecast;
  
  for (j in 1:k) {
    for (t in 1:h) {
      forecast[t, j] = neg_binomial_2_rng(infections[N + t, j], phi);
    }
  }
}
