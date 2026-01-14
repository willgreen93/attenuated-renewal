data {
  int<lower=1> N;                      // Number of observed days
  int<lower=1> h;                      // Number of forecast days
  int<lower=1> L;                      // Length of generation time distribution
  int<lower=1> k;                      // Number of strata
  int<lower=0> incidence_strat[N, k]; // Observed incidence per group
  vector<lower=0>[L] gen_weights;     // Generation time distribution
  int<lower=1> S0[k];                 // Initial susceptible population
}

parameters {
  vector<lower=0>[k] R0;              // Initial reproduction number per group
  vector<lower=0>[k] delta;           // Decay rate per group
  real<lower=1, upper=N> t_change;    // Time when decay begins (shared)
  real<lower=0> phi;                  // Overdispersion
  vector<lower=0>[k] seed_infections; // Initial infections
}

transformed parameters {
  matrix<lower=0>[N + h, k] infections;
  matrix<lower=0>[N + h, k] cumulative_infections;
  
  for (j in 1:k) {
    infections[1, j] = seed_infections[j];
    cumulative_infections[1, j] = seed_infections[j];
    
    for (t in 2:(N + h)) {
      real Rt;
      if (t <= t_change) {
        Rt = R0[j];
      } else {
        Rt = R0[j] * exp(-delta[j] * (t - t_change));
      }
      
      // Renewal equation
      real infectivity = 0;
      for (s in 1:L) {
        if (t > s) {
          infectivity += infections[t - s, j] * gen_weights[s];
        }
      }
      
      infections[t, j] = fmax(1e-6, Rt * infectivity * (1-cumulative_infections[t-1, j]/S0[j]));
      cumulative_infections[t, j] = cumulative_infections[t - 1, j] + infections[t, j];
    }
  }
}

model {
  R0 ~ lognormal(log(10), 0.4);
  //delta ~ normal(0.02, 0.005);
  delta ~ lognormal(-3.94, 0.246);
  t_change ~ uniform(1, N); // flat prior over observed period
  phi ~ lognormal(log(30), 0.5);
  seed_infections ~ lognormal(log(1), 0.5);

  for (j in 1:k) {
    for (t in 1:N) {
      incidence_strat[t, j] ~ neg_binomial_2(infections[t, j], phi);
    }
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
