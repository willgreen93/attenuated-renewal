data {
  int<lower=1> N;                      // Number of observed days
  int<lower=1> h;                      // Number of forecast days
  int<lower=1> L;                      // Length of generation time distribution
  int<lower=1> k;                      // Number of strata
  int<lower=0> incidence_strat[N, k];  // Observed incidence per group
  vector<lower=0>[L] gen_weights;      // Generation time distribution
  int<lower=1> S0[k];                  // Initial susceptible population
}

parameters {
  vector<lower=0>[k] R0;               // Initial (upper asymptote) reproduction number per group
  vector<lower=0>[k] delta;            // Steepness kappa per group (repurposed)
  real<lower=1, upper=N> t_change;     // Logistic midpoint t_mid (repurposed)
  real<lower=0> phi;                   // Overdispersion
  vector<lower=0>[k] seed_infections;  // Initial infections
  vector<lower=0>[k] R_floor;          // Lower asymptote per group
}

transformed parameters {
  matrix<lower=0>[N + h, k] infections;
  matrix<lower=0>[N + h, k] cumulative_infections;
  
  for (j in 1:k) {
    infections[1, j] = seed_infections[j];
    cumulative_infections[1, j] = seed_infections[j];
    
    for (t in 2:(N + h)) {
      // Sigmoid (logistic) fall-off for Rt:
        // Rt = R_floor + (R0 - R_floor) * inv_logit( -delta * (t - t_change) )
        real Rt = R_floor[j]
        + (R0[j] - R_floor[j]) * inv_logit( -delta[j] * (t - t_change) );
        
        // Renewal equation
        real infectivity = 0;
        for (s in 1:L) {
          if (t > s) infectivity += infections[t - s, j] * gen_weights[s];
        }
        
        infections[t, j] = fmax(1e-6,
                                Rt * infectivity
                                * (1 - cumulative_infections[t - 1, j] / S0[j]));
        cumulative_infections[t, j] = cumulative_infections[t - 1, j] + infections[t, j];
    }
  }
}

model {
  // Priors (tune as you like)
  R0 ~ lognormal(log(10), 0.4);
  //R_floor  ~ normal(0.5, 0.5);       // lower asymptote; adjust to your context
  R_floor  ~ lognormal(-1.04, 0.83);       // lower asymptote; adjust to your context
  //delta    ~ normal(0.05, 0.0125);     // steepness (per day); larger = sharper drop
  delta ~ lognormal(-3.02, 0.246);
  t_change ~ uniform(1, 365);          // midpoint within observed window
  phi ~ lognormal(log(30), 0.5);
  seed_infections ~ lognormal(log(1), 0.5);

  // Likelihood over observed window
  for (j in 1:k)
    for (t in 1:N)
      incidence_strat[t, j] ~ neg_binomial_2(infections[t, j], phi);
}

generated quantities {
  matrix[h, k] forecast;
  
  for (j in 1:k)
    for (t in 1:h)
      forecast[t, j] = neg_binomial_2_rng(infections[N + t, j], phi);
}
