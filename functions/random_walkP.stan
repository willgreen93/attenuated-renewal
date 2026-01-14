data {
  int<lower=1> k;                 // number of activity groups
  int<lower=0> N;                 // number of observed days
  int<lower=0> h;                 // number of future days to project
  int<lower=0> incidence_strat[N, k];   // observed cases by group and day
}

parameters {
  vector<lower=0>[k] sigma;       
  real<lower=0> infections[N + h, k];
  real<lower=0> phi;
}

model {
  // Priors for process noise
  //sigma ~ normal(0, 5);        // weakly informative prior
  sigma ~ exponential(0.25);        // weakly informative prior
  phi ~ lognormal(log(30), 0.5);
  
  for (j in 1:k) {
    infections[1, j] ~ normal(incidence_strat[1, j], sigma[j]) T[0,];
  }
  
  for (j in 1:k) {
    for (t in 2:(N + h)) {
      infections[t, j] ~ normal(infections[t-1, j], sigma[j]) T[0,];
    }
  }
  
  for (j in 1:k) {
    for (t in 1:N) {
      incidence_strat[t, j] ~ neg_binomial_2(infections[t, j], phi);
    }
  }
}

generated quantities {
  real forecast[h, k];
  for (j in 1:k) {
    for (t in 1:h) {
      forecast[t, j] = neg_binomial_2_rng(infections[N + t, j], phi);
    }
  }
}
