data {
  int<lower=1> k;                      // number of strata
  int<lower=0> N;                      // observed time points
  int<lower=0> h;                      // future time points
  int<lower=0> incidence_strat[N, k]; // observed incidence
}

parameters {
  vector<lower=0>[k] sigma;                   // RW std dev per stratum
  real<lower=0> phi;
  matrix[N + h, k] log_infections;            // latent log-infections
}

transformed parameters {
  matrix<lower=0>[N + h, k] infections;

  for (j in 1:k) {
    for (t in 1:(N + h)) {
      infections[t, j] = exp(log_infections[t, j]);
    }
  }
}

model {
  phi ~ lognormal(log(30), 0.5);
  //sigma ~ normal(0, 1);  // weakly informative
  sigma ~ exponential(sqrt(pi()/2));

  for (j in 1:k) {
    // Prior for initial value
    log_infections[1, j] ~ normal(log1p(incidence_strat[1, j]), sigma[j]);

    // Random walk over log-infections
    for (t in 2:(N + h)) {
      log_infections[t, j] ~ normal(log_infections[t - 1, j], sigma[j]);
    }

    // Observation likelihood
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
