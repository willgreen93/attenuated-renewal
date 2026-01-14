library(lhs)

set.seed(1)

n_samples <- 100

# 1. Latin Hypercube for 11 independent variables
lhs_matrix <- randomLHS(n_samples, 13)
colnames(lhs_matrix) <- c("num_nodes", "assortativity_kernel", "ass_v_param", "spatial_kernel", "transmission_prob", 
                          "isolation_time_mean", "isolation_time_sd", "directionality", 
                          "incubation_period", "rec_daily", "infectious_time",
                          "dd_param", "dd_upper")

sample_simplex <- function(u) {
  # u is a vector of length k-1 in (0,1)
  k <- length(u) + 1
  sticks <- numeric(k)
  prod_left <- 1
  for (i in seq_len(k - 1)) {
    sticks[i] <- u[i] * prod_left
    prod_left <- prod_left * (1 - u[i])
  }
  sticks[k] <- prod_left
  return(sticks)
}

alpha_riv <- rep(1, 3) * 60  # target around 1/3, tight bounds

lhs_riv_raw <- randomLHS(n_samples, 3)

lhs_riv <- t(apply(lhs_riv_raw, 1, function(u) {
  gamma_draws <- qgamma(u, shape = alpha_riv, rate = 1)
  gamma_draws / sum(gamma_draws)
}))

colnames(lhs_riv) <- c("riv1", "riv2", "riv3")

target <- c(0.30, 0.45, 0.10, 0.10, 0.03)
kappa <- 100  # total concentration: higher = tighter around mean
alpha <- target * kappa

lhs_raw <- randomLHS(n_samples, 5)  # 5 independent uniform [0,1] draws

lhs_main <- t(apply(lhs_raw, 1, function(u) {
  gamma_draws <- qgamma(u, shape = alpha, rate = 1)
  gamma_draws / sum(gamma_draws)
}))

colnames(lhs_main) <- c("p0", "p1", "p2", "p3", "p4")

# 4. Combine all
lhs_full <- cbind(lhs_matrix, lhs_riv, lhs_main)

# 5. Transform variables to their real-world ranges (example ranges)
lhs <- as.data.frame(lhs_full) %>%
  dplyr::mutate(
    num_nodes = round(scales::rescale(num_nodes, to = c(25000, 50000))),
    assortativity_kernel = scales::rescale(assortativity_kernel, to = c(-0.5, 0.5)),
    ass_v_param = scales::rescale(ass_v_param, to = c(0.5, 0.85)),
    spatial_kernel = scales::rescale(spatial_kernel, to = c(0, 3)),
    transmission_prob = scales::rescale(transmission_prob, to = c(0.3, 0.7)),
    isolation_time_mean = round(scales::rescale(isolation_time_mean, to = c(45))),
    isolation_time_sd = scales::rescale(isolation_time_sd, to = c(0)),
    directionality = scales::rescale(directionality, to = c(1, 5)),
    incubation_period = round(scales::rescale(incubation_period, to = c(5, 14))),
    rec_daily = scales::rescale(rec_daily, to = c(0.05, 0.2)),
    infectious_time = round(scales::rescale(infectious_time, to = c(15, 45))),
    dd_param = scales::rescale(dd_param, to = c(-2, -1.6)),
    dd_upper = round(scales::rescale(dd_upper, to = c(800, 1200)))
  )

# lhs homogeneous
n_samples <- 100

# 1. Latin Hypercube for 11 independent variables
n_samples_homo <- 30
lhs_matrix_homo <- randomLHS(n_samples_homo, 13)
colnames(lhs_matrix_homo) <- c("num_nodes", "assortativity_kernel", "ass_v_param", "spatial_kernel", "transmission_prob", 
                               "isolation_time_mean", "isolation_time_sd", "directionality", 
                               "incubation_period", "rec_daily", "infectious_time",
                               "dd_param", "dd_upper")

alpha_riv_homo <- c(0,0,1) * 60  # target around 1/3, tight bounds

lhs_riv_raw_homo <- randomLHS(n_samples_homo, 3)

lhs_riv_homo <- t(apply(lhs_riv_raw_homo, 1, function(u) {
  gamma_draws <- qgamma(u, shape = alpha_riv_homo, rate = 1)
  gamma_draws / sum(gamma_draws)
}))

colnames(lhs_riv_homo) <- c("riv1", "riv2", "riv3")

target_homo <- c(1, 0, 0, 0, 0)
kappa_homo <- 100  # total concentration: higher = tighter around mean
alpha_homo <- target_homo * kappa_homo

lhs_raw_homo <- randomLHS(n_samples_homo, 5)  # 5 independent uniform [0,1] draws

lhs_main_homo <- t(apply(lhs_raw_homo, 1, function(u) {
  gamma_draws <- qgamma(u, shape = alpha_homo, rate = 1)
  gamma_draws / sum(gamma_draws)
}))

colnames(lhs_main_homo) <- c("p0", "p1", "p2", "p3", "p4")

# 4. Combine all
lhs_full_homo <- cbind(lhs_matrix_homo, lhs_riv_homo, lhs_main_homo)

lhs_homo <- as.data.frame(lhs_full_homo) %>%
  dplyr::mutate(
    num_nodes = round(scales::rescale(num_nodes, to = c(5000, 10000))),
    assortativity_kernel = scales::rescale(assortativity_kernel, to = c(0)),
    ass_v_param = scales::rescale(ass_v_param, to = c(0.67)),
    spatial_kernel = scales::rescale(spatial_kernel, to = c(0)),
    transmission_prob = scales::rescale(transmission_prob, to = c(0.3, 0.7)),
    isolation_time_mean = round(scales::rescale(isolation_time_mean, to = c(45))),
    isolation_time_sd = scales::rescale(isolation_time_sd, to = c(0)),
    directionality = scales::rescale(directionality, to = c(1)),
    incubation_period = round(scales::rescale(incubation_period, to = c(5, 14))),
    rec_daily = scales::rescale(rec_daily, to = c(0.1)),
    infectious_time = round(scales::rescale(infectious_time, to = c(15, 45))),
    dd_param = scales::rescale(dd_param, to = c(-2)),
    dd_upper = round(scales::rescale(dd_upper, to = c(100, 200)))
  )

### LSC for hetero

n_samples_hetero <- 30
lhs_matrix_hetero <- randomLHS(n_samples_hetero, 13)
colnames(lhs_matrix_hetero) <- c("num_nodes", "assortativity_kernel", "ass_v_param", "spatial_kernel", "transmission_prob", 
                                 "isolation_time_mean", "isolation_time_sd", "directionality", 
                                 "incubation_period", "rec_daily", "infectious_time",
                                 "dd_param", "dd_upper")

alpha_riv_hetero <- c(0,0,1) * 60  # target around 1/3, tight bounds

lhs_riv_raw_hetero <- randomLHS(n_samples_hetero, 3)

lhs_riv_hetero <- t(apply(lhs_riv_raw_hetero, 1, function(u) {
  gamma_draws <- qgamma(u, shape = alpha_riv_hetero, rate = 1)
  gamma_draws / sum(gamma_draws)
}))

colnames(lhs_riv_hetero) <- c("riv1", "riv2", "riv3")

target_hetero <- c(1, 0, 0, 0, 0)
kappa_hetero <- 100  # total concentration: higher = tighter around mean
alpha_hetero <- target_hetero * kappa_hetero

lhs_raw_hetero <- randomLHS(n_samples_hetero, 5)  # 5 independent uniform [0,1] draws

lhs_main_hetero <- t(apply(lhs_raw_hetero, 1, function(u) {
  gamma_draws <- qgamma(u, shape = alpha_hetero, rate = 1)
  gamma_draws / sum(gamma_draws)
}))

colnames(lhs_main_hetero) <- c("p0", "p1", "p2", "p3", "p4")

# 4. Combine all
lhs_full_hetero <- cbind(lhs_matrix_hetero, lhs_riv_hetero, lhs_main_hetero)

lhs_hetero <- as.data.frame(lhs_full_hetero) %>%
  dplyr::mutate(
    num_nodes = round(scales::rescale(num_nodes, to = c(25000, 50000))),
    assortativity_kernel = scales::rescale(assortativity_kernel, to = c(0)),
    ass_v_param = scales::rescale(ass_v_param, to = c(0.67)),
    spatial_kernel = scales::rescale(spatial_kernel, to = c(0)),
    transmission_prob = scales::rescale(transmission_prob, to = c(0.3, 0.7)),
    isolation_time_mean = round(scales::rescale(isolation_time_mean, to = c(45))),
    isolation_time_sd = scales::rescale(isolation_time_sd, to = c(0)),
    directionality = scales::rescale(directionality, to = c(1, 5)),
    incubation_period = round(scales::rescale(incubation_period, to = c(5, 14))),
    rec_daily = scales::rescale(rec_daily, to = c(0.1)),
    infectious_time = round(scales::rescale(infectious_time, to = c(15, 45))),
    dd_param = scales::rescale(dd_param, to = c(-1.6,-2)),
    dd_upper = round(scales::rescale(dd_upper, to = c(800, 1200)))
  )
