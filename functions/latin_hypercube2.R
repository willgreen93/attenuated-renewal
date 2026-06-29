library(lhs)

set.seed(1)

n_samples <- 100

# 1. Latin Hypercube for 11 independent variables
lhs_matrix <- randomLHS(n_samples, 12)
colnames(lhs_matrix) <- c("num_nodes", "assortativity_kernel", "ass_v_param", "spatial_kernel", "transmission_prob", "directionality", 
                          "incubation_period", "rec_daily", "infectious_time", "dd_param", "dd_upper", "vers")

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

target <- c(0.30, 0.45, 0.10, 0.10, 0.03) #https://pmc.ncbi.nlm.nih.gov/articles/PMC7089812/#SM1
kappa <- 100  # total concentration: higher = tighter around mean
alpha <- target * kappa

lhs_raw <- randomLHS(n_samples, 5)  # 5 independent uniform [0,1] draws

lhs_main <- t(apply(lhs_raw, 1, function(u) {
  gamma_draws <- qgamma(u, shape = alpha, rate = 1)
  gamma_draws / sum(gamma_draws)
}))

colnames(lhs_main) <- c("p0", "p1", "p2", "p3", "p4")

# 4. Combine all
lhs_full <- cbind(lhs_matrix, lhs_main)

# 5. Transform variables to their real-world ranges (example ranges)
lhs <- as.data.frame(lhs_full) %>%
  dplyr::mutate(
    tag=1:nrow(lhs_full),
    num_nodes = round(scales::rescale(num_nodes, to = c(10000, 50000))),
    assortativity_kernel = scales::rescale(assortativity_kernel, to = c(-0.5, 1)),
    ass_v_param = scales::rescale(ass_v_param, to = c(-0.5, 1)),
    spatial_kernel = scales::rescale(spatial_kernel, to = c(0, 3)),
    vers = scales::rescale(vers, to = c(0.4, 0.6)),
    transmission_prob = scales::rescale(transmission_prob, to = c(0.2, 0.5)),
    directionality = scales::rescale(directionality, to = c(1, 5)),
    incubation_period = round(scales::rescale(incubation_period, to = c(5, 14))),
    rec_daily = scales::rescale(rec_daily, to = c(0.05, 0.2)),
    infectious_time = round(scales::rescale(infectious_time, to = c(15, 45))),
    dd_param = scales::rescale(dd_param, to = c(-2, -1.6)),
    dd_upper = round(scales::rescale(dd_upper, to = c(100, 500))), 
    folder="MSM_like"
  )

# HETEOR AND HOMO
lhs_hetero_raw <- randomLHS(n=30, 7)
colnames(lhs_hetero_raw) <- c("num_nodes", "transmission_prob", "incubation_period", "rec_daily", "infectious_time", "dd_param", "dd_upper")

lhs_hetero <- as.data.frame(lhs_hetero_raw) %>%
  dplyr::mutate(
    tag=1:nrow(lhs_hetero_raw),
    num_nodes = round(scales::rescale(num_nodes, to = c(10000, 50000))),
    transmission_prob = scales::rescale(transmission_prob, to = c(0.2, 0.5)),
    incubation_period = round(scales::rescale(incubation_period, to = c(5, 14))),
    infectious_time = round(scales::rescale(infectious_time, to = c(15, 45))),
    dd_param = scales::rescale(dd_param, to = c(-2, -1.6)),
    dd_upper = round(scales::rescale(dd_upper, to = c(100, 500))), 
    folder="hetero")

lhs_homo   <- lhs_hetero %>% mutate(folder="homo", dd_param=NA, dd_upper=scales::rescale(dd_upper, to=c(40,100)), num_nodes=round(scales::rescale(num_nodes, to=c(5000,10000)), 0)) 

lhs_hh <- rbind(lhs_hetero, lhs_homo)
