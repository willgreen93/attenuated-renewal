library(igraph)
library(matrixStats)
library(GGally)
library(ggplot2)

files <- list.files("mixing_matricies2/MSM_like", full.names=T)
info <- file.info(files)
mat_list <- files[order(info$mtime, decreasing=T)]


fit_trunc_powerlaw_alpha <- function(deg, kmax = max(deg), alpha_init = 2) {
  # basic checks
  if (any(deg < 1) || any(deg != floor(deg))) {
    deg <- deg[deg>0]
  }
  if (any(deg > kmax)) {
    stop("All observed degrees must be <= kmax.")
  }
  
  n <- length(deg)
  log_deg_sum <- sum(log(deg))
  k_vals <- 1:kmax
  
  neg_loglik <- function(alpha) {
    # constrain alpha to be positive
    if (alpha <= 0) return(Inf)
    
    log_norm <- log(sum(k_vals^(-alpha)))
    alpha * log_deg_sum + n * log_norm
  }
  
  fit <- optim(
    par = alpha_init,
    fn = neg_loglik,
    method = "L-BFGS-B",
    lower = 1e-6,
    upper = 20
  )
  
  alpha_hat <- fit$par
  
  # approximate SE from observed Hessian if available
  h <- optimHess(alpha_hat, neg_loglik)
  se_hat <- if (is.finite(h[1, 1]) && h[1, 1] > 0) sqrt(1 / h[1, 1]) else NA_real_
  
  list(
    alpha_hat = alpha_hat,
    se = se_hat,
    logLik = -fit$value,
    convergence = fit$convergence,
    message = fit$message
  )
}

compute_versatile_assortativity <- function(g, role_vector, versatile_code = 3) {
  # Attach roles
  V(g)$role <- role_vector
  
  # Get edge endpoints
  edge_ends <- ends(g, E(g), names = FALSE)
  
  r1 <- V(g)$role[edge_ends[, 1]]
  r2 <- V(g)$role[edge_ends[, 2]]
  
  # --- Edge-based version ---
  is_V_edge <- (r1 == versatile_code | r2 == versatile_code)
  is_VV_edge <- (r1 == versatile_code & r2 == versatile_code)
  
  prop_edge_based <- if (sum(is_V_edge) > 0) {
    sum(is_VV_edge) / sum(is_V_edge)
  } else {
    NA_real_
  }
  
  # --- Partner-based version (recommended) ---
  edges_V <- which(is_V_edge)
  
  partners <- ifelse(
    r1[edges_V] == versatile_code,
    r2[edges_V],
    r1[edges_V]
  )
  
  prop_partner_based <- if (length(partners) > 0) {
    mean(partners == versatile_code)
  } else {
    NA_real_
  }
  
  # --- Random mixing expectation ---
  p_V <- mean(role_vector == versatile_code)
  
  ratio_vs_random <- prop_partner_based / p_V
  
  # Output
  list(
    prop_edge_based = prop_edge_based,
    prop_partner_based = prop_partner_based,
    expected_random = p_V,
    ratio_vs_random = ratio_vs_random
  )
}

log_assortativity <- function(g) {
  edgelist <- igraph::as_edgelist(g)
  deg <- igraph::degree(g)
  
  k1 <- log(deg[edgelist[,1]])
  k2 <- log(deg[edgelist[,2]])
  
  cor(k1, k2)
}

compute_edge_distance_stats <- function(g, x, y) {
  # g: igraph object
  # x, y: numeric vectors of node coordinates (same order as V(g))
  
  if (length(x) != vcount(g) || length(y) != vcount(g)) {
    stop("Length of x and y must match number of nodes in graph")
  }
  
  coords <- cbind(x, y)
  
  # edge endpoints
  edge_ends <- ends(g, E(g), names = FALSE)
  
  # coordinates
  x1 <- coords[edge_ends[, 1], 1]
  y1 <- coords[edge_ends[, 1], 2]
  x2 <- coords[edge_ends[, 2], 1]
  y2 <- coords[edge_ends[, 2], 2]
  
  # distances
  edge_dist <- sqrt((x1 - x2)^2 + (y1 - y2)^2)
  
  # return summary
  list(
    mean_edge_distance = mean(edge_dist),
    sd_edge_distance = sd(edge_dist)
  )
}

#fit_trunc_powerlaw_alpha(deg=degree(sim$g$G), sim$g$dd_upper, alpha_init=2)

for(i in 1:length(mat_list)){
  print(i)
  
  sim <- readRDS(mat_list[i])
  graph1 <- igraph::union(sim$g_main, sim$g_ot1)
  
  df_raw <- data.frame(#i=i, 
                       #num_nodes=sim$num_nodes,
                       alpha_param = sim$graph_params$alpha, 
                       alpha_measured=fit_trunc_powerlaw_alpha(deg=degree(graph1), kmax=max(degree(graph1)), alpha_init=2)$alpha_hat,
                       #assortativity_ot=assortativity_degree(sim$g_ot1),
                       #assortativity_main=assortativity_degree(sim$g_main),
                       assortativity_param=sim$graph_params$degree_ass,
                       assortativity_measured=assortativity_degree(graph1),
                       #log_assortativity=log_assortativity(graph1),
                       clustering=transitivity(sim$G, type = "global"),
                       max_degree=max(degree(graph1)),
                       #vers_assortativity2=compute_versatile_assortativity(g=graph1, role_vector=sim$graph_param$roles)$prop_partner_based,
                       vers_assortativity_param=sim$graph_params$role_ass_param,
                       vers_assortativity_measured=compute_versatile_assortativity(g=graph1, role_vector=sim$graph_param$roles)$prop_edge_based,
                       #vers_assortativity2=compute_versatile_assortativity(sim$G, sim$role_vector)$ratio_vs_random,
                       spatial_param=sim$graph_params$spatial_ass,
                       mean_edge_distance=compute_edge_distance_stats(graph1, sim$graph_params$coords[,1], sim$graph_params$coords[,2])[[1]],
                       sd_edge_distance=compute_edge_distance_stats(graph1, sim$graph_param$coords[,1], sim$graph_params$coords[,2])[[2]],
                       biggest_cluster=max(components(sim$G)$csize)/sim$graph_params$n,
                       second_biggest_cluster=sort(components(sim$G)$csize, decreasing=T)[2]/sim$graph_params$n)
  
  if(i==1) df <- df_raw
  else df <- rbind(df, df_raw)
             
}

colQuantiles(as.matrix(df), probs=c(0.025, 0.5, 0.975))

#names <- c(alpha_param="Degree exponent", assortativity_G="Degree assortativity", clustering="Clustering coefficient", max_degree="Maxium degree", vers_assortativity=)

col_names <- c(
  "alpha_param"               = "alpha[input]",
  "alpha_measured"            = "alpha[measured]",
  "assortativity_param"       = "zeta[input]",
  "assortativity_measured"    = "zeta[measured]",
  "vers_assortativity_param"  = "delta[input]",
  "vers_assortativity_measured" = "delta[measured]",
  "spatial_param"             = "epsilon[input]",
  "mean_edge_distance"        = "bar(d)",
  "clustering"                = "C"
)

df_plot <- df[, names(col_names)]
names(df_plot) <- col_names

my_lower <- function(data, mapping) {
  ggplot(data, mapping) +
    geom_point(alpha = 0.5, size = 1.2, colour = "grey40") +
    scale_x_continuous(n.breaks = 3) +
    scale_y_continuous(n.breaks = 3) +
    theme_bw() +
    theme(panel.grid = element_blank())
}

my_upper <- function(data, mapping) {
  x     <- GGally::eval_data_col(data, mapping$x)
  y     <- GGally::eval_data_col(data, mapping$y)
  ct    <- cor.test(x, y)
  corr  <- ct$estimate
  stars <- symnum(ct$p.value, cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                  symbols = c("***", "**", "*", ""))
  col   <- col_numeric(c("steelblue", "white", "firebrick"), domain = c(-1, 1))(corr)
  
  ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label    = paste0(round(corr, 3), stars),
             colour   = col,
             size     = 4,
             fontface = "bold") +
    xlim(0, 1) + ylim(0, 1) +
    theme_void()
}

p_network <- GGally::ggpairs(df_plot,
                             labeller = label_parsed,
                             lower = list(continuous = my_lower),
                             diag  = list(continuous = my_diag),
                             upper = list(continuous = my_upper)) +
  theme(strip.placement   = "outside",
        strip.background  = element_rect(fill = "white", colour = "white"),
        panel.background  = element_rect(fill = "white", colour = "black"),
        panel.grid        = element_blank(),
        strip.text        = element_text(size = 9),
        axis.text.x       = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y       = element_text(size = 7),
        legend.position   = "none")


p_network

ggsave("figures/network_properties.png", p_network, width=10, height=10)

ggplot(df, aes(x=biggest_cluster)) +
  geom_histogram() + 
  theme_bw()
