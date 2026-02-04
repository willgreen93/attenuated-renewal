library(dplyr)
library(purrr)
library(ggplot2)
library(igraph)
library(EpiEstim)
library(profvis)
library(truncnorm)

word_to_num <- c("one" = 1, "two" = 2, "three" = 3, "four" = 4)

main_ot_partners_generator2 <- function(type="hetero", num_nodes=1000, dd_upper=100, dd_param=-1.81, main_partners_prop=c(zero = 0.3, one = 0.45, two = 0.125, three = 0.098, four = 0.027), monogamy=0) {
  if(type %in% c("hetero", "sim")) desired_degrees <- sample(1:dd_upper, num_nodes, replace = TRUE, prob = (1:dd_upper)^dd_param)
  if(type=="homo") desired_degrees <- unname(rep(dd_upper, num_nodes))
  
  if (monogamy != 0) {
    one_partner_indices <- which(desired_degrees == 1)
    indices_to_remove <- sample(one_partner_indices, floor(monogamy * num_nodes * main_partners_prop["one"]))
    desired_degrees <- desired_degrees[-indices_to_remove]
    num_nodes <- length(desired_degrees)
    
    main_partners_prop_hold <- main_partners_prop
    denom <- 1 - monogamy * main_partners_prop_hold["one"]
    main_partners_prop["one"] <- ((1 - monogamy) * main_partners_prop_hold["one"]) / denom
    main_partners_prop[c("zero", "two", "three", "four")] <- main_partners_prop_hold[c("zero", "two", "three", "four")] / denom
  }
  
  main_partners <- rep(0, num_nodes)
  assigned <- rep(FALSE, num_nodes)
  
  for (k in rev(names(main_partners_prop))) {
    eligible <- which(desired_degrees >= word_to_num[k] & !assigned)
    if (length(eligible) == 0) next
    
    pk <- main_partners_prop[as.character(k)] * num_nodes / length(eligible)
    selected <- rbinom(length(eligible), size = 1, prob = pk) == 1
    main_partners[eligible[selected]] <- word_to_num[k]
    assigned[eligible[selected]] <- TRUE
  }
  
  ot_partners <- desired_degrees - main_partners
  
  tot_partners <- main_partners + ot_partners
  order_idx <- order(tot_partners, decreasing = TRUE)
  
  return(list(main_partners = main_partners[order_idx],
              ot_partners = ot_partners[order_idx],
              tot_partners = tot_partners[order_idx]))
}

generate_coordinates <- function(n) list(x = runif(n, 0, 1), y = runif(n, 0, 1))

compute_weights <- function(tot_degree_vec, i, degrees_left, assortativity_kernel, role_vec, role_matrix, x, y, spatial_kernel, adj_list_ot=NULL, nodes) {
  node_whole_network <- as.character(nodes[i])
  role <- role_vec[i]
  
  num_nodes <- length(tot_degree_vec)
  
  target_indicies <- (i+1):num_nodes
  
  ot_vec <- rep(1L,length(target_indicies))
  
  if(!is.null(adj_list_ot)){
    ot_node_connections <- unique(match(adj_list_ot[[node_whole_network]], nodes[target_indicies])) # indexes of connections in the main network
    ot_vec[ot_node_connections] <- 0
  }
  
  assort_weights <- exp(-assortativity_kernel * abs(log(tot_degree_vec[i]) - log(tot_degree_vec[target_indicies])))
  role_weights <- role_matrix[role_vec[i], role_vec[target_indicies]]
  spatial_weights <- exp(-spatial_kernel * sqrt((x[i] - x[target_indicies])^2 + (y[i] - y[target_indicies])^2))
  
  out <- degrees_left[target_indicies] * assort_weights * role_weights * spatial_weights * tot_degree_vec[target_indicies] * ot_vec
  candidates <- 1:(num_nodes-i)
  
  weights <- out/sum(out)
  candidates <- c((i+1):num_nodes)
  
  #[sample(candidates, size = degrees_left[i], prob = weights[candidates], replace = FALSE)]
  
  if(length(candidates)==1) contacts <- candidates
  else contacts <- sample(candidates, size = degrees_left[i], prob = weights, replace = FALSE)
  
  return(contacts)
}

compute_weights2 <- function(tot_degree_vec, i, degrees_left, assortativity_kernel, role_vec, role_matrix, x, y, spatial_kernel, adj_list_ot=NULL, nodes) {
  node_whole_network <- as.character(nodes[i])
  role <- role_vec[i]
  
  num_nodes <- length(tot_degree_vec)
  
  target_indicies <- (i+1):num_nodes
  
  if(role==3) compatible_indicies <- target_indicies
  if(role==1) compatible_indicies <- target_indicies[role_vec[target_indicies] %in% c(2,3)]
  if(role==2) compatible_indicies <- target_indicies[role_vec[target_indicies] %in% c(1,3)]
  
  indicies_degrees_available <- target_indicies[degrees_left[target_indicies] != 0]
  
  non_ot_indicies <- target_indicies
  
  if(!is.null(adj_list_ot)){
    ot_node_connections <- unique(match(adj_list_ot[[node_whole_network]], nodes[target_indicies])) # indexes of connections in the main network
    non_ot_indicies <- target_indicies[!target_indicies %in% ot_node_connections]
  }
  
  candidates <- intersect(intersect(compatible_indicies, indicies_degrees_available), non_ot_indicies)
  
  assort_weights <- exp(-assortativity_kernel * abs(log(tot_degree_vec[i]) - log(tot_degree_vec[candidates])))
  role_weights <- role_matrix[role_vec[i], role_vec[candidates]]
  spatial_weights <- exp(-spatial_kernel * sqrt((x[i] - x[candidates])^2 + (y[i] - y[candidates])^2))
  
  out <- assort_weights * role_weights * spatial_weights * tot_degree_vec[candidates]
  
  weights <- out/sum(out)
  
  if(length(candidates)==1) contacts <- candidates
  else contacts <- sample(candidates, size = degrees_left[i], prob = weights, replace = FALSE)
  
  return(contacts)
}

compute_weights2_fast <- function(tot_degree_vec, i, degrees_left, assortativity_kernel, role_vec, role_matrix, x, y, spatial_kernel, adj_list_ot = NULL, nodes) {
  node_name <- as.character(nodes[i])
  role_i <- role_vec[i]
  deg_i <- tot_degree_vec[i]
  
  num_nodes <- length(tot_degree_vec)
  target_indices <- (i + 1):num_nodes
  if (length(target_indices) == 0) return(integer(0))  # no candidates left
  
  # Role compatibility mask
  role_targets <- role_vec[target_indices]
  if (role_i == 1) role_mask <- role_targets %in% c(2, 3)
  else if (role_i == 2) role_mask <- role_targets %in% c(1, 3)
  else role_mask <- rep(TRUE, length(target_indices))
  
  # Degree availability mask
  degree_mask <- degrees_left[target_indices] > 0
  
  # OT exclusion mask
  if (!is.null(adj_list_ot)) {
    ot_neighbors <- adj_list_ot[[node_name]]
    ot_mask <- !nodes[target_indices] %in% ot_neighbors
  } else {
    ot_mask <- rep(TRUE, length(target_indices))
  }
  
  valid_mask <- role_mask & degree_mask & ot_mask
  if (!any(valid_mask)) return(integer(0))
  
  candidates <- target_indices[valid_mask]
  dx <- x[i] - x[candidates]
  dy <- y[i] - y[candidates]
  
  deg_cand <- tot_degree_vec[candidates]
  role_cand <- role_vec[candidates]
  
  assort_weights <- exp(-assortativity_kernel * abs(log(deg_i) - log(deg_cand)))
  role_weights <- role_matrix[role_i, role_cand]
  spatial_weights <- exp(-spatial_kernel * sqrt(dx^2 + dy^2))
  
  weights <- assort_weights * role_weights * spatial_weights * deg_cand
  weights <- weights / sum(weights)
  
  num_contacts <- min(degrees_left[i], length(candidates))
  
  if(length(candidates)==1) return(candidates)
  else return(sample(candidates, size = num_contacts, prob = weights, replace = FALSE))
}

full_adj_list <- function(adj_list_named, node_set, num_nodes) {
  full <- vector("list", num_nodes)
  for (node in node_set) {
    idx <- as.character(node)
    if (!is.null(adj_list_named[[idx]])) {
      full[[node]] <- adj_list_named[[idx]]
    }
  }
  return(full)
}

build_graph_adj_list <- function(num_nodes_g, nodes, degrees_raw, tot_degree_vec, role_vec, role_matrix, x, y, assortativity_kernel, spatial_kernel, adj_list_ot=NULL, cm) {
  set.seed(1)
  adj_list <- vector("list", length = num_nodes_g)
  degrees <- degrees_raw
  
  for (i in seq_len(num_nodes_g)) {
    if(i %% 1000 == 0) print(i)
    if(degrees[i] == 0) next
    #if(cm == 1) contacts <- compute_weights(tot_degree_vec, i, degrees_left=degrees, assortativity_kernel, role_vec, role_matrix, x, y, spatial_kernel, adj_list_ot=adj_list_ot, nodes)
    #if(cm == 2) contacts <- compute_weights2(tot_degree_vec, i, degrees_left=degrees, assortativity_kernel, role_vec, role_matrix, x, y, spatial_kernel, adj_list_ot=adj_list_ot, nodes)
    contacts <- compute_weights2_fast(tot_degree_vec, i, degrees_left=degrees, assortativity_kernel, role_vec, role_matrix, x, y, spatial_kernel, adj_list_ot=adj_list_ot, nodes)
    
    adj_list[[i]] <- sort(unique(c(adj_list[[i]], contacts)))
    for (j in contacts) {
      adj_list[[j]] <- unique(c(adj_list[[j]], i))
    }
    
    degrees[i] <- 0
    degrees[contacts] <- pmax(0, degrees[contacts] - 1)
    if(sum(degrees) < 0.005*sum(degrees_raw)) break
    #print(adj_list[[i]])
  }
  
  names(adj_list) <- nodes
  adj_list_final <- lapply(adj_list, function(idxs) nodes[idxs])
  
  return(adj_list_final)
}

assortative_graph_generator2 <- function(type, num_nodes, riv, assortativity_kernel, spatial_kernel, main_partners_prop, ass_v_param, dd_param, dd_upper, matrix_tag, seed_add, role, monogamy, cm, ot_max=365) {
  set.seed(seed_add)
  partners <- main_ot_partners_generator2(type, num_nodes, dd_upper, dd_param, main_partners_prop, monogamy)
  
  main_deg <- partners$main_partners
  ot_deg <- partners$ot_partners
  tot_deg <- partners$tot_partners
  
  if(sum(ot_deg)%%2 != 0) ot_deg <- ot_deg + replace(rep(0, length(ot_deg)), sample(length(ot_deg), 1), 1)
  if(sum(main_deg)%%2 != 0) main_deg <- main_deg + replace(rep(0, length(main_deg)), sample(length(main_deg), 1), 1)
  
  role_vector <- if (role) rep(c(3, 2, 1), length.out = num_nodes) else rep(3, num_nodes)
  
  role_matrix_base <- matrix(c(0, 1/3, 0, 1/3, 0, 0, 0, 0, 1/3), nrow = 3)
  role_matrix_adj <- matrix(c(0, -1/6, 1/6, -1/6, 0, 1/6, 1/6, 1/6, -1/3), nrow = 3)
  role_matrix <- role_matrix_base + ass_v_param * role_matrix_adj
  
  # Identify nodes for each network
  main_nodes <- which(main_deg > 0)
  ot_nodes <- which(ot_deg > 0)
  
  coords <- generate_coordinates(num_nodes)
  coords_main <- list(x=coords$x[main_nodes], y=coords$y[main_nodes])
  coords_ot <- list(x=coords$x[ot_nodes], y=coords$y[ot_nodes])
  
  # Build adj list for OT partners
  adj_list_ot <- build_graph_adj_list(num_nodes_g=length(ot_nodes), nodes=ot_nodes, degrees_raw=ot_deg[ot_deg > 0], tot_degree_vec=tot_deg[ot_nodes], 
                                      role_vec=role_vector[ot_nodes], role_matrix=role_matrix, x=coords_ot$x, y=coords_ot$y, 
                                      assortativity_kernel=assortativity_kernel, spatial_kernel=spatial_kernel, adj_list_ot=NULL, cm=cm)
  
  #for(i in 1:length(ot_nodes)) if(length(adj_list_ot[[as.character(i)]]) != ot_deg[i]) print(i)
  
  # Build adj list for main partners
  adj_list_main <- build_graph_adj_list(num_nodes_g=length(main_nodes), nodes=main_nodes, degrees_raw=main_deg[main_deg > 0], tot_degree_vec=tot_deg[main_nodes], 
                                        role_vec=role_vector[main_nodes], role_matrix=role_matrix, coords_main$x, coords_main$y, 
                                        assortativity_kernel=assortativity_kernel, spatial_kernel=spatial_kernel, adj_list_ot=adj_list_ot, cm=cm)
  
  #for(i in 1:length(main_nodes)) if(length(adj_list_main[[as.character(i)]]) != main_deg[i]) print(i)
  
  
  # Build full adjacency lists
  list_main <- full_adj_list(adj_list_main, main_nodes, num_nodes)
  list_ot   <- full_adj_list(adj_list_ot, ot_nodes, num_nodes)
  
  g_main <- igraph::simplify(igraph::graph_from_adj_list(list_main, mode = "all"))
  g_ot   <- igraph::simplify(igraph::graph_from_adj_list(list_ot, mode = "all"))
  
  # Optionally add edge attributes
  E(g_main)$source <- "main"
  E(g_ot)$source   <- "ot"
  E(g_ot)$day <- sample(1:ot_max, ecount(g_ot), replace = TRUE)
  
  # Combine graphs
  G <- igraph::union(g_main, g_ot)
  
  output <- list(num_nodes = num_nodes, g_main = g_main, g_ot = g_ot, G = G, tot_partners = tot_deg, main_partners = main_deg, ot_partners = ot_deg, riv = riv, role_vector = role_vector, role = role, monogamy=monogamy,
                 ass_v_param = ass_v_param, assortativity_kernel = assortativity_kernel, spatial_kernel = spatial_kernel, dd_param = dd_param, dd_upper = dd_upper, matrix_tag = matrix_tag,
                 main_partners_prop = main_partners_prop)
  
  saveRDS(output, file = paste0("mixing_matricies/G_", num_nodes, "_assortativity_kernel=", assortativity_kernel, "_spatial_kernel=", spatial_kernel, "_role=", role, "_role_ass=", round(ass_v_param, 2),
                                "_", matrix_tag, "_dd=", dd_param, "_d_lim=", dd_upper, "_monogamy=", monogamy, ".rds"))
  
  return(output)
}

assortative_graph_generator3 <- function(type, num_nodes, riv, assortativity_kernel, spatial_kernel, main_partners_prop, ass_v_param, dd_param, dd_upper, matrix_tag, seed_add, role, monogamy, cm, ot_max=365) {
  set.seed(seed_add)
  partners <- main_ot_partners_generator2(type, num_nodes, dd_upper, dd_param, main_partners_prop, monogamy=0)
  
  main_deg <- partners$main_partners
  ot_deg <- partners$ot_partners
  tot_deg <- partners$tot_partners
  
  if(sum(ot_deg)%%2 != 0) ot_deg <- ot_deg + replace(rep(0, length(ot_deg)), sample(length(ot_deg), 1), 1)
  if(sum(main_deg)%%2 != 0) main_deg <- main_deg + replace(rep(0, length(main_deg)), sample(length(main_deg), 1), 1)
  
  role_vector <- if (role) sample(c(1, 2, 3), size = num_nodes, replace = TRUE, prob = riv) else rep(3, num_nodes)
  
  role_matrix_base <- matrix(c(0, 1/3, 0, 1/3, 0, 0, 0, 0, 1/3), nrow = 3)
  role_matrix_adj <- matrix(c(0, -1/6, 1/6, -1/6, 0, 1/6, 1/6, 1/6, -1/3), nrow = 3)
  role_matrix <- role_matrix_base + ass_v_param * role_matrix_adj
  
  # Identify nodes for each network
  main_nodes <- which(main_deg > 0)
  ot_nodes <- which(ot_deg > 0)
  
  coords <- generate_coordinates(num_nodes)
  coords_main <- list(x=coords$x[main_nodes], y=coords$y[main_nodes])
  coords_ot <- list(x=coords$x[ot_nodes], y=coords$y[ot_nodes])
  
  # Build adj list for OT partners
  adj_list_ot <- build_graph_adj_list(num_nodes_g=length(ot_nodes), nodes=ot_nodes, degrees_raw=ot_deg[ot_deg > 0], tot_degree_vec=tot_deg[ot_nodes], 
                                      role_vec=role_vector[ot_nodes], role_matrix=role_matrix, x=coords_ot$x, y=coords_ot$y, 
                                      assortativity_kernel=assortativity_kernel, spatial_kernel=spatial_kernel, adj_list_ot=NULL, cm=cm)
  
  #for(i in 1:length(ot_nodes)) if(length(adj_list_ot[[as.character(i)]]) != ot_deg[i]) print(i)
  
  # Build adj list for main partners
  adj_list_main <- build_graph_adj_list(num_nodes_g=length(main_nodes), nodes=main_nodes, degrees_raw=main_deg[main_deg > 0], tot_degree_vec=tot_deg[main_nodes], 
                                        role_vec=role_vector[main_nodes], role_matrix=role_matrix, coords_main$x, coords_main$y, 
                                        assortativity_kernel=assortativity_kernel, spatial_kernel=spatial_kernel, adj_list_ot=adj_list_ot, cm=cm)
  
  #for(i in 1:length(main_nodes)) if(length(adj_list_main[[as.character(i)]]) != main_deg[i]) print(i)
  
  
  # Build full adjacency lists
  list_main <- full_adj_list(adj_list_main, main_nodes, num_nodes)
  list_ot   <- full_adj_list(adj_list_ot, ot_nodes, num_nodes)
  
  g_main <- igraph::simplify(igraph::graph_from_adj_list(list_main, mode = "all"))
  g_ot   <- igraph::simplify(igraph::graph_from_adj_list(list_ot, mode = "all"))
  
  # Optionally add edge attributes
  E(g_main)$source <- "main"
  E(g_ot)$source   <- "ot"
  E(g_ot)$day <- sample(1:ot_max, ecount(g_ot), replace = TRUE)
  
  # Combine graphs
  G <- igraph::union(g_main, g_ot)
  
  output <- list(num_nodes = num_nodes, g_main = g_main, g_ot = g_ot, G = G, tot_partners = tot_deg, main_partners = main_deg, ot_partners = ot_deg, riv = as.vector(table(role_vector)/sum(table(role_vector))), role_vector = role_vector, role = role, monogamy=monogamy,
                 ass_v_param = ass_v_param, assortativity_kernel = assortativity_kernel, spatial_kernel = spatial_kernel, dd_param = dd_param, dd_upper = dd_upper, matrix_tag = matrix_tag,
                 main_partners_prop = main_partners_prop)
  
  saveRDS(output, file = paste0("mixing_matricies/G_", num_nodes, "_assortativity_kernel=", assortativity_kernel, "_spatial_kernel=", spatial_kernel, "_role=", role, "_role_ass=", round(ass_v_param, 2),
                                "_", matrix_tag, "_dd=", dd_param, "_d_lim=", dd_upper, "_monogamy=", monogamy, ".rds"))
  
  return(output)
}


get_adj_matrix_between_nodes <- function(g, infectious_nodes, susceptible_nodes, attr = NULL) {
  # Combine node set for subgraph
  all_nodes <- sort(unique(c(infectious_nodes, susceptible_nodes)))
  subgraph <- igraph::simplify(igraph::induced_subgraph(g, all_nodes))
  
  # Create adjacency matrix (optionally using edge attribute)
  adj_matrix <- igraph::as_adjacency_matrix(subgraph, attr = attr, sparse = FALSE)
  colnames(adj_matrix) <- rownames(adj_matrix) <- all_nodes
  
  # Extract relevant rows and columns
  adj_sub <- adj_matrix[as.character(infectious_nodes), 
                        as.character(susceptible_nodes), 
                        drop = FALSE]
  
  return(adj_sub)
}

make_incidence_matrix <- function(df, t_col = "t_E_start", var, t_max) {
  var_sym <- rlang::ensym(var)
  t_sym <- rlang::sym(t_col)
  
  var_name <- rlang::as_string(var_sym)
  t_name <- t_col
  
  # Build skeleton manually
  skeleton <- expand.grid(
    t_E_start = 1:t_max,
    var_val = 0:max(df[[var_name]], na.rm = TRUE)
  )
  names(skeleton) <- c(t_name, var_name)
  
  # Summarise incidence
  incidence_data <- df %>%
    dplyr::filter(!is.na(t_E)) %>%
    dplyr::group_by(!!t_sym, !!var_sym) %>%
    dplyr::summarise(incidence = dplyr::n(), .groups = "drop")
  
  # Fill in missing cells with 0
  incidence_joined <- dplyr::left_join(skeleton, incidence_data,
                                       by = c(t_name, var_name))
  incidence_joined[is.na(incidence_joined)] <- 0
  
  # Convert to matrix
  stats::xtabs(incidence ~ ., data = incidence_joined) %>% as.matrix()
}

subgraph_generator <- function(graph, infectious_nodes, neighbors_of_infecteds_unlisted_susceptible, attr=NULL){
  subgraph <- (induced_subgraph(graph, c(infectious_nodes, neighbors_of_infecteds_unlisted_susceptible)))
  adj_matrix_raw <- as_adjacency_matrix(subgraph, attr=attr)
  
  colnames(adj_matrix_raw) <- sort(c(infectious_nodes, neighbors_of_infecteds_unlisted_susceptible))
  rownames(adj_matrix_raw) <- sort(c(infectious_nodes, neighbors_of_infecteds_unlisted_susceptible))
  
  adj_matrix <- as.matrix(adj_matrix_raw)[which(colnames(adj_matrix_raw) %in% c(infectious_nodes)), which(colnames(adj_matrix_raw) %in% neighbors_of_infecteds_unlisted_susceptible), drop=FALSE]
  return(adj_matrix)
}

outbreak_simulator_fast4 <- function(timesteps=300, infectious_time=25, initial_exposed=1, transmission_prob=0.5, isolation_time_mean=25, isolation_time_sd=0, directionality=1, incubation_period=5, matrix_tag=1, tag="", monogamy=0, rec_daily=0.1, g){
  set.seed(2)
  
  if(is.null(g)) g <- readRDS(paste0("mixing_matricies/G_", num_nodes, "_infectious_time=", infectious_time, "_assortativity_kernel=", assortativity_kernel, "_spatial_kernel=", spatial_kernel, "_role=", role, "_role_ass=", round(ass_v_param,2), "_", matrix_tag, "_dd=", dd_param, "_d_lim=", dd_upper,"_monogamy=", monogamy, ".rds"))
  
  g_ot <- g$g_ot
  g_main <- g$g_main
  G <- g$G
  num_nodes <- length(V(g$G))
  role_vec <- g$role_vector
  
  param_list <- list(num_nodes = num_nodes, timesteps=timesteps, assortativity = g$assortativity_kernel, spatial = g$spatial_kernel, monogamy = g$monogamy, dd_upper = g$dd_upper, dd_param = g$dd_param,  
                     transmission_prob = transmission_prob, infectious_time = infectious_time, directionality = directionality, ass_v_param = g$ass_v_param, matrix_tag = matrix_tag, incubation_period = incubation_period, tag = tag, 
                     isolation_time_mean=isolation_time_mean, isolation_time_sd=isolation_time_sd, tot_degrees=sum(degree(g$G)), tot_ot_degrees=sum(degree(g$g_ot)), tot_main_degrees=sum(degree(g$g_main)), rec_daily=rec_daily,
                     main0=g$main_partners_prop[1], main1=g$main_partners_prop[2], main2=g$main_partners_prop[3], main3=g$main_partners_prop[4], main4=g$main_partners_prop[5], 
                     receptive=g$riv[1], insertive=g$riv[2], versatile=g$riv[3])
  
  main_degrees <- degree(g_main)
  
  g_main_adj_list <- as_adj_list(g_main)
  g_ot_adj_list <- as_adj_list(g_ot)
  G_adj_list <- as_adj_list(G)
  
  if(isolation_time_sd==0) infectious_time_indiv <- rep(isolation_time_mean, num_nodes)
  else infectious_time_indiv <- floor(rtruncnorm(n=num_nodes, mean=isolation_time_mean, sd=isolation_time_sd, a=0, b=infectious_time+0.5))
  
  degrees <- degree(G)
  
  degrees_annual <- g$tot_partners
  
  V(G)$name <- as.character(seq_len(num_nodes))
  V(G)$state <- "S"  # All nodes start as susceptible
  
  exposed_time <- rep(NA, num_nodes)  # Store the time when a node became exposed (useful for moving to I)
  infected_time <- rep(NA, num_nodes)
  recovered_time <- rep(NA, num_nodes)
  
  if(all(degrees<10)) initial_exposed_nodes <- sample(V(G), size=1)
  else initial_exposed_nodes <- sample(V(G)[which(degrees>10)], size=1, prob=degrees[which(degrees>10)])
    
  V(G)$state[as.numeric(initial_exposed_nodes)] <- "E"  # Set initial infected nodes to exposed
  exposed_neighbors <- sample(neighbors(G, initial_exposed_nodes), initial_exposed-1)
  V(G)$state[exposed_neighbors] <- "E" # Set some neighbors to also have been exposed
  
  exposed_time[as.numeric(c(initial_exposed_nodes, exposed_neighbors))] <- 0
  
  infector <- rep(NA, num_nodes)
  infected_by <- rep(NA, num_nodes)
  infector_role <- rep(NA, num_nodes)
  
  incidence_vec <- rep(NA, timesteps)
  state_matrix <- matrix(nrow=timesteps, ncol=4) %>% magrittr::set_colnames(c("S", "E", "I", "R"))
  
  t_b <- (2*transmission_prob*directionality)/(directionality+1)
  b_t <- (2*transmission_prob)/(directionality+1)
  
  btv_mat <- matrix(c(0,t_b,t_b,b_t,0,b_t,b_t,t_b,transmission_prob), nrow=3, ncol=3)
  colnames(btv_mat) <- c("b", "t", "v")
  rownames(btv_mat) <- c("b", "t", "v")
  
  state_vec <- V(G)$state
  #exposed_time <- V(G)$exposed_time
  #infected_time <- V(G)$infectious_time
  
  p <- proc.time()
  for (t in 1:timesteps){
    #if(t==10) break
    
    newly_infectious_nodes <- which(state_vec == "E" & t - exposed_time >= incubation_period)
    state_vec[newly_infectious_nodes] <- "I"
    infected_time[newly_infectious_nodes] <- t
    
    susceptible_nodes <- which(state_vec == "S")
    infectious_nodes <- which(state_vec == "I")
    
    neighbors_of_infecteds_main <- unique(unlist(g_main_adj_list[infectious_nodes]))
    neighbors_of_infecteds_ot <- unique(
      ends(g_ot, E(g_ot)[day == t], names = FALSE) |> 
        (\(m) m[m[,1] %in% infectious_nodes | m[,2] %in% infectious_nodes, , drop=FALSE])() |> 
        (\(m) c(m[m[,1] %in% infectious_nodes, 2], m[m[,2] %in% infectious_nodes, 1]))()
    )
    
    #neighbors_of_infecteds <- G_adj_list[infectious_nodes] 
    #names(neighbors_of_infecteds) <- infectious_nodes
    neighbors_of_infecteds_unlisted <- unique(c(neighbors_of_infecteds_main, neighbors_of_infecteds_ot))
    neighbors_of_infecteds_unlisted_susceptible <- sort(neighbors_of_infecteds_unlisted[which(state_vec[neighbors_of_infecteds_unlisted]=="S")])
    
    if(length(neighbors_of_infecteds_unlisted_susceptible)!=0){
      neigh_list <- adjacent_vertices(G, neighbors_of_infecteds_unlisted_susceptible, mode="all")
      neigh <- unique(unlist(neigh_list))
      infectious_neighbors_of_susceptibles <- as.vector(neigh[state_vec[neigh] == "I"])
      
      adj_matrix_main <- subgraph_generator(graph=g_main, infectious_nodes=infectious_neighbors_of_susceptibles, neighbors_of_infecteds_unlisted_susceptible=neighbors_of_infecteds_unlisted_susceptible)
      adj_matrix_main_timestep <- adj_matrix_main
      if(is.null(nrow(adj_matrix_main))==FALSE) adj_matrix_main_timestep[adj_matrix_main_timestep==1] = rbinom(n=sum(adj_matrix_main_timestep==1), size=1, prob=rec_daily)
      adj_matrix_main_timestep
      
      adj_matrix_ot <- subgraph_generator(graph=g_ot, infectious_nodes=infectious_neighbors_of_susceptibles, neighbors_of_infecteds_unlisted_susceptible=neighbors_of_infecteds_unlisted_susceptible, attr="day")
      #adj_matrix_ot[adj_matrix_ot==0] <- Inf
      adj_matrix_ot_timestep <- (adj_matrix_ot==t)*1
      
      main_ot_matrix <- matrix(nrow=nrow(adj_matrix_ot), ncol=ncol(adj_matrix_ot))
      main_ot_matrix[as.matrix(adj_matrix_ot!=0)] <- "ot"
      main_ot_matrix[as.matrix(adj_matrix_main!=0)] <- "main"
      rownames(main_ot_matrix) <- infectious_neighbors_of_susceptibles
      colnames(main_ot_matrix) <- neighbors_of_infecteds_unlisted_susceptible
      
      trans_matrix <- btv_mat[g$role_vector[infectious_neighbors_of_susceptibles],g$role_vector[neighbors_of_infecteds_unlisted_susceptible]]
      
      risk_onward_mat <- as.matrix((adj_matrix_ot_timestep + adj_matrix_main_timestep)*trans_matrix)
      #if(length(infectious_nodes)!=1) risk_onward_mat <- as.matrix((adj_matrix_ot_timestep + adj_matrix_main_timestep)*trans_matrix)
      if(ncol(risk_onward_mat)==1) colnames(risk_onward_mat) <- neighbors_of_infecteds_unlisted_susceptible
      
      infection_risk <- 1-matrixStats::colProds(1-as.matrix(risk_onward_mat))
      new_infections <- sort(neighbors_of_infecteds_unlisted_susceptible)[which(runif(length(infection_risk)) < infection_risk)]
      
      state_vec[new_infections] <- "E"
      exposed_time[new_infections] <- t
      incidence_vec[t] <- length(new_infections)
      
      if(sum(state_vec=="I") != 0 & length(new_infections)!=0){
        possible_infecteds_raw <- as.matrix(risk_onward_mat[,which(colnames(risk_onward_mat) %in% new_infections), drop=FALSE])
        possible_infecteds <- t(t(possible_infecteds_raw)/colSums(possible_infecteds_raw))
        
        infector[new_infections] <- infectious_neighbors_of_susceptibles[apply(possible_infecteds, 2, function(prob_col) sample(seq_along(prob_col), size = 1, prob = prob_col))]
        if(NA %in% infector[new_infections]) break
        infected_by[new_infections] <- main_ot_matrix[cbind(as.character(infector[new_infections]), as.character(new_infections))]
        infector_role[new_infections] <- role_vec[infector[new_infections]]
      }  
    } 
    
    recoveries <- which(state_vec == "I" & t-infected_time==infectious_time_indiv)
      
    state_vec[recoveries] <- "R"
    recovered_time[recoveries] <- t
    
    state_matrix[t,] <- as.numeric(table(factor(state_vec, levels = c("S", "E", "I", "R"))))
    
    cat("Timestep:", t, "Infectious:", sum(state_vec == "I"), "Exposed:", sum(state_vec == "E"), "\n")
    
    if((state_matrix[t,"I"]+state_matrix[t,"E"])==0 & t<timesteps){
      state_matrix[(t+1):timesteps,] = matrix(rep(state_matrix[t, ], timesteps - t), nrow = timesteps - t, byrow = TRUE)
      break
    }
    
  }
  
  V(G)$state <- state_vec
  
  (proc.time()-p)["elapsed"]
  
  t_peak <- which(incidence_vec==max(incidence_vec, na.rm=T))[1]
  
  generation_interval_matrix <- data.frame(person=1:num_nodes, infector=infector, infected_by=infected_by, 
                                           role=g$role_vector, infector_role=infector_role,
                                           t_E_start=exposed_time, t_E=exposed_time-t_peak, 
                                           t_I_start=infected_time, t_I=infected_time-t_peak,
                                           t_R=recovered_time-t_peak, 
                                           n_connections=degrees, n_connections_annual=degrees_annual,
                                           n_connections_rec=main_degrees, #n_connections_weighted=degrees_weighted,
                                           t_E_p=exposed_time[infector]-t_peak, t_I_p=infected_time[infector]-t_peak, 
                                           t_R_p=recovered_time[infector]-t_peak, 
                                           n_connections_p=degrees[infector]) %>%
    mutate(generation_time=t_E-t_E_p, serial_interval=t_I-t_I_p) %>%
    mutate(!!!param_list) %>%
    filter(is.na(t_E_start)==FALSE) 
  
  state_matrix[is.na(state_matrix)] <- 0
  state_matrix <- state_matrix %>% as.data.frame() %>% mutate(t=1:timesteps) %>% mutate(!!!param_list) 

  incidence <- data.frame(t_E=1:timesteps, incidence=incidence_vec) %>% 
    mutate(t_E_peak=t_E-which(incidence_vec==max(incidence_vec, na.rm=T))[1], incidence=ifelse(is.na(incidence), 0, incidence)) %>% 
    mutate(!!!param_list) %>%
    mutate(roll_incidence = zoo::rollmean(incidence, 7, fill = NA, align = "right")) 
  
  incidence_summary <- incidence %>% 
    summarise(num_nodes=mean(num_nodes),
              t_peak_roll_incidence = t_E[which.max(roll_incidence)],
              t_passed_5 = t_E[which(roll_incidence > 5)[1]],
              max_roll_incidence = max(roll_incidence, na.rm = TRUE),
              t_5_peak = t_peak_roll_incidence - t_passed_5,
              final_roll_incidence = roll_incidence[which.max(t_E)])
    
  total_incidence <- nrow(generation_interval_matrix)
  
  incidence_roles <- generation_interval_matrix %>% group_by(role) %>% summarise(n=n()) %>%
    tidyr::spread(role, n) %>%
    mutate(total=sum(`1`, `2`, `3`)) %>%
    mutate(p1=`1`/total, p2=`2`/total, p3=`3`/total) %>%
    mutate(proportion_network_infected=total_incidence/num_nodes)
  
  dependent_variables <- cbind(incidence_roles, incidence_summary) %>%
    mutate(!!!param_list)
  
  incidence_matrix_annual    <- make_incidence_matrix(df=generation_interval_matrix, var = "n_connections_annual", t_max = timesteps)
  incidence_matrix_recurrent <- make_incidence_matrix(df=generation_interval_matrix, var = "n_connections_rec",     t_max = timesteps)
  #incidence_matrix_nw        <- make_incidence_matrix(df=generation_interval_matrix, var = "n_connections",        t_max = timesteps)
  
  n_connections_comparison <- generation_interval_matrix %>% group_by(t_E) %>% summarise(n=n(), n_connections_mean = mean(n_connections), n_connections_sd=sd(n_connections), n_connections_se = n_connections_sd/sqrt(n)) %>% mutate(!!!param_list)
  
  output <- list(param_list=param_list, 
                 G=G, 
                 g_ot=g_ot, 
                 g_main=g_main, 
                 degrees_annual=degrees_annual,
                 degrees=as.vector(degrees),
                 degrees_recurrent=main_degrees,
                 incidence=incidence,
                 incidence_matrix_annual=incidence_matrix_annual,
                 #incidence_matrix_nw=incidence_matrix_nw,
                 incidence_matrix_recurrent=incidence_matrix_recurrent,
                 dependent_variables=dependent_variables,
                 role_vector=g$role_vector,
                 initial_exposed_nodes=initial_exposed_nodes, 
                 generation_interval_matrix=generation_interval_matrix, 
                 n_connections_comparison=n_connections_comparison, 
                 state_matrix=state_matrix)
  
  #saveRDS(output, paste0("simulation/epidemic_", num_nodes, "_ass_", g$assortativity_kernel, "_spatial=", g$spatial_kernel ,"_role=", g$role, "_monogamy=", monogamy, "_role_ass=", round(g$ass_v_param,2), 
  #                       "_iso_mean=",isolation_time_mean, "_iso_sd=", isolation_time_sd,"_dd=", g$dd_param, "_d_lim=", g$dd_upper, "_beta_", transmission_prob, "_infectious_time=", infectious_time, "_directionality_", directionality, "_incu_", incubation_period, "_tag_", tag, "_matrix_", matrix_tag, "_5.rds"))
  
  return(output)
  
}

g <- assortative_graph_generator2(num_nodes=100, riv=c(1/3, 1/3, 1/3), assortativity_kernel=0, spatial_kernel=0, main_partners_prop=c(zero = 0, one = 1, two = 0, three = 0, four = 0), ass_v_param=2/3, dd_param=-1.81, dd_upper=15, matrix_tag=1, seed_add=1, role=T, monogamy=0, cm=3)
igraph::plot.igraph(g$G, layout = layout_with_fr, vertex.size = 6, vertex.label = NA,  vertex.color = g$role_vector, vertex.frame.color=NA)
mtext("A", side = 3, line = -1.5, adj = 0, cex = 1, font = 2)

g1 <- assortative_graph_generator2(num_nodes=100, riv=c(1/3, 1/3, 1/3), assortativity_kernel=0.5, spatial_kernel=0, main_partners_prop=c(zero = 1, one = 0, two = 0, three = 0, four = 0), ass_v_param=2/3, dd_param=-1.81, dd_upper=15, matrix_tag=1, seed_add=1, role=T, monogamy=0, cm=3)
igraph::plot.igraph(g1$G, layout = layout_with_fr, vertex.size = 6, vertex.label = NA,  vertex.color = g$role_vector, vertex.frame.color=NA)
mtext("B", side = 3, line = -1.5, adj = 0, cex = 1, font = 2)

g2 <- assortative_graph_generator2(num_nodes=100, riv=c(1/3, 1/3, 1/3), assortativity_kernel=0, spatial_kernel=0, main_partners_prop=c(zero = 1, one = 0, two = 0, three = 0, four = 0), ass_v_param=0.15, dd_param=-1.81, dd_upper=15, matrix_tag=1, seed_add=1, role=T, monogamy=0, cm=3)
igraph::plot.igraph(g2$G, layout = layout_with_fr, vertex.size = 6, vertex.label = NA,  vertex.color = g$role_vector, vertex.frame.color=NA)
mtext("C", side = 3, line = -1.5, adj = 0, cex = 1, font = 2)

dev.off()

#ggsave(filename = "figures/assortative_networks.png", plot = panel, width = 12, height = 4, dpi = 300)

result <- outbreak_simulator_fast4(timesteps=350, infectious_time=25, initial_exposed=1, transmission_prob=0.5, isolation_time_mean=25, isolation_time_sd=0, directionality=1, incubation_period=5, matrix_tag=1, tag="", monogamy=0, rec_daily=0.1, g=A) # Fill with your actual arguments

graph_gen_outbreak_sim2 <- function(type="hetero", num_nodes=1000, riv=c(1/3, 1/3, 1/3), infectious_time=25, assortativity_kernel=0, main_partners_prop=c(zero=0.3, one=0.45, two=0.123, three=0.098, four=0.027), ass_v_param=2/3, dd_param=-1.81, dd_upper=100, matrix_tag=1, spatial_kernel=0, 
                                    timesteps=300, transmission_prob=0.5, isolation_time_mean=30, isolation_time_sd=0, directionality=1, incubation_period=5, rec_daily=0.1, dir="", initial_exposed=5, ot_max=365){
  
  g <- assortative_graph_generator3(type=type, num_nodes, riv, assortativity_kernel, spatial_kernel, main_partners_prop, ass_v_param, dd_param, dd_upper, matrix_tag=1, seed_add=1, role=T, monogamy=0, cm=3, ot_max=ot_max)
  
  outbreak <- outbreak_simulator_fast4(timesteps=timesteps, infectious_time=infectious_time, initial_exposed=initial_exposed, transmission_prob, isolation_time_mean, isolation_time_sd, directionality, 
                                       incubation_period, matrix_tag=1, tag=1, monogamy=0, rec_daily, g=g)
  
  output <- list(g=g, outbreak=outbreak)
  
  saveRDS(output, file=paste0("simulation/",dir,"epidemic_", num_nodes, "_ass=", assortativity_kernel, "_infectious_time=", infectious_time, "_ass_v=", round(ass_v_param,2), 
                              "_dd_param=", dd_param, "_dd_upper=", dd_upper, "_spatial=", spatial_kernel, "_timesteps=", timesteps, "_transmission_prob=", round(transmission_prob,2), "_isol=", isolation_time_mean, "_isol_sd=", isolation_time_sd, 
                              "_directionality=", round(directionality,2), "_incu=", incubation_period, "_rec=", round(rec_daily,2), ".rds"))
  
  return(output)
  
}

i <- 1
g <- assortative_graph_generator3(type="hetero", num_nodes=100, riv=c(lhs$riv1[i], lhs$riv2[i], lhs$riv3[i]), 
                                 assortativity_kernel=lhs$assortativity_kernel[i], spatial_kernel=lhs$spatial_kernel[i], 
                                 main_partners_prop=c(zero = lhs$p0[i], one = lhs$p1[i], two = lhs$p2[i], three = lhs$p3[i], four = lhs$p4[i]), 
                                 ass_v_param = lhs$ass_v_param[i], dd_param = round(lhs$dd_param[i],2), dd_upper = lhs$dd_upper[i],
                                 matrix_tag = 1, seed_add=1, role=T, monogamy=0, cm=3, ot_max=365)

g2 <- assortative_graph_generator3(type="homo", num_nodes=100, riv=c(lhs$riv1[i], lhs$riv2[i], lhs$riv3[i]), 
                                  assortativity_kernel=lhs$assortativity_kernel[i], spatial_kernel=lhs$spatial_kernel[i], 
                                  main_partners_prop=c(zero = lhs$p0[i], one = lhs$p1[i], two = lhs$p2[i], three = lhs$p3[i], four = lhs$p4[i]), 
                                  ass_v_param = lhs$ass_v_param[i], dd_param = round(lhs$dd_param[i],2), dd_upper = 4,
                                  matrix_tag = 1, seed_add=1, role=T, monogamy=0, cm=3, ot_max=365)

g3 <- assortative_graph_generator3(type="sim", num_nodes=100, riv=c(lhs$riv1[i], lhs$riv2[i], lhs$riv3[i]), 
                                   assortativity_kernel=0, spatial_kernel=lhs$spatial_kernel[i], 
                                   main_partners_prop=c(zero = lhs$p0[i], one = lhs$p1[i], two = lhs$p2[i], three = lhs$p3[i], four = lhs$p4[i]), 
                                   ass_v_param = lhs$ass_v_param[i], dd_param = round(lhs$dd_param[i],2), dd_upper = 50,
                                   matrix_tag = 1, seed_add=1, role=T, monogamy=0, cm=3, ot_max=365)

while (dev.cur() > 1) dev.off()

png("figures/homo_hetero_msm.png", width = 1800, height = 600, res = 150)
par(mfrow = c(1, 3), mar = c(1, 1, 3, 1))

igraph::plot.igraph(g2$G, layout = layout_with_fr, vertex.size = 6, vertex.label = NA,  vertex.color = 3, vertex.frame.color=NA)
mtext("Homogeneous", side = 3, line = -1.5, adj = 0, cex = 1, font = 2)

igraph::plot.igraph(g$G, layout = layout_with_fr, vertex.size = 6, vertex.label = NA,  vertex.color = 3, vertex.frame.color=NA)
mtext("Heterogeneous", side = 3, line = -1.5, adj = 0, cex = 1, font = 2)

main_edges <- as_edgelist(g3$g_main)
main_keys  <- apply(main_edges, 1, function(x) paste(sort(x), collapse = "|"))
G_edges <- as_edgelist(g3$G)
G_keys  <- apply(G_edges, 1, function(x) paste(sort(x), collapse = "|"))
is_main <- G_keys %in% main_keys
E(g3$G)$width <- ifelse(is_main, 1.5, 1) 

igraph::plot.igraph(g3$G, layout = L, vertex.size = 6, vertex.label = NA,  vertex.color = g$role_vector, vertex.frame.color=NA, edge.lty = "solid", edge.color="black", edge.width = E(g3$G)$width)
mtext("MSM-like", side = 3, line = -1.5, adj = 0, cex = 1, font = 2)

while (dev.cur() > 1) dev.off()

png("figures/main_ot_msm.png", width = 1800, height = 600, res = 150)
par(mfrow = c(1, 3), mar = c(1, 1, 3, 1))
L <- layout_with_fr(g3$G)

igraph::plot.igraph(g3$g_ot, layout = L, vertex.size = 6, vertex.label = NA,  vertex.color = g$role_vector, vertex.frame.color=NA, edge.lty = "solid", edge.color="black", edge.width = 1)
mtext("One-time network", side = 3, line = -1.5, adj = 0, cex = 1, font = 2)

igraph::plot.igraph(g3$g_main, layout = L, vertex.size = 6, vertex.label = NA,  vertex.color = g$role_vector, vertex.frame.color=NA, edge.lty = "solid", edge.color="black", edge.width = 1.5)
mtext("Recurrent network", side = 3, line = -1.5, adj = 0, cex = 1, font = 2)

igraph::plot.igraph(g3$G, layout = L, vertex.size = 6, vertex.label = NA,  vertex.color = g$role_vector, vertex.frame.color=NA, edge.lty = "solid", edge.color="black", edge.width = E(g3$G)$width)
mtext("MSM-like", side = 3, line = -1.5, adj = 0, cex = 1, font = 2)

dev.off()

ob1 <- graph_gen_outbreak_sim2(num_nodes=10000, riv=c(1/3, 1/3, 1/3), infectious_time=25, assortativity_kernel=0, main_partners_prop=c(zero=0, one=0, two=0, three=0, four=1), ass_v_param=2/3, dd_param=-1.81, dd_upper=100, matrix_tag=1, spatial_kernel=0, timesteps=500, transmission_prob=1, isolation_time_mean=30, isolation_time_sd=0, directionality=1, incubation_period=5, rec_daily=0.1)

p <- proc.time()
for(i in 45:nrow(lhs)){ 
  graph_gen_outbreak_sim2(
    type="hetero", 
    num_nodes = lhs$num_nodes[i],
    riv = c(lhs$riv1[i], lhs$riv2[i], lhs$riv3[i]),
    assortativity_kernel = round(lhs$assortativity_kernel[i],2),
    infectious_time = lhs$infectious_time[i],
    main_partners_prop = c(zero = lhs$p0[i], one = lhs$p1[i], two = lhs$p2[i], three = lhs$p3[i], four = lhs$p4[i]),
    ass_v_param = lhs$ass_v_param[i],
    dd_param = round(lhs$dd_param[i],2),
    dd_upper = lhs$dd_upper[i],
    matrix_tag = 1,
    spatial_kernel = round(lhs$spatial_kernel[i],2),
    timesteps = 500,
    transmission_prob = lhs$transmission_prob[i],
    isolation_time_mean = lhs$isolation_time_mean[i],
    isolation_time_sd = lhs$isolation_time_sd[i],
    directionality = lhs$directionality[i],
    incubation_period = lhs$incubation_period[i],
    rec_daily = lhs$rec_daily[i],
    dir="sim_2/"
  )
  print(paste0(round((proc.time()-p)["elapsed"]/i,1), " seconds"))
}

p <- proc.time()
for(i in 1:nrow(lhs_homo)){ 
  cat("\n")
  print(i)
  graph_gen_outbreak_sim2(
    type="homo",
    num_nodes = lhs_homo$num_nodes[i],
    riv = c(lhs_homo$riv1[i], lhs_homo$riv2[i], lhs_homo$riv3[i]),
    assortativity_kernel = round(lhs_homo$assortativity_kernel[i],2),
    infectious_time = lhs_homo$infectious_time[i],
    main_partners_prop = c(zero = lhs_homo$p0[i], one = lhs_homo$p1[i], two = lhs_homo$p2[i], three = lhs_homo$p3[i], four = lhs_homo$p4[i]),
    ass_v_param = lhs_homo$ass_v_param[i],
    dd_param = round(lhs_homo$dd_param[i],2),
    dd_upper = lhs_homo$dd_upper[i],
    matrix_tag = 1,
    spatial_kernel = round(lhs_homo$spatial_kernel[i],2),
    timesteps = 1000,
    transmission_prob = lhs_homo$transmission_prob[i],
    isolation_time_mean = lhs_homo$isolation_time_mean[i],
    isolation_time_sd = lhs_homo$isolation_time_sd[i],
    directionality = lhs_homo$directionality[i],
    incubation_period = lhs_homo$incubation_period[i],
    rec_daily = lhs_homo$rec_daily[i],
    dir="sim_homo/",
    initial_exposed=3, 
    ot_max=1000
  )
  print(paste0(round((proc.time()-p)["elapsed"]/i,1), " seconds"))
}

p <- proc.time()
for(i in 1:nrow(lhs_hetero)){ 
  cat("\n")
  print(i)
  og1 <- graph_gen_outbreak_sim2(
    type="hetero", 
    num_nodes = lhs_hetero$num_nodes[i],
    riv = c(lhs_hetero$riv1[i], lhs_hetero$riv2[i], lhs_hetero$riv3[i]),
    assortativity_kernel = round(lhs_hetero$assortativity_kernel[i],2),
    infectious_time = lhs_hetero$infectious_time[i],
    main_partners_prop = c(zero = lhs_hetero$p0[i], one = lhs_hetero$p1[i], two = lhs_hetero$p2[i], three = lhs_hetero$p3[i], four = lhs_hetero$p4[i]),
    ass_v_param = lhs_hetero$ass_v_param[i],
    dd_param = round(lhs_hetero$dd_param[i],2),
    dd_upper = lhs_hetero$dd_upper[i],
    matrix_tag = 1,
    spatial_kernel = round(lhs_hetero$spatial_kernel[i],2),
    timesteps = 1000,
    transmission_prob = lhs_hetero$transmission_prob[i],
    isolation_time_mean = lhs_hetero$isolation_time_mean[i],
    isolation_time_sd = lhs_hetero$isolation_time_sd[i],
    directionality = lhs_hetero$directionality[i],
    incubation_period = lhs_hetero$incubation_period[i],
    rec_daily = lhs_hetero$rec_daily[i],
    dir="sim_het/",
    initial_exposed=3
  )
  print(paste0(round((proc.time()-p)["elapsed"]/i,1), " seconds"))
}
