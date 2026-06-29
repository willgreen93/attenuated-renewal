library(igraph)
library(rstan)
library(ggplot2)
library(dplyr)
library(reshape2)
library(tidyr)
library(zoo)
library(matrixStats)

main_partners_generator <- function(main_partners_prop, degrees_raw){
  set.seed(1)
  n <- length(degrees_raw)
  
  main_degrees_raw <- rep(0, n)
  
  for(k in 4:1){
    eligible <- which(degrees_raw>=k & main_degrees_raw==0)
    n_k <- min(length(eligible), main_partners_prop[paste(k)]*n)
    
    if (n_k > 0) {
      chosen <- sample(eligible, n_k, replace = FALSE)
      main_degrees_raw[chosen] <- k
    }
  }
  
  return(main_degrees_raw)
}

generate_partner_graph <- function(network="ot", all_degrees, degrees_raw, x, y, roles, role_ass, spatial_ass, degree_ass, stop_prop = 0.01, forbidden_adj = NULL, min_graph_day=0){
  set.seed(1+min_graph_day)
  n <- length(degrees_raw)
  degrees <- degrees_raw
  adj_matrix <- as.list(rep(0,n))
  
  while(sum(degrees) > sum(degrees_raw)*stop_prop){
    node_to_assign <- sample(1:n, size=1, prob=degrees^2)
    
    edges_to_assign <- degrees[node_to_assign]
    
    spatial_p <- exp(-spatial_ass*sqrt((x - x[node_to_assign])^2 + (y - y[node_to_assign])^2))
    assortativity_p <- exp(-degree_ass*abs(log(degrees_raw[node_to_assign])-log(all_degrees)))
    role_p <- role_ass[roles[node_to_assign],roles]
    degree_p <- degrees; degree_p[node_to_assign] <- 0
    
    degree_p[adj_matrix[[node_to_assign]]] <- 0
    
    if(network=="main") degree_p[forbidden_adj[[node_to_assign]]] <- 0
    
    probs_raw <- spatial_p * role_p * assortativity_p * degree_p
    probs <- probs_raw/sum(probs_raw)
    if(NaN %in% probs) print(c(roles[node_to_assign],roles[degree_p>0]))
    
    partners <- sample(seq_along(degrees), size = edges_to_assign, replace = FALSE, prob = probs)
    
    adj_matrix[[node_to_assign]] <- c(adj_matrix[[node_to_assign]], partners)
    for(i in partners) adj_matrix[[i]] <- c(adj_matrix[[i]], node_to_assign)
    
    degrees[node_to_assign] <- degrees[node_to_assign]-edges_to_assign
    degrees[partners] <- degrees[partners]-1
    #print(c(network=network, node_to_assign = floor(node_to_assign), degrees = floor(degrees[node_to_assign]), prop_assigned = round(1-sum(degrees)/sum(degrees_raw),2)), quote=F)
    #cat("\n")
  }
  
  adj_list <- lapply(adj_matrix, \(x) x[x != 0])
  g <- graph_from_adj_list(adj_list, mode = "all")
  
  E(g)$source <- network
  E(g)$day <- sample((1:365)+min_graph_day, size=length(E(g)), replace = TRUE)
  
  return(list(network=network, g=g, adj_list=adj_list))
  
}

#g_ot   <- generate_partner_graph(network="ot",   ot_degrees_raw,   x, y, roles, role_ass, spatial_ass, degree_ass, stop_prop = 0.01, forbidden_adj = NULL)
#g_main <- generate_partner_graph(network="main", main_degrees_raw, x, y, roles, role_ass, spatial_ass, degree_ass, stop_prop = 0.01, forbidden_adj = g_ot$adj_list)

MSM_like_graph_generator <- function(n, alpha=NULL, kmax, main_partners_prop, roles_prop, role_ass_param, spatial_ass, degree_ass, tag, stop_prop, folder="MSM_like"){
  print(tag)
  set.seed(1)
  x <- runif(n, 0, 1)
  y <- runif(n, 0, 1)
  roles <- rep(NA, n)
  
  roles_raw <- sample(c(1,2,3), n, replace=T, prob=roles_prop)
  
  if(!is.na(alpha)) degrees_raw <- sample(1:kmax, size = n, replace = TRUE, prob = (1:kmax)^(alpha))
  if(is.na(alpha))  degrees_raw <- round(rep(kmax, n), 0)
  
  roles <- roles_raw
  
  main_degrees_raw <- main_partners_generator(main_partners_prop, degrees_raw)
  ot_degrees_raw <- degrees_raw - main_degrees_raw
  
  role_ass <- matrix(c(0, 2/3, 1/3, 2/3, 0, 1/3, 1/3, 1/3, 1/3), nrow=3) + role_ass_param*matrix(c(0, 1/6, -1/6, 1/6, 0, -1/6, -1/6, -1/6, 1/3), nrow=3)
  
  g_ot1   <- generate_partner_graph(network="ot", all_degrees=degrees_raw, degrees_raw=ot_degrees_raw, x=x, y=y, roles=roles, role_ass=role_ass, spatial_ass=spatial_ass, degree_ass=degree_ass, stop_prop = 0.02, forbidden_adj = NULL)
  g_ot2   <- generate_partner_graph(network="ot", all_degrees=degrees_raw, degrees_raw=ot_degrees_raw, x=x, y=y, roles=roles, role_ass=role_ass, spatial_ass=spatial_ass, degree_ass=degree_ass, stop_prop = 0.02, forbidden_adj = g_ot1$adj_list, min_graph_day=365)
  g_ot <- igraph::union(g_ot1$g, g_ot2$g)
  E(g_ot)$day <- ifelse(is.na(E(g_ot)$day_1), E(g_ot)$day_2, E(g_ot)$day_1)
  
  g_main  <- generate_partner_graph(network="main", all_degrees=degrees_raw, degrees_raw=main_degrees_raw, x, y, roles, role_ass, spatial_ass, degree_ass, stop_prop = 0.02, forbidden_adj = g_ot$adj_list)
  
  G <- igraph::union(g_main$g, g_ot)
  
  params <- list(n=n, class=folder, alpha=alpha, kmax=kmax, roles_prop=roles_prop, role_ass_param=role_ass_param, spatial_ass=spatial_ass, 
  degree_ass=degree_ass, stop_prop=stop_prop, coords=cbind(x,y), roles=roles, degrees_raw=degrees_raw, main_degrees_raw=main_degrees_raw,
  ot_degrees_raw=ot_degrees_raw, role_ass=role_ass)
  
  output <- list(tag=tag, G=G, g_ot=g_ot, g_ot1=g_ot1$g, g_ot2=g_ot2$g, g_main=g_main$g, adj_list_ot=as_adj_list(g_ot), adj_list_main=g_main$adj_list, graph_params=params)
  
  saveRDS(output, paste0("mixing_matricies2/", folder, "/",folder,"_", tag=tag,".rds"))
  
  return(output)

}

hh_graph_generator <- function(n, alpha, kmax, tag, stop_prop, folder){
  if(folder=="homo") alpha <- NA
  output <- MSM_like_graph_generator(n = n, alpha = alpha, kmax = kmax, 
                                     main_partners_prop = c('0'=1, '1'=0, '2'=0, '3'=0, '4'=0), roles_prop = c('1'=0, '2'=0, '3'=1), 
                                     role_ass_param = 0, spatial_ass = 0, degree_ass = 0, tag = tag, stop_prop = 0.02, folder=folder)
  return(output)
}

fill_matrix_from_graph <- function(M, graph, prob) {
  el <- ends(graph, E(graph), names = FALSE)
  
  from <- el[,1]
  to   <- el[,2]
  
  i_to_s <- from %in% rownames(M) & to %in% colnames(M)
  s_to_i <- to   %in% rownames(M) & from %in% colnames(M)
  
  inf_nodes <- c(from[i_to_s], to[s_to_i])
  sus_nodes <- c(to[i_to_s],   from[s_to_i])
  
  index_list <- cbind(match(inf_nodes, rownames(M)),
                      match(sus_nodes, colnames(M)))
  
  M[index_list] <- rbinom(n=nrow(index_list), size=1, prob=prob)
  
  M
}

outbreak_simulator_fast5 <- function(timesteps=300, infectious_period=25, initial_exposed=1, transmission_prob=0.5, directionality=1, incubation_period=5, 
                                     tag="", rec_daily=0.1, g, seed=1){
  set.seed(seed)
  
  G <- g$G
  g_ot <- g$g_ot
  g_main <- g$g_main
  num_nodes <- g$graph_params$n
  
  adj_list_main <- g$adj_list_main
  adj_list_ot   <- g$adj_list_ot
  
  state_vector <- rep("S",num_nodes)
  role_vector <- g$graph_params$roles
  
  exposed_time <- rep(NA, num_nodes)
  infectious_time <- rep(NA, num_nodes)
  recovered_time <- rep(NA, num_nodes)
  infector <- rep(NA, num_nodes)
  
  exposed_vector <- rep(0, timesteps)
  incidence_vector <- rep(0, timesteps)
  
  degrees <- degree(G)
  
  initial_exposed_nodes <- sample(V(G)[which(degrees>50)], size=1, prob=degrees[which(degrees>50)])
  
  state_vector[initial_exposed_nodes] <- "E"  # All nodes start as susceptible
  
  #initial_infected_nodes <- c(13, 37)
  #state_vector[initial_infected_nodes] <- "I"
  
  exposed_time[initial_exposed_nodes] <- 0
  #infectious_time[initial_infected_nodes] <- 0
  
  t_b <- (2*transmission_prob*directionality)/(directionality+1)
  b_t <- (2*transmission_prob)/(directionality+1)
  btv_mat <- matrix(c(0,t_b,t_b,b_t,0,b_t,b_t,t_b,transmission_prob), nrow=3, ncol=3)
  
  for(t in 1:timesteps){
    #if(t==50) break
    #print(c(t, sum(state_vector=="I", na.rm=T)))
    newly_infectious_nodes <- which(state_vector == "E" & t - exposed_time >= incubation_period)
    incidence_vector[t] <- length(newly_infectious_nodes)
    #if(length(newly_infectious_nodes) > 1) break
    state_vector[newly_infectious_nodes] <- "I"
    infectious_time[newly_infectious_nodes] <- t
    #print(infectious_time[newly_infectious_nodes])
    
    if(!("I" %in% state_vector)) next
    
    infectious_nodes <- which(state_vector=="I")
    
    neighbors_main <- sort(unlist(adj_list_main[infectious_nodes]))
    susceptible_neighbors_main <- intersect(neighbors_main, which(state_vector=="S"))
    
    neighbors_ot <- unique(ends(g_ot, E(g_ot)[day == t], names = FALSE) |> 
        (\(m) m[m[,1] %in% infectious_nodes | m[,2] %in% infectious_nodes, , drop=FALSE])() |> 
        (\(m) c(m[m[,1] %in% infectious_nodes, 2], m[m[,2] %in% infectious_nodes, 1]))())
    susceptible_neighbors_ot <- intersect(neighbors_ot, which(state_vector=="S"))
    
    susceptible_neighbors <- sort(unique(c(susceptible_neighbors_ot, susceptible_neighbors_main)))
    if(length(susceptible_neighbors)==0) next
    
    M <- matrix(nrow=length(infectious_nodes), ncol=length(susceptible_neighbors))
    rownames(M) <- infectious_nodes
    colnames(M) <- susceptible_neighbors
    
    M_ot <- fill_matrix_from_graph(M, g_ot, 1)
    M <- fill_matrix_from_graph(M_ot, g_main, rec_daily)
    M[is.na(M)] <- 0
    
    trans_prob <- btv_mat[role_vector[infectious_nodes], role_vector[susceptible_neighbors], drop=F]
    rownames(trans_prob) <- infectious_nodes
    colnames(trans_prob) <- susceptible_neighbors
    
    infection_probability <- trans_prob*M
    
    new_infecteds <- susceptible_neighbors[rbinom(n=length(susceptible_neighbors), size=1, prob=1-colProds(1-infection_probability))==1]
    #print(role_vector[new_infecteds] %>% table())
    
    if(length(new_infecteds)>0) for(m in new_infecteds){
      if(length(infectious_nodes)==1) infector[m] <- infectious_nodes
      else infector[m] <- sample(as.numeric(rownames(infection_probability)), size=1, prob=infection_probability[,paste(m)])
    } 
    
    state_vector[new_infecteds] <- "E"
    exposed_time[new_infecteds] <- t
    
    #print(exposed_time[new_infecteds])
    
    exposed_vector[t] <- length(new_infecteds)
    
    recoveries <- which(state_vector == "I" & t-infectious_time==infectious_period)
    state_vector[recoveries] <- "R"
    recovered_time[recoveries] <- t
    
    if(t %% 50 == 0) cat("Timestep:", t, "Infectious:", sum(state_vector == "I"), "Exposed:", sum(state_vector == "E"), "Recovered:", sum(state_vector=="R"), "\n")
    #cat("\n")
  }
  
  generation_interval_matrix <- data.frame(index=1:num_nodes, 
             role=role_vector,
             exposed_time = exposed_time, 
             infectious_time = infectious_time,
             recovered_time = recovered_time,
             degrees_annual = degree(g$g_ot1),
             degree_recurrent = degree(g_main),
             infector=infector, 
             infector_role=role_vector[as.numeric(infector)],
             infector_exposed=exposed_time[infector],
             infector_infectious=infectious_time[infector],
             generation_interval = exposed_time-exposed_time[infector]) %>%
    filter(!is.na(infector)) %>%
    rowwise() %>%
    mutate(infected_by=ifelse(are_adjacent(g_ot, index, infector), "ot", ifelse(are_adjacent(g_main, index, infector), "main", NA)))
  
  incidence_df <- data.frame(index=1:timesteps,
                             exposed=exposed_vector,
                             incidence=incidence_vector)
  
  inc_role_raw <- generation_interval_matrix %>% group_by(role, infectious_time)             %>% summarise(count=n()) %>% mutate(role=factor(role)) 
  
  inc_role     <- inc_role_raw %>% tidyr::spread(role, count, fill=0) 
  inc_inf_by <-   generation_interval_matrix %>% group_by(infected_by, infectious_time)      %>% summarise(count=n()) %>% tidyr::spread(infected_by, count, fill=0)      
  inc_deg_rec <-  generation_interval_matrix %>% group_by(degree_recurrent, infectious_time) %>% summarise(count=n()) %>% tidyr::spread(degree_recurrent, count, fill=0) 
  
  ggplot(inc_role_raw %>% filter(!is.na(infectious_time)) %>% mutate(role=factor(role)), aes(x=infectious_time, y=count, color=role)) + geom_line() + theme_bw()
  
  outbreak_params <- list(timesteps=timesteps, infectious_period=infectious_period, initial_exposed=initial_exposed, transmission_prob=transmission_prob, 
                          directionality=directionality, incubation_period=incubation_period, rec_daily=rec_daily, tag=tag, seed=seed)
  
  return(list(generation_interval_matrix=generation_interval_matrix, 
              incidence_vector=incidence_vector, exposed_vector=exposed_vector, 
              inc_role=inc_role, inc_inf_by=inc_inf_by, inc_deg_rec=inc_deg_rec, outbreak_params=outbreak_params))
  
}

for(i in 1:100){
  print(i)
  p <- proc.time()
  graph <- MSM_like_graph_generator(n = lhs$num_nodes[i], alpha = lhs$dd_param[i], kmax = lhs$dd_upper[i], main_partners_prop = c('0'=lhs$p0[i], '1'=lhs$p1[i], '2'=lhs$p2[i], '3'=lhs$p3[i], '4'=lhs$p4[i]),
                                    roles_prop = c('1'=(1-lhs$vers[i])/2, '2'=(1-lhs$vers[i])/2, '3'=lhs$vers[i]), role_ass_param = lhs$ass_v_param[i],
                                    spatial_ass = lhs$spatial_kernel[i], degree_ass = lhs$assortativity_kernel[i], tag = lhs$tag[i], stop_prop = 0.02, folder="MSM_like")
  
  tot_inf <- 99
  seed <- 1
  
  while(tot_inf < 500){
    outbreak <- outbreak_simulator_fast5(timesteps = 730, infectious_period = lhs$infectious_time[i], initial_exposed = 1, transmission_prob = lhs$transmission_prob[i], directionality = lhs$directionality[i],
                                         incubation_period = lhs$incubation_period[i], tag = tag, rec_daily = lhs$rec_daily[i], g = graph, seed=seed)
    
    tot_inf <- sum(outbreak$incidence_vector)
    print(tot_inf)
    seed <- seed + 1
    if (seed > 10) tot_inf = 101 
  }
  
  output <- list(graph=graph, outbreak=outbreak)
  print((proc.time()-p)["elapsed"])
  saveRDS(output, file=paste0("simulation2/MSM_like/epidemic_", i, ".rds"))
}


fit_func2 <- function(simulation, model_name, cutoff, epi_phase=NA, dir="", iter=1000, chains=4, adapt_delta=0.9, tag=NA, keep_fit=F, seed){
  set.seed(5)
  model <- model_list[[model_name]]
  
  if(length(simulation)>1){
    infectious_time <- simulation$outbreak$outbreak_params$infectious_period
    incubation_period <- simulation$outbreak$outbreak_params$incubation_period
    incidence_stratified <- simulation$outbreak$incidence_vector
    weekly_inc <- rollsum(incidence_stratified, k=7, fill=NA)
    
    first_over10 <- which(weekly_inc>10)[1] 
    
    t0.1 <- first_over10[1]-30-1
    if(t0.1 < cutoff & t0.1>0){
      incidence_stratified <- c(rep(0,7),simulation$outbreak$incidence_vector[-c(1:(t0.1))])
      cutoff <- cutoff-t0.1+7
    }
  
    k <- 1
    S0 <- array(as.integer(simulation$graph$graph_params$n), dim = 1)
  }
  
  print(cutoff, quote=F)
  
  if(length(simulation)==1){
    #incidence_stratified <- list()
    first_nonzero <- which(dat[simulation] != 0)[1]
    incidence_stratified <- c(rep(0,7),dat[simulation][first_nonzero:nrow(dat[simulation]),])
    infectious_time <- 14
    incubation_period <- 7
    k <- 1
    S0 <- array(as.integer(MSM_pop[simulation]), dim=1)
    
    weekly_inc <- rollsum(incidence_stratified, k=7, fill=NA)
    
    first_over10 <- which(weekly_inc>10)[1] 
    
    t0.1 <- first_over10[1]-30-1
    if(t0.1 < cutoff & t0.1>0){
      incidence_stratified <- c(rep(0,7),simulation$outbreak$incidence_vector[-c(1:(t0.1))])
      cutoff <- cutoff-t0.1+7
    }
    
  }
  
  input_list <- list(N=cutoff, 
                     h=28, 
                     k=1, 
                     incidence_strat=matrix(incidence_stratified[1:cutoff], ncol=k),
                     S0=S0, 
                     model=model_name, 
                     max_lag=(incubation_period+infectious_time), L=(incubation_period+infectious_time), 
                     gen_weights=c(rep(0,incubation_period), rep(1/infectious_time, infectious_time)), 
                     epi_phase=epi_phase)
  
  if(k == 1) R0 = array(4, dim=1)
  if(k != 1) R0 = rep(4, k)
  
  init_raw <- list(R0 = R0)
  init <- replicate(chains, init_raw, simplify = FALSE)
  
  print(paste0("Model name = ", model_name, "     Epi phase = ", epi_phase, "     Tag = ", tag), quote=F)
  
  seed <- seed
  
  p <- proc.time()
  fit <- sampling(model, input_list, iter=iter, chains=chains, cores=4, seed=seed, init=init, control = list(adapt_delta = adapt_delta)) 
  
  max_Rhat <- max(summary(fit)$summary[,"Rhat"])
  divergences <- sapply(rstan::get_sampler_params(fit, inc_warmup = FALSE), function(x) sum(x[, "divergent__"]))
  total_divergences <- sum(divergences)
  
  print(paste0("Rhat = ", round(max_Rhat, 3), " Divergences = ", total_divergences), quote=F)
  
  summary <- rbind(summary(fit)$summary[1:6,], summary(fit)$summary %>% tail(1)) 
  
  time_taken <- round((proc.time()-p)["elapsed"],2)
  
  print(paste0("time taken = ", time_taken), quote=F)
  if(total_divergences>1) cat("\n\n\n\n\n\n\n")
  
  output2 <- list(infections=apply(rstan::extract(fit)$infections, c(1,2), sum), 
                  forecast=apply(rstan::extract(fit)$forecast, c(1,2), sum), 
                  epi_phase=epi_phase,#range_limits=range_limits,
                  input_list=input_list, tag=tag, model_name=model_name, #strat=strat, 
                  cutoff=cutoff, summary=summary,
                  incidence_full = matrix(incidence_stratified[1:(cutoff+input_list$h)], ncol=k), 
                  max_Rhat = max_Rhat, divergences = divergences, total_divergences = total_divergences, seed=seed, adapt_delta=adapt_delta, iter=iter)
  
  if(length(simulation)==1){
    output2$simulation <- NA
    output2$num_nodes <- 100000
    folder <- "cities"
  }
  
  if(length(simulation)>1){
    output2$simulation <- simulation
    output2$num_nodes <- 
    output2$param_list <- simulation$param_list
    folder <- "sim"
  }
  
  output2$time_taken <- time_taken
  output2$file_name <- paste0("fits/", folder,"/fit_", model_name,"_k=", k, "_epi_phase=", epi_phase,"_", tag, ".rds")
  
  plotting_func(output2)
  
  if(keep_fit==T) output2$fit <- fit
  
  saveRDS(output2, output2$file_name)
  
  return(output2)
}

#####
lrwP <- stan_model("functions/log_random_walkP.stan")
R_random5 <- stan_model("functions/R_random5.stan")
R_biased5 <- stan_model("functions/R_biased5.stan")
sigmoidP9 <- stan_model("functions/sigmoidP9.stan")
expI4 <- stan_model("functions/expI4.stan")
expI4.2 <- stan_model("functions/expI4.2.stan")

model_list <- list(lrwP=lrwP, expI4=expI4, sigmoidP9=sigmoidP9, R_random5=R_random5, R_biased5=R_biased5) #sigmoidP=sigmoidP, sigmoidPS=sigmoidPS, 

model_list <- list(expI4.2=expI4.2) #sigmoidP=sigmoidP, sigmoidPS=sigmoidPS, 

#####

plotting_func <- function(fit_output){
  input <- fit_output$input_list
  #fit <- fit_output$fit
  incidence_full <- fit_output$incidence_full
  
  df1 <- cbind(data.frame(x=1:input$N, matrixStats::colQuantiles(apply(fit_output$infections,c(1,2),sum), probs=c(0.025, 0.5, 0.975))[1:input$N,], incidence=c(rowSums(as.matrix(input$incidence_strat))))) %>% magrittr::set_colnames(c("x", "lower", "median", "upper", "incidence"))
  df2 <- cbind(data.frame(x=(input$N+1):(input$N+input$h), matrixStats::colQuantiles(apply(fit_output$forecast,  c(1,2),sum), probs=c(0.025, 0.5, 0.975)), incidence=NA)) %>% magrittr::set_colnames(c("x", "lower", "median", "upper", "incidence"))
  
  df <- rbind(df1, df2) %>% mutate(incidence=rowSums(incidence_full))
  
  p <- ggplot(df, aes(x=x)) + geom_point(aes(y=incidence)) + geom_line(aes(y=median)) + geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) + theme_bw()
  
  inf <- fit_output$infections
  forecast <- fit_output$forecast
  
  q_array <- apply(array(inf, dim=c(dim(inf), 1)), c(2, 3), quantile, probs = c(0.025, 0.5, 0.975))
  
  df2 <- melt(q_array, varnames = c("quantile", "time", "stratum")) %>%
    tidyr::pivot_wider(names_from = quantile, values_from = value) %>%
    rename(lower = `2.5%`, median = `50%`, upper = `97.5%`)
  n_time <- input$N + input$h
  n_strata <- dim(inf)[3]
  
  inc <- rbind(incidence_full)
  df2$incidence <- as.vector(inc)
  df2$stratum <- factor(df2$stratum)
  
  q <- ggplot(df2, aes(x = time, fill=stratum, color=stratum)) +
    geom_point(aes(y = incidence)) +
    geom_line(aes(y = median)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    facet_wrap(~ stratum, scales="free_y") +
    theme_bw()
  
  #print(p)
  #print(q)
  
  r <- ggpubr::ggarrange(p, q, nrow=2)
  print(r)
  
}

sim_files <- c(list.files("simulation2/homo/", full.names=T))
nums <- as.numeric(sub(".*epidemic_(\\d+)\\.rds$", "\\1", sim_files))
fit_files <- sim_files[order(nums)]

sim_files2 <- c(list.files("simulation2/hetero", full.names=T), list.files("simulation2/homo", full.names=T))
fit_files2 <- sim_files2[order(file.info(sim_files2)$mtime)]

for(j in fit_files[1:10]){
  simulation <- readRDS(j)
  plot(simulation$outbreak$incidence_vector)
}

for(j in fit_files){
  print(j)
  simulation <- readRDS(j)
  
  inc <- simulation$outbreak$incidence_vector
  cum_inc <- cumsum(inc)
  tot_inc <- sum(inc)
  peak_inc <- max(inc)
  peak_inc_t <- which(inc==max(inc))
  first_lower_3 <- which(inc<(peak_inc*0.03) & (1:length(inc)) > peak_inc_t)[1]
  
  lower <- which(cum_inc > tot_inc * 0.05)[1]
  if( grepl("hetero", j)) upper <- first_lower_3 
  if(!grepl("hetero", j)) upper <- which(cum_inc > tot_inc * 0.95)[1]-28 
  
  times <- round(seq(lower, upper, length.out=10))
  
  tag <- simulation$graph$tag
  print(tag)
  dir <- sub("^[^/]+/([^/]+)/.*", "\\1", j)
  
  plot(simulation$outbreak$incidence_vector[1:(times[10]+28)])
  
  for(i in times[3]){
    for(q in "lrwP"){#names(model_list)){
      print(i)
      a7 <- fit_func2(simulation, model_name=q, cutoff=i, epi_phase=which(times==i), dir=dir, iter=1000, chains=4, tag=paste0(tag,"T"), keep_fit=F, seed=31)
    }
  }
}

adjust_weekday <- function(df, date_col = "date", count_col = "confirm") {
  df$dow <- factor(weekdays(df[[date_col]]),
                   levels = c("Monday","Tuesday","Wednesday","Thursday",
                              "Friday","Saturday","Sunday"))
  
  m    <- glm(df[[count_col]] ~ dow, family = quasipoisson(link = "log"), data = df)
  co   <- coef(m)
  levs <- levels(df$dow)
  
  eff <- setNames(numeric(length(levs)), levs)
  eff["Monday"] <- exp(co["(Intercept)"])
  for (lv in levs[-1]) eff[lv] <- exp(co["(Intercept)"] + co[paste0("dow", lv)])
  
  rel <- eff / mean(eff)
  
  print(round(data.frame(multiplier = rel), 3))
  
  df[[paste0(count_col, "_adj")]] <- round(df[[count_col]] / rel[df$dow])
  df
}

NYC_data_raw <- read.csv("data/NYC_cases.csv") %>% mutate(date=as.Date(diagnosis_date, format="%d/%m/%Y")) %>% select(date, count) %>% rename(confirm=count)
NYC_data <- adjust_weekday(NYC_data_raw) %>% rename(NYC=confirm_adj, NYC_raw=confirm) #%>% select(date, NYC)

SF_data_raw <- read.csv("data/SF_cases.csv") %>% mutate(date=as.Date(episode_date, format="%d/%m/%Y")) %>% mutate(confirm=new_cases) %>% select(date, confirm)
SF_data <- adjust_weekday(SF_data_raw) %>% rename(SF=confirm_adj, SF_raw=confirm) #%>% select(date, SF)

adj_plot_df <- full_join(SF_data, NYC_data, by="date") %>% select(date, SF_raw, SF, NYC_raw, NYC) %>%
  tidyr::gather(measure, value, 2:5) %>%
  mutate(place=ifelse(grepl("SF", measure), "SF", "NYC"),
         measure2=ifelse(grepl("raw", measure), "Unadjusted", "Adjusted")) %>%
  filter(date < as.Date("2023-01-01")) %>%
  mutate(measure2 = factor(measure2, levels = (unique(measure2)))) 

adj_plot <- ggplot(adj_plot_df, aes(x = date, y = value)) +
  geom_line(aes(color = measure2)) +
  facet_grid(measure2 ~ place, scales = "free",
             labeller = labeller(place = c(
               "NYC" = "New York City",
               "SF"  = "San Francisco",
               "LON" = "London",
               "BCN" = "Barcelona",
               "MAD" = "Madrid"
             ))) +
  scale_x_date(date_labels = "%b %Y",
               date_breaks = "2 months") +
  labs(color = NULL, y = "Cases", x = "") +
  theme_bw() +
  theme(panel.grid       = element_blank(),
        strip.background = element_rect(fill = "white", colour = NA),
        strip.text       = element_text(face = "bold", size = 10),
        axis.text.x      = element_text(angle = 45, hjust = 1, size = 8))

ggsave("figures/NYC_SF_adjustment.png", adj_plot, width=6, height=6)

NYC_data <- df %>% dplyr::select(date, confirm_adj) %>% rename(NYC=confirm_adj)
#SF_data <- read.csv("data/SF_cases.csv") %>% mutate(date=as.Date(episode_date, format="%d/%m/%Y"), SF=new_cases) %>% select(date, SF)  %>% mutate(day_type = ifelse(weekdays(date) %in% c("Saturday", "Sunday"), "weekend", "weekday"))
SPN_data <- readxl::read_xlsx("data/mpox_BCN_MAD.xlsx", sheet="Incidence") %>% mutate(date=as.Date(Date), MAD=Madrid, BCN=Barcelona) %>% select(date, MAD, BCN)
skeleton <- data.frame(date=seq(min(c(NYC_data$date, SF_data$date, SPN_data$date)), max(c(NYC_data$date, SF_data$date, SPN_data$date)), by="day"))

dat <- full_join(skeleton, NYC_data, by="date") %>%
  full_join(., SF_data, by="date") %>%
  full_join(., SPN_data, by="date") %>%
  mutate(across(c(NYC, SF, MAD, BCN), ~ replace_na(., 0))) %>%
  filter(row_number() < 365) %>%
  select(-dow.x) %>% select(-dow.y)

ggplot(dat[1:200,] %>% tidyr::gather("city", "value", 2:5), aes(x=date, y=value)) +
  geom_line(aes(color=city)) +
  theme_bw() +
  theme(panel.grid = element_blank())

pop <- c(BCN=1.73e6,  #https://portaldades.ajuntament.barcelona.cat/en/statistics/yzlntdm2fs
         MAD=3.4e6,   #https://www.citypopulation.de/en/spain/madrid/madrid/28079__madrid/
         NYC=7.94e6,  #https://worldpopulationreview.com/us-cities
         SF=0.768e6,
         LON=9.8e6)   #https://worldpopulationreview.com/cities/united-kingdom/london

#prop_MSM <- 0.038

#MSM_pop <- pop*prop_MSM
#
#MSM_pop["LON"] <- 0.038*pop["LON"]*0.5*0.5        # Inequalities in Sexual Health. Update on HIV and STIs in men who have sex with men in London 
#MSM_pop["NYC"] <- round(397399*0.5)               # Estimating Population Sizes of Men Who Have Sex with Men in the United States # chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://depts.washington.edu/hivtcg/presentations/uploads/35/estimating_population_sizes_of_men_who_have_sex_with_men_in_the_united_states.pdf?utm_source=chatgpt.com
#MSM_pop["SF"] <- round(145972*0.5)                # Estimating Population Sizes of Men Who Have Sex with Men in the United States # chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://depts.washington.edu/hivtcg/presentations/uploads/35/estimating_population_sizes_of_men_who_have_sex_with_men_in_the_united_states.pdf?utm_source=chatgpt.com
#MSM_pop["BCN"] <- round(pop["BCN"]*0.5*0.5*0.053) # https://pmc.ncbi.nlm.nih.gov/articles/PMC4594901/#:~:text=A%20total%20of%201560%20cases,was%20only%2011.9/100%2C000%20inhabitants.
#MSM_pop["MAD"] <- round(pop["MAD"]*0.5*0.5*0.053) # https://pmc.ncbi.nlm.nih.gov/articles/PMC4594901/#:~:text=A%20total%20of%201560%20cases,was%20only%2011.9/100%2C000%20inhabitants.

MSM_pop <- c(LON=98330, # https://journals.sagepub.com/doi/pdf/10.1258/ijsa.2010.010181 # Men who have sex with men: estimating the size of at-risk populations in London primary care trust
             SF=70000, # https://link.springer.com/article/10.1007/s10461-018-2321-0
             NYC=235000) # http://nature.com/articles/s41591-025-03526-9

MSM_pop["BCN"] <- round(MSM_pop["LON"]/pop["LON"]*pop["BCN"],-3)
MSM_pop["MAD"] <- round(MSM_pop["LON"]/pop["LON"]*pop["MAD"],-3)

for(j in c("BCN", "MAD", "NYC", "SF")){
  simulation <- j
  
  cum_inc_raw <- cumsum(dat[,j])
  first_nonzero <- which(cum_inc_raw != 0)[1]
  cum_inc <- c(rep(0,7), cum_inc_raw[first_nonzero:length(cum_inc_raw)])
  
  tot_inc <- max(cum_inc)
  lower <- which(cum_inc > tot_inc * 0.05)[1]
  upper <- which(cum_inc > tot_inc * 0.95)[1]-28
  times <- round(seq(lower, upper, length.out=10))
  
  tag <- paste0(j,3)
  
  for(i in times){
    for(p in c("lrwP", "R_random5", "R_biased5", "sigmoidP9", "expI4")){
      max_Rhat <- 1.06
      total_divergences <- 1
      seed <-30
      
      while(is.nan(max_Rhat) | max_Rhat > 1.05 | total_divergences > 0){
        a0 <- fit_func2(simulation, seed=seed, model_name=p, cutoff=i, epi_phase=which(times==i), dir="CIT", tag=tag, adapt_delta=0.95, iter=3000)
        max_Rhat <- a0$max_Rhat
        total_divergences <- a0$total_divergences
        print(c(max_Rhat, total_divergences))
        seed <- seed + 1
        print(seed)
      }
      
    }
  }
}

