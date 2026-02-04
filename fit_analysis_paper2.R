library(dplyr)
library(scoringutils)
library(ggplot2)
library(boot)
library(scales)
library(tidyr)
library(glmnet)
library(purrr)

metric_levels <- c("crps", "overprediction", "dispersion", "underprediction")
model_levels <- c("lrwP", "R_random5", "R_biased5", "expdelayP", "incdelayP", "powdelayP", "sigmoidP")
place_levels <- c("BCN", "LON", "MAD", "NYC", "SF", "HOMO", "HET", "SIM")

id_cols <- c("model", "cutoff", "epi_phase", "place", "sr")
full_names <- c(lrwP="LRW", R_random5="RRAND", R_biased5="RBIAS", expdelay="EXPT", incdelayP="EXPI", powdelay="POWI", sigmoidP="SIGT")

files_div <- readRDS("fits/outputs/Rhat_sim.rds") %>% filter(divergences>=1) %>% mutate(scenario=as.character(scenario))

skeleton_plot_sim      <- readRDS("fits/outputs/skeleton_plot_sim_all.rds")      %>% mutate(model=ifelse(model=="sigmoidPS", "sigmoidP", model))
skeleton_plot_sim_clean <- skeleton_plot_sim %>% anti_join(files_div %>% select(model, scenario, sr, epi_phase) %>% distinct(), by = c("model", "place"="scenario", "sr", "epi_phase"))

skeleton_last_week_sim <- readRDS("fits/outputs/skeleton_last_week_sim_all.rds") %>% mutate(model=ifelse(model=="sigmoidPS", "sigmoidP", model))
skeleton_last_week_sim_clean <- skeleton_last_week_sim %>% anti_join(files_div %>% select(model, scenario, sr, epi_phase) %>% distinct(), by = c("model", "place"="scenario", "sr", "epi_phase"))

skeleton_plot_cit <- readRDS("fits/outputs/skeleton_plot_cit.rds")
skeleton_last_week_cit <- readRDS("fits/outputs/skeleton_last_week_cit.rds")

skeleton_plot <- rbind(skeleton_plot_sim_clean, skeleton_plot_cit)
skeleton_last_week_raw <- rbind(skeleton_last_week_sim_clean, skeleton_last_week_cit)

min_n <- skeleton_last_week_raw %>% group_by(epi_phase, place, sr, model) %>% summarise(n = n(), .groups = "drop") %>% summarise(min_n = min(n)) %>% pull(min_n)
skeleton_last_week <- skeleton_last_week_raw %>% group_by(epi_phase, place, sr, model) %>% arrange(sample, .by_group = TRUE) %>% slice(round(seq(n()/min_n, n(), length.out = min_n))) %>% ungroup() %>% arrange(sr, place, model, epi_phase, sample)

skeleton_last_week %>% group_by(epi_phase, place, sr, model) %>% summarise(n = n(), .groups = "drop")

fc <- as_forecast_sample(data = skeleton_last_week, observed = "true_value", sample = "sample", forecast_unit = id_cols) %>% mutate(grouping = case_when(sr == "fit_sim3" ~ "SIM", sr == "fit_homo" ~ "HOMO", sr == "fit_het"  ~ "HET", sr == "CIT"      ~ place, TRUE ~ NA_character_))
fc_log <- transform_forecasts(fc, offset = 1, append = FALSE, label = "log") 

metrics     <- get_metrics(fc,     select = c("crps", "overprediction", "underprediction", "dispersion"))
metrics_log <- get_metrics(fc_log, select = c("crps", "overprediction", "underprediction", "dispersion"))

scoring_lin0 <- as.data.frame(score(fc,     metrics = metrics)) %>% mutate(measure="Linear CRPS")
scoring_log0 <- as.data.frame(score(fc_log, metrics = metrics)) %>% mutate(measure="Log CRPS")

scoring_lin_LON <- read.csv("fits/outputs/scoring_lin_LON.csv") %>% filter(model != "sigmoidPS") %>% select(-X) %>% mutate(grouping="LON")
scoring_log_LON <- read.csv("fits/outputs/scoring_log_LON.csv") %>% filter(model != "sigmoidPS") %>% select(-X) %>% mutate(grouping="LON")

scoring_lin <- rbind(scoring_lin0, scoring_lin_LON)
scoring_log <- rbind(scoring_log0, scoring_log_LON)

scoring_results_raw0 <- rbind(scoring_lin, scoring_log) %>% 
  pivot_longer(cols = c(crps, overprediction, underprediction, dispersion), names_to = "metric", values_to = "value") %>%
  mutate(model=factor(model, levels=model_levels),
         grouping = factor(grouping, levels=c(place_levels)))

scoring_results_raw <- scoring_results_raw0 %>%
  group_by(model, sr, epi_phase, grouping, metric, measure) %>%
  summarise(value=mean(value)) %>%
  mutate(metric = factor(metric, levels=metric_levels)) %>% 
  mutate(model2=factor(full_names[model], levels=full_names))

scoring_results_agg <- scoring_results_raw0 %>%
  group_by(model, grouping, measure, metric) %>%
  summarise(value=mean(value)) 

scoring_results_agg %>% filter(metric=="crps") %>% mutate(value=round(value, ifelse(measure=="Linear CRPS", 0, 2))) %>% tidyr::spread(grouping, value) %>% select(!metric) %>% arrange(measure, model) 

R_case <- function(simul_vec){
  for(k in simul_vec){
    print(k)
    b <- readRDS(paste0("fits/sim/fit_sim3/fit_sigmoidP_k=1_strat=annual_likelihood=1_epi_phase=7_",k,".rds"))
    c <- rstan::extract(b$fit)
    
    gT_dist <- b$input_list$gen_weights[b$input_list$gen_weights!=0]
    tmax <- length(gT_dist)
    
    mean_R_df <- b$simulation$generation_interval_matrix %>%
      group_by(person, t_I_start) %>%
      left_join(b$simulation$generation_interval_matrix %>% group_by(infector) %>% summarise(onward_infections=n()), by = c("person" = "infector")) %>%
      select(t_I_start, person, onward_infections) %>%
      mutate(onward_infections = ifelse(is.na(onward_infections), 0, onward_infections)) %>%
      group_by(t_I_start) %>%
      summarise(mean_R=mean(onward_infections)) %>%
      mutate(rollmean_R = zoo::rollmean(mean_R, 7, fill=NA)) %>%
      filter(t_I_start < b$input_list$N+b$input_list$h-tmax)
    
    ggplot(mean_R_df, aes(x=t_I_start, y=rollmean_R)) +
      geom_point() +
      theme_minimal()
    
    j_vec <- c(0,seq_len(b$input_list$N+b$input_list$h))  # for example
    
    R_mat <- sapply(j_vec, function(jk) c$R_floor + (c$R0 - c$R_floor) * plogis(-c$delta * (jk-as.vector(c$t_change))))
    #R_vec <- R_mat[1,]
    
    imm_mat <- (1-cbind(0,c$cumulative_infections[,,1]/b$param_list$num_nodes))
    #imm_vec <- imm_mat[1,]
    
    
    #R_case <- vector(length=length(R_vec))
    #for(i in 1:length(R_case)) R_case[i] <- sum(R_vec[i:(i+tmax-1)]*gT_dist*imm_vec[i:(i+tmax-1)])
    
    R_case_mat <- matrix(nrow=nrow(R_mat), ncol=ncol(R_mat)-tmax)
    for(i in 1:nrow(R_case_mat)){
      R_vec <- R_mat[i,]
      imm_vec <- imm_mat[i,]
      for(j in 1:ncol(R_case_mat)) R_case_mat[i,j] <- sum(R_vec[j:(j+tmax-1)]*gT_dist*imm_vec[j:(j+tmax-1)])
    }
    
    col_ci <- apply(R_case_mat, 2, quantile, probs = c(0.025, 0.5, 0.975))[,mean_R_df$t_I_start+1]
    
    if(k==simul_vec[1]) R_plotting <- mean_R_df %>% mutate(R_lower=col_ci[1,], R_median=col_ci[2,], R_upper=col_ci[3,]) %>% mutate(simulation=k)
    else R_plotting <- rbind(R_plotting, mean_R_df %>% mutate(R_lower=col_ci[1,], R_median=col_ci[2,], R_upper=col_ci[3,]) %>% mutate(simulation=k))
    
  }
  
  p <- ggplot(R_plotting, aes(x=t_I_start)) +
    geom_line(aes(y=R_median)) +
    geom_ribbon(aes(ymin=R_lower, ymax=R_upper), alpha=0.2) +
    geom_point(aes(y=rollmean_R)) +
    theme_minimal() +
    facet_wrap(~simulation, scales="free") +
    theme(panel.grid = element_blank()) +
    labs(x="Time (days)", y=("Reproduction Number"))
  
  print(p)
  
  return(p)
}

fig <- R_case(simul_vec=seq(11,19,1))
ggsave("figures/R_case_comp.png", fig, width=9, height=9)

ts_comb <- ggplot(scoring_results_raw %>% filter(metric != "crps", grouping %in% c("SIM")) %>% mutate(epi_phase=as.numeric(as.character(epi_phase)))) +
  geom_bar(aes(x = epi_phase, y = value, fill = metric),
           position = "stack",
           stat = "identity") +
  facet_grid(measure ~ model2, switch = "x", scale="free_y") +
  theme_minimal() + 
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm"),
        panel.grid=element_blank(), 
        legend.position="none") +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  labs(x="", title="CRPS measured on the linear forecasting scale")

ggsave("figures/ts_comb.png", ts_comb, width=12, height=7.5)

cit_epi_phase <- ggplot(scoring_results_raw %>% filter(measure=="Linear CRPS", metric != "crps", sr %in% c("CIT")) %>% mutate(epi_phase=as.numeric(as.character(epi_phase)))) +
  geom_bar(aes(x = epi_phase, y = value, fill = metric),
           position = "stack",
           stat = "identity", color=NA) +
  facet_grid(grouping~model2, switch = "x", scales="free_y") +
  theme_minimal() + 
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm"),
        panel.grid=element_blank(),
        legend.position="none") +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  labs(x="", title="CRPS measured on the linear forecasting scale")

cit_epi_phase

ggsave("figures/crps_epi_phase_cit.png", cit_epi_phase, width=10, height=10)

comparison_func_real <- function(skeleton_plot, model_list, epi_phase_list, place_list, cutoff_graph=F, sr_list="SIM"){
  if(cutoff_graph==T) skeleton_plot <- skeleton_plot %>% filter(cutoff<=day+10)
  
  skeleton_plot_fil <- skeleton_plot %>% filter(model %in% model_list, epi_phase %in% epi_phase_list, place %in% place_list, sr %in% sr_list) %>%
    mutate(model2=factor(full_names[model], levels=full_names))
  
  p1 <- ggplot(skeleton_plot_fil %>% filter(true_value>0 | day > cutoff), aes(x=day, fill=factor(epi_phase))) +
    geom_point(aes(y=true_value)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
    geom_line(aes(y=median)) +
    geom_vline(aes(xintercept=cutoff), linetype="dashed") +
    theme_minimal() +
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = NA, color = "white"),
          panel.spacing = unit(-.01,"cm"),
          panel.grid=element_blank(),
          legend.position="none") +
    #scale_y_log10() +
    facet_grid(model~place, scales="free") +
    theme(panel.grid=element_blank()) 
  
  p1
  
  p2 <- ggplot(skeleton_plot_fil %>% filter(true_value>=1, median>=1 | day > cutoff, day > 20), aes(x=day, fill=factor(epi_phase))) +
    geom_point(aes(y=true_value), size=0.5) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
    geom_line(aes(y=median)) +
    geom_vline(aes(xintercept=cutoff), linetype="dashed") +
    theme_minimal() + 
    #scale_fill_brewer(palette = "Set3") +
    scale_fill_manual(values = c(
      "#007BFF",  # bright blue
      "#9B00FF",  # bright purple
      "#FF1493"   # bright pink
    )) +
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = NA, color = NA),
          panel.spacing = unit(0.2,"cm"),
          panel.grid=element_blank(),
          legend.position="none",
          axis.line = element_line(),
          axis.ticks = element_line()) +
    scale_y_log10() +
    facet_grid(model2~place, scales="free") +
    theme(panel.grid=element_blank()) +
    labs(x="Day", y="Incidence") 
  
  p2
  #q <- scoring_results_real_raw %>% filter(model %in% model_list) 
  
  print(ggpubr::ggarrange(p1, p2, nrow=1))
  print(p1)
  print(p2)
  return(list(p1, p2))
}

place_comp <- comparison_func_real(skeleton_plot, model_list=c("lrwP", "R_biased5", "R_random5", "sigmoidP"), epi_phase_list=c(2, 6, 10), place_list=c("BCN", "SF", "MAD", "NYC"), cutoff_graph=F, sr_list="CIT")
#place_comp <- comparison_func_real(skeleton_plot, model_list=c("sigmoidP"), epi_phase_list=c(6), place_list=c(4), cutoff_graph=F, sr_list="CIT")

ggsave("figures/place_comp.png", place_comp[[2]], width=15, height=9)

g <- rbind(scoring_lin, scoring_log) %>% filter(sr=="CIT") %>% group_by(epi_phase, place, measure) %>% 
  mutate(rank=rank(crps)) %>% group_by(model, measure) %>% 
  summarise(av_rank=mean(rank), low_rank = quantile(rank, 0.25), high_rank=quantile(rank, 0.75)) %>%
  rename(score_type=measure) %>% arrange(av_rank)

offset <- 0.1

cit_lol <- ggplot(g, aes(color = score_type)) +
  geom_errorbarh(aes(xmin = low_rank, xmax = high_rank,
                     y = as.numeric(reorder(model, av_rank)) + ifelse(score_type == "Linear CRPS", -offset, +offset)),
                 height = 0.15, alpha  = 0.6, linewidth = 1) +
  geom_point(aes(x = av_rank, y = as.numeric(reorder(model, av_rank)) +
                   ifelse(score_type == "Linear CRPS", -offset, +offset)), size = 3) +
  scale_y_continuous(breaks = 1:length( levels(reorder(g$model, g$av_rank))), labels =  levels(reorder(g$model, g$av_rank))) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "Average rank across cities and epidemic phases", y = "Model")

cit_lol

ggsave("figures/prop_win_cit.png", cit_lol, width=5, height=5)









skeleton_plot_place <- readRDS("fits/outputs/skeleton_plot_cit.rds") %>% mutate(sr="real", epi_phase=as.numeric(as.character(epi_phase)), cutoff2=factor(cutoff)) %>% rename(place=scenario) %>% filter(!place %in% c("BG", "LON", "LON_SIM"))
skeleton_last_week_place <- readRDS("fits/outputs/skeleton_last_week_cit.rds") %>% mutate(sr="real", epi_phase=as.numeric(as.character(epi_phase)), cutoff2=factor(cutoff)) %>% rename(place=scenario) %>% filter(!place %in% c("BG", "LON", "LON_SIM"))

skeleton_plot_LON <- read.csv("fits/outputs/skeleton_plot_LON_censored.csv") 
skeleton_last_week_LON <- read.csv("fits/outputs/skeleton_last_week_LON_censored.csv") 

skeleton_plot_sim <- readRDS("fits/outputs/skeleton_plot_sim_all15.rds") #%>% mutate(place=as.character(scenario), sr="SIM", epi_phase=as.numeric(as.character(epi_phase))) %>% select(day, true_value, cutoff, model, epi_phase, total_infection, lower, median, upper, cutoff2, place, sr) 
skeleton_last_week_sim <- readRDS("fits/outputs/skeleton_last_week_sim_all15.rds") #%>% mutate(place=as.character(scenario), sr="SIM", epi_phase=as.numeric(as.character(epi_phase))) %>% select(sample, model, cutoff, epi_phase, total_infection, prediction, true_value, place, sr) 

skeleton_plot_homo <- readRDS("fits/outputs/skeleton_plot_homo.rds") 
skeleton_last_week_homo <- readRDS("fits/outputs/skeleton_last_week_homo.rds") 

skeleton_plot_het <- readRDS("fits/outputs/skeleton_plot_het2.rds") 
skeleton_last_week_het <- readRDS("fits/outputs/skeleton_last_week_het2.rds") 

skeleton_plot <- bind_rows(skeleton_plot_place, skeleton_plot_sim, skeleton_plot_homo, skeleton_plot_het) %>% filter(!model %in% c("lindelay", "incdelay", "gompertz"))
skeleton_last_week <- bind_rows(skeleton_last_week_place, skeleton_last_week_sim, skeleton_last_week_homo, skeleton_last_week_het) %>% filter(!model %in% c("lindelay", "incdelay", "gompertz"))

comparison_func_real <- function(skeleton_plot, model_list, epi_phase_list, place_list, cutoff_graph=F, sr_list="SIM"){
  if(cutoff_graph==T) skeleton_plot <- skeleton_plot %>% filter(cutoff<=day+10)
  
  skeleton_plot_fil <- skeleton_plot %>% filter(model %in% model_list, epi_phase %in% epi_phase_list, place %in% place_list, sr %in% sr_list) %>%
    mutate(model2=factor(full_names[model], levels=full_names))
  
  p1 <- ggplot(skeleton_plot_fil %>% filter(true_value>0 | day > cutoff), aes(x=day, fill=factor(epi_phase))) +
    geom_point(aes(y=true_value)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
    geom_line(aes(y=median)) +
    geom_vline(aes(xintercept=cutoff), linetype="dashed") +
    theme_minimal() +
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = NA, color = "white"),
          panel.spacing = unit(-.01,"cm"),
          panel.grid=element_blank(),
          legend.position="none") +
    #scale_y_log10() +
    facet_grid(model~place, scales="free") +
    theme(panel.grid=element_blank()) 
  
  p2 <- ggplot(skeleton_plot_fil %>% filter(true_value>=1, median>=1 | day > cutoff, day > 20), aes(x=day, fill=factor(epi_phase))) +
    geom_point(aes(y=true_value), size=0.5) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
    geom_line(aes(y=median)) +
    geom_vline(aes(xintercept=cutoff), linetype="dashed") +
    theme_minimal() + 
    #scale_fill_brewer(palette = "Set3") +
    scale_fill_manual(values = c(
      "#007BFF",  # bright blue
      "#9B00FF",  # bright purple
      "#FF1493"   # bright pink
    )) +
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = NA, color = NA),
          panel.spacing = unit(0.2,"cm"),
          panel.grid=element_blank(),
          legend.position="none",
          axis.line = element_line(),
          axis.ticks = element_line()) +
    scale_y_log10() +
    facet_grid(model2~place, scales="free") +
    theme(panel.grid=element_blank()) +
    labs(x="Day", y="Incidence") 
  
  p2
  #q <- scoring_results_real_raw %>% filter(model %in% model_list) 
  
  print(ggpubr::ggarrange(p1, p2, nrow=1))
  print(p1)
  print(p2)
  return(list(p1, p2))
}

place_comp <- comparison_func_real(skeleton_plot, model_list=c("rw", "R_biased", "R_random", "sigmoid"), epi_phase_list=c(2, 6, 10), place_list=c("BCN", "SF", "MAD", "NYC"), cutoff_graph=F, sr_list="real")
place_comp2 <- comparison_func_real(skeleton_plot %>% filter(place %in% c(1, 4, 7, 11)) %>% mutate(place=as.numeric(factor(place))), model_list=c("rw", "R_biased", "R_random", "sigmoid"), epi_phase_list=c(2, 6, 10), place_list=c(1, 2,3,4), cutoff_graph=F, sr_list="SIM")

place_comp2

ggsave("figures/place_comp_lin.png", place_comp[[1]], width=15, height=9)
ggsave("figures/place_comp.png", place_comp[[2]], width=15, height=9)
ggsave("figures/place_comp2.png", place_comp2[[2]], width=12, height=6)

fc <- as_forecast_sample(data = dplyr::distinct(skeleton_last_week %>% rename(predicted = prediction)), observed = "true_value", sample = "sample", forecast_unit = id_cols) %>% mutate(grouping=ifelse(sr=="real", place, sr))
fc_log <- transform_forecasts(fc, offset = 1, append = FALSE, label = "log") 

metrics <- get_metrics(fc, select = c("crps", "overprediction", "underprediction", "dispersion"))
metrics_log <- get_metrics(fc_log, select = c("crps", "overprediction", "underprediction", "dispersion"))

scoring_lin_LON <- read.csv("fits/outputs/scoring_lin_LON.csv")[,-c(1)]
scoring_log_LON <- read.csv("fits/outputs/scoring_log_LON.csv")[,-c(1)]

scoring_lin <- rbind(as.data.frame(score(fc, metrics = metrics)) %>% mutate(measure="Linear CRPS"), scoring_lin_LON)
scoring_log <- rbind(as.data.frame(score(fc_log, metrics = metrics_log)) %>% mutate(measure="Log CRPS"), scoring_log_LON)

f_prop <- rbind(scoring_lin, scoring_log) %>% filter(sr=="SIM") %>% 
  group_by(epi_phase, place, measure) %>% mutate(rank=rank(crps)) %>% 
  group_by(model, rank, measure, epi_phase) %>% filter(rank==1) %>% 
  summarise(count=n()) %>% 
  arrange(measure, count) %>% 
  group_by(measure, epi_phase) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()

library(ggridges)  

ridge_plot <- ggplot(f_prop, aes(x = epi_phase, y = model, height = prop, group = model, fill = measure)) +
  geom_ridgeline(scale = 1, alpha = 0.7) +
  facet_grid(~ measure) +
  theme_minimal() +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm"),
        panel.grid=element_blank(),
        legend.position="none") +
  labs(x = "Epidemic phase", y = "Model") +
  scale_x_continuous(breaks = scales::pretty_breaks())

ggsave("figures/prop_win_sim.png", ridge_plot, width=5, height=5)

g <- rbind(scoring_lin, scoring_log) %>% filter(sr=="real") %>% group_by(epi_phase, place, measure) %>% 
  mutate(rank=rank(crps)) %>% group_by(model, measure) %>% 
  summarise(av_rank=mean(rank), low_rank = quantile(rank, 0.25), high_rank=quantile(rank, 0.75)) %>%
  rename(score_type=measure) %>% arrange(av_rank)

offset <- 0.1

cit_lol <- ggplot(g, aes(color = score_type)) +
  geom_errorbarh(aes(xmin = low_rank, xmax = high_rank,
                     y = as.numeric(reorder(model, av_rank)) + ifelse(score_type == "Linear CRPS", -offset, +offset)),
                 height = 0.15, alpha  = 0.6, linewidth = 1) +
  geom_point(aes(x = av_rank, y = as.numeric(reorder(model, av_rank)) +
                   ifelse(score_type == "Linear CRPS", -offset, +offset)), size = 3) +
  scale_y_continuous(breaks = 1:length( levels(reorder(g$model, g$av_rank))), labels =  levels(reorder(g$model, g$av_rank))) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "Average rank across cities and epidemic phases", y = "Model")

cit_lol

ggsave("figures/prop_win_cit.png", cit_lol, width=5, height=5)

scoring_results_raw0 <- rbind(scoring_lin, scoring_log) %>% 
  #mutate(model_type = model_type_dict[model]) %>%
  pivot_longer(cols = c(crps, overprediction, underprediction, dispersion), names_to = "metric", values_to = "value") %>%
  mutate(model=factor(model, levels=model_levels),
         grouping = factor(grouping, levels=c(place_levels)))

full_names <- c(rw="Random walk", lrw="Log random walk", R_random3="R random walk", R_biased3="R biased random", expdelay="Exponential time", incdelay3="Exponential attack rate", powdelay="Power law attack rate", sigmoid="Sigmoid")

scoring_results_raw <- scoring_results_raw0 %>%
  group_by(model, sr, epi_phase, grouping, metric, measure) %>%
  summarise(value=mean(value)) %>%
  mutate(metric = factor(metric, levels=metric_levels)) %>% #%>%
  #mutate(value = ifelse(metric=="underprediction", -value, value))
  mutate(model2=factor(full_names[model], levels=full_names))

ts_lin <- ggplot(scoring_results_raw %>% filter(measure=="Linear CRPS", metric != "crps", grouping %in% c("SIM"))) +
  geom_bar(aes(x = epi_phase, y = value, fill = metric),
           position = "stack",
           stat = "identity") +
  facet_grid(~ model, switch = "x") +
  theme_minimal() + 
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm"),
        panel.grid=element_blank(), 
        legend.position="none") +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  labs(x="", title="CRPS measured on the linear forecasting scale")

metric_names <- c(overprediction="Overprediction", underprediction="Underprediction", dispersion="Dispersion")

ts_log <- ggplot(scoring_results_raw %>% mutate(Metric=metric_names[as.character(metric)]) %>% filter(measure=="Log CRPS", metric != "crps", grouping %in% c("SIM"))) +
  geom_bar(aes(x = epi_phase, y = value, fill = Metric),
           position = "stack",
           stat = "identity") +
  facet_grid(~ model2, switch = "x") +
  theme_minimal() + 
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        panel.spacing = unit(0.2,"cm"),
        panel.grid=element_blank(),
        legend.position="right",
        axis.line = element_line(),
        axis.ticks = element_line()) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  labs(x="Epidemic phase", y="CRPS", fill=NULL)

ts_log

ts_comb <- ggpubr::ggarrange(ts_lin, ts_log, ncol = 1, common.legend = TRUE, legend = "none", align = "h")
ggsave("figures/ts_log.png", ts_log, width=13, height=5)
ggsave("figures/ts_comb.png", ts_comb, width=12, height=7.5)

scoring_results_agg <- scoring_results_raw0 %>%
  group_by(model, grouping, measure, metric) %>%
  summarise(value=mean(value)) 

scoring_results_agg %>% filter(metric=="crps") %>% mutate(value=round(value, ifelse(measure=="lin", 1, 2))) %>% tidyr::spread(grouping, value) %>% select(!metric) %>% arrange(measure, model) 
scoring_results_agg %>% filter(metric=="crps", measure=="Log CRPS") %>% mutate(value=round(value, ifelse(measure=="lin", 1, 2))) %>% tidyr::spread(grouping, value) %>% select(!metric) %>% arrange(measure, model) %>% select(model, measure, LON, BCN, MAD, NYC, SF, HOMO, HET, SIM) %>% mutate(model=full_names[model])

cit_epi_phase <- ggplot(scoring_results_raw %>% filter(grouping != "LON", measure=="Linear CRPS", metric != "crps", sr %in% c("CIT"))) +
  geom_bar(aes(x = epi_phase, y = value, fill = metric),
           position = "stack",
           stat = "identity", color=NA) +
  facet_grid(grouping~model2, switch = "x", scales="free_y") +
  theme_minimal() + 
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm"),
        panel.grid=element_blank(),
        legend.position="none",
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA)) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  labs(x="", y="CRPS", title="CRPS measured on the linear forecasting scale")

ggsave("figures/crps_epi_phase_cit.png", cit_epi_phase, width=10, height=10)

gh <- ggpubr::ggarrange(ts_comb, ridge_plot, ncol=2)
ggsave("figures/figure2.png", gh, width=10, height=5)



R_inst <- function(p){
  b <- readRDS(paste0("fits/sim/fit_sim/fit_sigmoid_k=1_strat=annual_likelihood=1_epi_phase=10_",p,".rds"))
  c <- rstan::extract(b$fit)
  
  mean_R_df <- b$simulation$generation_interval_matrix %>%
    group_by(person, t_I_start) %>%
    left_join(b$simulation$generation_interval_matrix %>% group_by(infector) %>% summarise(onward_infections=n()), by = c("person" = "infector")) %>%
    select(t_I_start, person, onward_infections) %>%
    mutate(onward_infections = ifelse(is.na(onward_infections), 0, onward_infections)) %>%
    group_by(t_I_start) %>%
    summarise(mean_R=mean(onward_infections)) %>%
    mutate(rollmean_R = zoo::rollmean(mean_R, 7, fill=NA)) %>%
    filter(t_I_start < b$input_list$N+b$input_list$h)
  
  ggplot(mean_R_df, aes(x=t_I_start, y=rollmean_R)) +
    geom_point() +
    theme_minimal()
  
  j_vec <- c(0,seq_len(b$input_list$N+b$input_list$h))  # for example
  
  R_mat <- sapply(j_vec, function(jk) c$R_floor + (c$R0 - c$R_floor) * plogis(-c$delta * (jk-as.vector(c$t_change))))
  
  imm_mat <- (1-cbind(0,c$cumulative_infections[,,1]/b$param_list$num_nodes))
  
  col_ci <- apply(R_mat*imm_mat, 2, quantile, probs = c(0.025, 0.5, 0.975))[,mean_R_df$t_I_start+1]
  
  R_plotting <- mean_R_df %>% mutate(R_lower=col_ci[1,], R_median=col_ci[2,], R_upper=col_ci[3,])
  
  p <- ggplot(R_plotting, aes(x=t_I_start)) +
    geom_line(aes(y=R_median)) +
    geom_ribbon(aes(ymin=R_lower, ymax=R_upper), alpha=0.2) +
    geom_point(aes(y=rollmean_R)) +
    theme_minimal()
  
  print(p)
}

for(i in 10:20) R_inst(i)

R_case <- function(p){
  b <- readRDS(paste0("fits/sim/fit_sim/fit_sigmoid_k=1_strat=annual_likelihood=1_epi_phase=10_",p,".rds"))
  c <- rstan::extract(b$fit)
  
  gT_dist <- b$input_list$gen_weights[b$input_list$gen_weights!=0]
  tmax <- length(gT_dist)
  
  mean_R_df <- b$simulation$generation_interval_matrix %>%
    group_by(person, t_I_start) %>%
    left_join(b$simulation$generation_interval_matrix %>% group_by(infector) %>% summarise(onward_infections=n()), by = c("person" = "infector")) %>%
    select(t_I_start, person, onward_infections) %>%
    mutate(onward_infections = ifelse(is.na(onward_infections), 0, onward_infections)) %>%
    group_by(t_I_start) %>%
    summarise(mean_R=mean(onward_infections)) %>%
    mutate(rollmean_R = zoo::rollmean(mean_R, 7, fill=NA)) %>%
    filter(t_I_start < b$input_list$N+b$input_list$h-tmax)
  
  ggplot(mean_R_df, aes(x=t_I_start, y=rollmean_R)) +
    geom_point() +
    theme_minimal()
  
  j_vec <- c(0,seq_len(b$input_list$N+b$input_list$h))  # for example
  
  R_mat <- sapply(j_vec, function(jk) c$R_floor + (c$R0 - c$R_floor) * plogis(-c$delta * (jk-as.vector(c$t_change))))
  #R_vec <- R_mat[1,]
  
  imm_mat <- (1-cbind(0,c$cumulative_infections[,,1]/b$param_list$num_nodes))
  #imm_vec <- imm_mat[1,]
  
  
  #R_case <- vector(length=length(R_vec))
  #for(i in 1:length(R_case)) R_case[i] <- sum(R_vec[i:(i+tmax-1)]*gT_dist*imm_vec[i:(i+tmax-1)])
  
  R_case_mat <- matrix(nrow=nrow(R_mat), ncol=ncol(R_mat)-tmax)
  for(i in 1:nrow(R_case_mat)){
    R_vec <- R_mat[i,]
    imm_vec <- imm_mat[i,]
    for(j in 1:ncol(R_case_mat)) R_case_mat[i,j] <- sum(R_vec[j:(j+tmax-1)]*gT_dist*imm_vec[j:(j+tmax-1)])
  }
  
  col_ci <- apply(R_case_mat, 2, quantile, probs = c(0.025, 0.5, 0.975))[,mean_R_df$t_I_start+1]
  
  R_plotting <- mean_R_df %>% mutate(R_lower=col_ci[1,], R_median=col_ci[2,], R_upper=col_ci[3,])
  
  p <- ggplot(R_plotting, aes(x=t_I_start)) +
    geom_line(aes(y=R_median)) +
    geom_ribbon(aes(ymin=R_lower, ymax=R_upper), alpha=0.2) +
    geom_point(aes(y=rollmean_R)) +
    theme_minimal()
  
  print(p)
}

for(i in 10:20) R_case(i)

R_inst <- function(p){
  b <- readRDS(paste0("fits/sim/fit_sim3/fit_sigmoidP_k=1_strat=annual_likelihood=1_epi_phase=10_",p,".rds"))
  c <- rstan::extract(b$fit)
  
  mean_R_df <- b$simulation$generation_interval_matrix %>%
    group_by(person, t_I_start) %>%
    left_join(b$simulation$generation_interval_matrix %>% group_by(infector) %>% summarise(onward_infections=n()), by = c("person" = "infector")) %>%
    select(t_I_start, person, onward_infections) %>%
    mutate(onward_infections = ifelse(is.na(onward_infections), 0, onward_infections)) %>%
    group_by(t_I_start) %>%
    summarise(mean_R=mean(onward_infections)) %>%
    mutate(rollmean_R = zoo::rollmean(mean_R, 7, fill=NA)) %>%
    filter(t_I_start < b$input_list$N+b$input_list$h)
  
  ggplot(mean_R_df, aes(x=t_I_start, y=rollmean_R)) +
    geom_point() +
    theme_minimal()
  
  j_vec <- c(0,seq_len(b$input_list$N+b$input_list$h))  # for example
  
  R_mat <- sapply(j_vec, function(jk) c$R_floor + (c$R0 - c$R_floor) * plogis(-c$delta * (jk-as.vector(c$t_change))))
  
  imm_mat <- (1-cbind(0,c$cumulative_infections[,,1]/b$param_list$num_nodes))
  
  col_ci <- apply(R_mat*imm_mat, 2, quantile, probs = c(0.025, 0.5, 0.975))[,mean_R_df$t_I_start+1]
  
  R_plotting <- mean_R_df %>% mutate(R_lower=col_ci[1,], R_median=col_ci[2,], R_upper=col_ci[3,])
  
  p <- ggplot(R_plotting, aes(x=t_I_start)) +
    geom_line(aes(y=R_median)) +
    geom_ribbon(aes(ymin=R_lower, ymax=R_upper), alpha=0.2) +
    geom_point(aes(y=rollmean_R)) +
    theme_minimal()
  
  print(p)
}

for(i in 10:20) R_inst(i)
