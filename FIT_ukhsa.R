library(scoringutils)
library(cmdstanr)
library(dplyr)
library(rstan)
library(ggplot2)
library(rstan)
library(igraph)
library(magrittr)
library(EpiEstim)
library(zoo)
library(matrixStats)
library(tidyr)
library(reshape2)
library(glmnet)
library(boot)
library(epinowcast)
library(lubridate)
library(runjags)

### THIS NEEDS TO BE CHANGED
setwd("~/attenuated-renewal")

pop <- c(LON=9.8e6)   #https://worldpopulationreview.com/cities/united-kingdom/london

prop_MSM <- 0.038

MSM_pop <- pop*prop_MSM*0.5*0.5
MSM_pop["LON"] <- 0.038*pop["LON"]*0.5*0.5 # Inequalities in Sexual Health. Update on HIV and STIs in men who have sex with men in London 

## THIS SHOULD BE A VECTOR OF DAILY CASE COUNTS BY SYMPTOM ONSET DATE WITH THE FIRST ELEMENT BEING THE FIRST CASE COUNT > 0
case_counts <- c(2, 0, 0, 5, 0, 0, 2, 0, 0, 2, 0, 2, 3, 3, 3, 8, 8, 4, 
                 15, 15, 11, 23, 20, 33, 27, 11, 22, 20, 12, 47, 35, 22, 
                 38, 44, 34, 47, 39, 51, 20, 24, 26, 42, 34, 32, 23, 17, 
                 13, 17, 14, 19, 8, 5, 11, 7, 10, 8, 5, 5, 1, 8, 4, 6, 
                 5, 4, 3, 3, 2, 3, 1, 0, 1, 2, 1, 0, 2, 2, 0, 0, 2, 0, 
                 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

case_counts <- c(rbind(case_counts, case_counts))
case_counts[2] <- 0

# CHANGE THESE BASED ON WHERE YOU SAVE THE STAN FILES
lrwP <- stan_model("functions/log_random_walkP.stan")
R_random5 <- stan_model("functions/R_random5.stan")
R_biased5 <- stan_model("functions/R_biased5.stan")
expdelayP <- stan_model("functions/expdelayP.stan")
incdelayP <- stan_model("functions/incdelayP.stan")
powdelayP <- stan_model("functions/powdelayP.stan")
sigmoidP <- stan_model("functions/sigmoidP.stan")
sigmoidPS <- stan_model("functions/sigmoidPS.stan")

model_list <- list(lrwP=lrwP, expdelayP=expdelayP, incdelayP=incdelayP, powdelayP=powdelayP, R_random5=R_random5, R_biased5=R_biased5, sigmoidP=sigmoidP, sigmoidPS=sigmoidPS)

## CHANGE THIS TO THE DIRECTORY YOU WANT TO SAVE THE OUTPUTS TO
saving_dir <- "fits/LON"

fit_func_UKHSA <- function(case_counts, model_name, range_limits, strat, cutoff, tag, keep_fit=F, epi_phase=NA, decay_type=NA, dir=saving_dir){
  model <- model_list[[model_name]]
  adapt_delta <- 0.8
  
  incidence_stratified <- list()
  first_nonzero <- which(case_counts != 0)[1]
  incidence_stratified$incidence_strat <- data.frame(rounded_cases=c(rep(0,7),case_counts[first_nonzero:length(case_counts)]))
  
  incidence_stratified$n_people <- MSM_pop["LON"]
  infectious_time <- 14
  incubation_period <- 7
  k <- 1
  
  input_list <- list(N=cutoff, incidence_strat=matrix(incidence_stratified$incidence_strat[1:cutoff,], ncol=k), h=28, k=k,
                     S0=array(as.integer(incidence_stratified$n_people), dim = k), model=model_name, 
                     max_lag=(incubation_period+infectious_time), L=(incubation_period+infectious_time), 
                     gen_weights=c(rep(0,incubation_period), rep(1/infectious_time, infectious_time)), 
                     epi_phase=epi_phase, decay_type=decay_type)
  
  if(k == 1) R0 = array(4, dim=1)
  if(k != 1) R0 = rep(4, k)
  
  iter <- 1000
  seed <- 3
  chains <- 4
  
  init <- lapply(1:chains, function(i) list(R0 = R0 + rnorm(1, 0, 0.5)))
  
  fit <- sampling(model, input_list, iter=iter, chain=chains, cores=8, seed=seed, init=init, control = list(adapt_delta = adapt_delta))
  
  max_Rhat <- max(summary(fit)$summary[,"Rhat"])
  divergences <- sapply(rstan::get_sampler_params(fit, inc_warmup = FALSE), function(x) sum(x[, "divergent__"]))
  total_divergences <- sum(divergences)
  
  while(max_Rhat > 1.05 |  total_divergences >= 1){
    print(paste0("Rhat = ", max_Rhat, "     divergences = ", total_divergences), quote=F)
    
    if(max_Rhat > 1.05) iter <- iter*2
    if(total_divergences > 1) adapt_delta <- 0.99
    
    print(paste0("Retrying with ", iter, " iterations and adapt_delta = ", adapt_delta), quote=F)
    fit <- sampling(model, input_list, iter=iter, chain=chains, cores=8, seed=seed, init=init, control = list(adapt_delta = adapt_delta))
    
    max_Rhat <- max(summary(fit)$summary[,"Rhat"])
    divergences <- sapply(rstan::get_sampler_params(fit, inc_warmup = FALSE), function(x) sum(x[, "divergent__"]))
    total_divergences <- sum(divergences)
    seed <- seed + 1
    
    if(iter > 32000 | seed > 10) break
 
  }
  
  print(paste0("Rhat = ", max_Rhat, "     divergences = ", total_divergences), quote=F)
  
  output2 <- list(infections=apply(rstan::extract(fit)$infections, c(1,2), sum), 
                  forecast=apply(rstan::extract(fit)$forecast, c(1,2), sum), 
                  range_limits=range_limits, epi_phase=epi_phase,
                  input_list=input_list, tag=tag, model_name=model_name, strat=strat, cutoff=cutoff,
                  incidence_full = matrix(incidence_stratified$incidence_strat[1:(cutoff+input_list$h),], ncol=k),
                  max_Rhat = max_Rhat, divergences = divergences, total_divergences=total_divergences)
  
  output2$simulation <- NA
  output2$num_nodes <- incidence_stratified$n_people
  folder <- "cities"
  
  if (keep_fit==TRUE) output2$fit <- fit
  
  saveRDS(output2, paste0(dir,"/fit_", model_name,"_k=", k, "_strat=", strat, "_epi_phase=", epi_phase,"_", tag, ".rds"))
  
  return(output2)
}

cum_inc_raw <- cumsum(case_counts)
first_nonzero <- which(cum_inc_raw != 0)[1]
cum_inc <- c(rep(0,7), cum_inc_raw[first_nonzero:length(cum_inc_raw)])

tot_inc <- max(cum_inc)
lower <- which(cum_inc > tot_inc * 0.05)[1]
upper <- which(cum_inc > tot_inc * 0.95)[1]-28
times <- round(seq(lower, upper, length.out=10))

for(i in rev(times)){
  for(p in c("lrwP", "expdelayP", "powdelayP", "incdelayP", "sigmoidP", "R_biased5", "R_random5", "sigmoidP", "sigmoidPS")){
    print(p, quote=F)
    a0 <- fit_func_UKHSA(case_counts, model_name=p, range_limits=c(0), strat="annual", cutoff=i, tag="LON", keep_fit=TRUE, epi_phase=which(times==i), decay_type=NA) 
    cat("\n")
  }
}

files_LON <- list.files(saving_dir, full.names = TRUE) 

info_LON <- file.info(files_LON)
files_LON <- files_LON[!info_LON$isdir]

plot_list_LON     <- vector("list", length(files_LON))
lastweek_list_LON <- vector("list", length(files_LON))
Rhat_list <- vector("list", length(files_LON))

thin_every <- 10

for (ii in seq_along(files_LON)){
  print(ii)
  i <- files_LON[ii]
  # message(ii, "/", length(files), "  ", basename(i))  # optional light logging
  
  fit_output <- readRDS(i)
  N  <- fit_output$input_list$N
  h  <- fit_output$input_list$h
  Nd <- N + h
  total_infection = sum(fit_output$incidence_full)
  
  # True infections vector (length N + h) and totals
  true_vec         <- fit_output$incidence_full[seq_len(Nd)]
  #total_infection  <- sum(true_vec)
  
  # Sum across strata (keep [samples x days]) using apply margins (keeps your original intent)
  # infections: [samples x days] (keep first N cols)
  inf_m  <- apply(fit_output$infections, c(1, 2), sum)[, seq_len(N), drop = FALSE]
  # forecast:   [samples x h]
  fore_m <- apply(fit_output$forecast,   c(1, 2), sum)
  
  # Combine past + forecast → [samples x (N+h)]
  pred_mat <- cbind(inf_m, fore_m)
  
  # Thin samples by row index (same as filter(sample %% 10 == 0))
  nsamp    <- nrow(pred_mat)
  thin_idx <- which(seq_len(nsamp) %% thin_every == 0)
  pred_thin <- pred_mat[thin_idx, , drop = FALSE]
  
  # ---------- Fast 'inf_plot' (per-day quantiles across samples) ----------
  # Compute column quantiles (lower/median/upper) across rows (samples)
  lower  <- matrixStats::colQuantiles(pred_thin, probs = 0.025, na.rm = TRUE)
  median <- matrixStats::colQuantiles(pred_thin, probs = 0.50,  na.rm = TRUE)
  upper  <- matrixStats::colQuantiles(pred_thin, probs = 0.975, na.rm = TRUE)
  
  # Add params once, then build the 'simul' string from those cols
  #pl <- fit_output$simulation$param_list  # named list of scalar params (assortativity, spatial, etc.)
  
  inf_plot <- tibble(
    day          = seq_len(Nd),
    true_value   = true_vec,
    cutoff       = N,
    model        = fit_output$model_name,
    epi_phase    = fit_output$epi_phase,
    total_infection = total_infection,
    lower        = as.numeric(lower),
    median       = as.numeric(median),
    upper        = as.numeric(upper),
    cutoff2      = factor(cutoff),
    place        = "LON",
    sr           = "CIT"
  )

  # ---------- Fast 'inf_last_week' (sum over days N+22 : N+28) ----------
  week_idx <- intersect((N + 22):(N + 28), seq_len(Nd))  # guard if h < 28
  # Per-sample sums for last week (rows = samples)
  pred_week_sum <- if (length(week_idx)) {
    rowSums(pred_thin[, week_idx, drop = FALSE])
  } else {
    rep(NA_real_, nrow(pred_thin))
  }
  true_week_sum <- if (length(week_idx)) sum(true_vec[week_idx]) else NA_real_
  
  inf_last_week <- tibble(
    sample          = thin_idx,  # preserve original sample ids
    model           = fit_output$model_name,
    cutoff          = N,
    epi_phase       = fit_output$epi_phase,
    total_infection = total_infection,
    predicted       = as.numeric(pred_week_sum),
    true_value      = true_week_sum,
    place           = "LON",
    sr              = "CIT"
  ) 
  
  Rhat_div <- tibble(
    scenario        = fit_output$tag,
    model           = fit_output$model_name,
    cutoff          = N,
    epi_phase       = fit_output$epi_phase,
    Rhat            = max(rstan::summary(fit_output$fit)$summary[,"Rhat"], na.rm=T),
    divergences     = sapply(rstan::get_sampler_params(fit_output$fit, inc_warmup = FALSE), function(x) sum(x[, "divergent__"]))
  )
  
  # Accumulate (no rbind inside loop)
  plot_list_LON[[ii]]     <- inf_plot
  lastweek_list_LON[[ii]] <- inf_last_week
  Rhat_list[[i]] <- Rhat_div
}

skeleton_plot_LON      <- bind_rows(plot_list_LON)
skeleton_last_week_LON <- bind_rows(lastweek_list_LON)  
Rhat_LON <- bind_rows(Rhat_list)

## NEED THESE!
skeleton_plot_LON_censored <- skeleton_plot_LON %>% select(!true_value)
skeleton_last_week_LON_censored <- skeleton_last_week_LON %>% select(!true_value)
##

id_cols <- c("model","cutoff","epi_phase","place", "sr")

fc_LON <- as_forecast_sample(data = dplyr::distinct(skeleton_last_week_LON), predicted="predicted", observed = "true_value", sample = "sample", forecast_unit = id_cols) %>% mutate(grouping=ifelse(sr=="real", place, sr))
fc_log_LON <- transform_forecasts(fc_LON, offset = 1, append = FALSE, label = "log") 

metrics <- get_metrics(fc_LON, select = c("crps", "overprediction", "underprediction", "dispersion"))
metrics_log <- get_metrics(fc_log_LON, select = c("crps", "overprediction", "underprediction", "dispersion"))

## NEED THESE! ##
scoring_lin_LON <- as.data.frame(score(fc_LON, metrics = metrics)) %>% mutate(measure="Linear CRPS")
scoring_log_LON <- as.data.frame(score(fc_log_LON, metrics = metrics_log)) %>% mutate(measure="Log CRPS")
##

