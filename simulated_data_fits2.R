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

setwd("~/attenuated-renewal")

ranges_generator <- function(simulation, range_limits, strat){
  if(strat=="annual") incidence_strat <- simulation$incidence_matrix_annual
  if(strat=="infectious_period") incidence_strat <- simulation$incidence_matrix_nw
  if(strat=="recurrent") incidence_strat <- simulation$incidence_matrix_recurrent
  if(strat=="weighted") incidence_strat <- simulation$incidence_matrix_weighted
  
  max_value <- max(as.numeric(colnames(incidence_strat)))
  
  range_limits <- range_limits[range_limits<=max_value]
  
  ranges <- Map(
    function(start, end) if (is.na(end)) start:max_value else start:(end - 1),
    range_limits,
    c(range_limits[-1], NA) #+ ifelse(strat=="recurrent", 1, 0)
  )
  
  if (strat %in% c("recurrent", "annual", "infectious_period", "weighted")) ranges <- lapply(ranges, function(x) x + 1)
  
  if(length(ranges)>1){
    names(ranges) <- c(
      paste0(range_limits[-length(range_limits)], "-", range_limits[-1] - 1),
      paste0(range_limits[length(range_limits)], "+")
    )
  }
  
  return(ranges)
}

incidence_strat_generator <- function(simulation, range_limits, strat){
  if(strat=="annual"){
    incidence_strat <- simulation$incidence_matrix_annual
    degrees <- simulation$degrees_annual
  } 
  if(strat=="infectious_period"){
    incidence_strat <- simulation$incidence_matrix_nw
    degrees <- simulation$degrees
  } 
  if(strat=="recurrent"){
    incidence_strat <- simulation$incidence_matrix_recurrent
    degrees <- simulation$degrees_recurrent
  } 
  if(strat=="weighted"){
    incidence_strat <- simulation$incidence_matrix_weighted
    degrees <- simulation$degrees_weighted
  } 
  if(strat=="infected_by"){
    days <- data.frame(t_E_start=simulation$incidence$t_E)
    incidence_strat_raw <- simulation$generation_interval_matrix %>%
      group_by(t_E_start, infected_by) %>%
      summarise(count=n()) %>%
      filter(!is.na(infected_by)) %>%
      tidyr::spread(infected_by, count)  
    
    incidence_strat <- left_join(days, incidence_strat_raw, by="t_E_start") %>% select(-t_E_start) %>%
      replace_na(list(main = 0, ot = 0))
  }
  
  ranges <- ranges_generator(simulation, range_limits, strat)
  
  summed_df <- as.matrix(as.data.frame(lapply(ranges, function(cols) rowSums(incidence_strat[, cols, drop = FALSE]))))
  colnames(summed_df) <- names(ranges)
  
  n_people <- sapply(ranges, function(range) sum(degrees %in% range))
  
  # edges <- simulation$generation_interval_matrix %>%
  #   group_by(t_E_start) %>%
  #   summarise(connections = sum(n_connections_weighted),
  #             exposed=n()) %>%
  #   mutate(cum_connections = cumsum(connections),
  #          cum_exposed = cumsum(exposed)) %>%
  #   complete(t_E_start = 0:max(t_E_start, na.rm=T)) %>%
  #   fill(everything(), .direction = "down") %>%
  #   mutate(average_edges=cum_connections/cum_exposed) %>%
  #   filter(t_E_start > 0) %>%
  #   select(average_edges)
  # 
  return(list(incidence_strat=summed_df,
              n_people=n_people, 
              n=length(ranges)))#,
  #average_edges=average_edges))
}

lrwP <- stan_model("functions/log_random_walkP.stan")
R_biased5 <- stan_model("functions/R_biased5.stan")
R_random5 <- stan_model("functions/R_random5.stan")
expdelayP <- stan_model("functions/expdelayP.stan")
incdelayP <- stan_model("functions/incdelayP.stan")
powdelayP <- stan_model("functions/powdelayP.stan")
sigmoidP <- stan_model("functions/sigmoidP.stan")
sigmoidPS <- stan_model("functions/sigmoidPS.stan")

model_list <- list(lrwP=lrwP, expdelayP=expdelayP, incdelayP=incdelayP, powdelayP=powdelayP, R_random5=R_random5, R_biased5=R_biased5, sigmoidP=sigmoidP, sigmoidPS=sigmoidPS)

NYC_data_raw <- read.csv("data/NYC_cases.csv") %>% mutate(date=as.Date(diagnosis_date, format="%d/%m/%Y")) %>% select(date, count) %>% rename(confirm=count)

df <- NYC_data_raw %>% mutate(dow = factor(weekdays(date), levels = c("Monday","Tuesday","Wednesday","Thursday","Friday","Saturday","Sunday")))
m <- glm(confirm ~ dow, family = quasipoisson(link = "log"), data = df)

# Get multiplicative weekday factors relative to Monday
co <- coef(m); levs <- levels(df$dow)
eff <- setNames(numeric(length(levs)), levs)
eff["Monday"] <- exp(co["(Intercept)"])
for (lv in levs[-1]) eff[lv] <- exp(co["(Intercept)"] + co[paste0("dow", lv)])

rel <- eff / eff["Monday"]                 # Monday = 1
df$confirm_adj <- round(df$confirm / rel[df$dow]) # de-weekended series

NYC_data <- df %>% select(date, confirm_adj) %>% rename(NYC=confirm_adj)
SF_data <- read.csv("data/SF_cases.csv") %>% mutate(date=as.Date(episode_date, format="%d/%m/%Y"), SF=new_cases) %>% select(date, SF)
SPN_data <- readxl::read_xlsx("data/mpox_BCN_MAD.xlsx", sheet="Incidence") %>% mutate(date=as.Date(Date), MAD=Madrid, BCN=Barcelona) %>% select(date, MAD, BCN)
skeleton <- data.frame(date=seq(min(c(NYC_data$date, SF_data$date, SPN_data$date)), max(c(NYC_data$date, SF_data$date, SPN_data$date)), by="day"))

dat <- full_join(skeleton, NYC_data, by="date") %>%
  full_join(., SF_data, by="date") %>%
  full_join(., SPN_data, by="date") %>%
  mutate(across(c(NYC, SF, MAD, BCN), ~ replace_na(., 0))) %>%
  filter(row_number() < 365)

ggplot(dat[1:200,] %>% tidyr::gather("city", "value", 2:5), aes(x=date, y=value)) +
  geom_line(aes(color=city)) +
  theme_bw() +
  theme(panel.grid = element_blank())

pop <- c(BCN=1.73e6,  #https://portaldades.ajuntament.barcelona.cat/en/statistics/yzlntdm2fs
         MAD=3.4e6,   #https://www.citypopulation.de/en/spain/madrid/madrid/28079__madrid/
         NYC=7.94e6,  #https://worldpopulationreview.com/us-cities
         SF=0.768e6,
         LON=9.8e6)   #https://worldpopulationreview.com/cities/united-kingdom/london

prop_MSM <- 0.038

MSM_pop <- pop*prop_MSM

MSM_pop["LON"] <- 0.038*pop["LON"]*0.5*0.5        # Inequalities in Sexual Health. Update on HIV and STIs in men who have sex with men in London 
MSM_pop["NYC"] <- round(397399*0.5)               # Estimating Population Sizes of Men Who Have Sex with Men in the United States # chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://depts.washington.edu/hivtcg/presentations/uploads/35/estimating_population_sizes_of_men_who_have_sex_with_men_in_the_united_states.pdf?utm_source=chatgpt.com
MSM_pop["SF"] <- round(145972*0.5)                # Estimating Population Sizes of Men Who Have Sex with Men in the United States # chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://depts.washington.edu/hivtcg/presentations/uploads/35/estimating_population_sizes_of_men_who_have_sex_with_men_in_the_united_states.pdf?utm_source=chatgpt.com
MSM_pop["BCN"] <- round(pop["BCN"]*0.5*0.5*0.053) # https://pmc.ncbi.nlm.nih.gov/articles/PMC4594901/#:~:text=A%20total%20of%201560%20cases,was%20only%2011.9/100%2C000%20inhabitants.
MSM_pop["MAD"] <- round(pop["MAD"]*0.5*0.5*0.053) # https://pmc.ncbi.nlm.nih.gov/articles/PMC4594901/#:~:text=A%20total%20of%201560%20cases,was%20only%2011.9/100%2C000%20inhabitants.

fit_func <- function(simulation, model_name, range_limits, strat, cutoff, tag, keep_fit=F, epi_phase=NA, decay_type=NA, dir="", iter=1000, chains=4, adapt_delta=0.9){
  model <- model_list[[model_name]]
  
  likelihood_type <- ifelse(length(range_limits)==1, 1, 2)
  
  if(length(simulation)>1){
    infectious_time <- simulation$param_list$infectious_time
    incubation_period <- simulation$param_list$incubation_period
    incidence_stratified <- incidence_strat_generator(simulation=simulation, range_limits=range_limits, strat=strat)  
    k <- incidence_stratified$n
    rec_daily <- simulation$param_list$rec_daily
  }
  
  if(length(simulation)==1){
    incidence_stratified <- list()
    first_nonzero <- which(dat[simulation] != 0)[1]
    incidence_stratified$incidence_strat <- data.frame(rounded_cases=c(rep(0,7),dat[simulation][first_nonzero:nrow(dat[simulation]),]))
    #incidence_stratified$incidence_strat <- country_data %>% filter(country==simulation) %>% mutate(rounded_cases=round(zoo::rollmean(new_cases,7, fill=NA))) %>% filter(!is.na(rounded_cases)) %>% select(rounded_cases) 
    incidence_stratified$n_people <- MSM_pop[simulation]
    infectious_time <- 14
    incubation_period <- 7
    k <- 1
  }
  
  input_list <- list(N=cutoff, incidence_strat=matrix(incidence_stratified$incidence_strat[1:cutoff,], ncol=k), h=28, k=k,
                     S0=array(as.integer(incidence_stratified$n_people), dim = k), model=model_name, 
                     max_lag=(incubation_period+infectious_time), L=(incubation_period+infectious_time), 
                     gen_weights=c(rep(0,incubation_period), rep(1/infectious_time, infectious_time)), 
                     epi_phase=epi_phase, decay_type=decay_type)
  
  if(k == 1) R0 = array(4, dim=1)
  if(k != 1) R0 = rep(4, k)
  
  init_raw <- list(R0 = R0)
  init <- replicate(chains, init_raw, simplify = FALSE)
  
  print(paste0("Model name = ", model_name, "     Epi phase = ", epi_phase, "     Tag = ", tag), quote=F)
  
  p <- proc.time()
  fit <- sampling(model, input_list, iter=iter, chains=chains, cores=4, seed=(15), init=init, control = list(adapt_delta = adapt_delta)) 
  
  max_Rhat <- max(summary(fit)$summary[,"Rhat"])
  divergences <- sapply(rstan::get_sampler_params(fit, inc_warmup = FALSE), function(x) sum(x[, "divergent__"]))
  total_divergences <- sum(divergences)
  
  print(paste0("Rhat = ", round(max_Rhat, 2), " Divergences = ", total_divergences), quote=F)
  
  if(max_Rhat>1.05) iter <- iter*2
  if(total_divergences > 1) adapt_delta = 0.95
  
  if(max_Rhat>1.05 | total_divergences >= 1){
    print(paste0("Rerunning with iter = ", iter, " and adapt_delta = ", adapt_delta), quote=F)
    fit <- sampling(model, input_list, iter=iter, chains=chains, cores=4, seed=(16), init=init, control = list(adapt_delta = adapt_delta)) 
    max_Rhat <- max(summary(fit)$summary[,"Rhat"])
    divergences <- sapply(rstan::get_sampler_params(fit, inc_warmup = FALSE), function(x) sum(x[, "divergent__"]))
    total_divergences <- sum(divergences)
  } 
  
  print(paste0("Rhat = ", round(max_Rhat, 2), " Divergences = ", total_divergences), quote=F)
  
  if(max_Rhat>1.05) iter <- iter*2
  if(total_divergences > 1) adapt_delta = 0.99
  
  if(max_Rhat>1.05 | total_divergences >= 1){
    print(paste0("Rerunning with iter = ", iter, " and adapt_delta = ", adapt_delta), quote=F)
    fit <- sampling(model, input_list, iter=iter, chains=chains, cores=4, seed=(17), init=init, control = list(adapt_delta = adapt_delta))
    max_Rhat <- max(summary(fit)$summary[,"Rhat"])
    divergences <- sapply(rstan::get_sampler_params(fit, inc_warmup = FALSE), function(x) sum(x[, "divergent__"]))
    total_divergences <- sum(divergences)
  }
  
  max_Rhat <- max(summary(fit)$summary[,"Rhat"])
  divergences <- sapply(rstan::get_sampler_params(fit, inc_warmup = FALSE), function(x) sum(x[, "divergent__"]))
  
  print(paste0("Rhat = ", round(max_Rhat, 2), " Divergences = ", total_divergences), quote=F)
  
  if(max_Rhat>1.05) iter <- iter*2
  if(total_divergences > 1) adapt_delta = 0.995
  
  if(max_Rhat>1.05 | total_divergences >= 1){
    print(paste0("Rerunning with iter = ", iter, " and adapt_delta = ", adapt_delta), quote=F)
    fit <- sampling(model, input_list, iter=iter, chains=chains, cores=4, seed=(18), init=init, control = list(adapt_delta = adapt_delta))
    max_Rhat <- max(summary(fit)$summary[,"Rhat"])
    divergences <- sapply(rstan::get_sampler_params(fit, inc_warmup = FALSE), function(x) sum(x[, "divergent__"]))
    total_divergences <- sum(divergences)
  }
  
  max_Rhat <- max(summary(fit)$summary[,"Rhat"])
  divergences <- sapply(rstan::get_sampler_params(fit, inc_warmup = FALSE), function(x) sum(x[, "divergent__"]))
  
  print(paste0("Rhat = ", round(max_Rhat, 2), " Divergences = ", total_divergences), quote=F)
  
  time_taken <- round((proc.time()-p)["elapsed"],2)
  
  print(paste0("time taken = ", time_taken), quote=F)
  if(total_divergences>1) cat("\n\n\n\n\n\n\n")
  
  output2 <- list(infections=apply(rstan::extract(fit)$infections, c(1,2), sum), 
                  forecast=apply(rstan::extract(fit)$forecast, c(1,2), sum), 
                  range_limits=range_limits,
                  input_list=input_list, tag=tag, model_name=model_name, strat=strat, cutoff=cutoff,
                  incidence_full = matrix(incidence_stratified$incidence_strat[1:(cutoff+input_list$h),], ncol=k), 
                  max_Rhat = max_Rhat, divergences = divergences, total_divergences = total_divergences)
  
  if(length(simulation)==1){
    output2$simulation <- NA
    output2$num_nodes <- 100000
    folder <- "cities"
  }
  
  if(length(simulation)>1){
    output2$simulation <- simulation
    output2$num_nodes <- length(simulation$degrees)
    output2$param_list <- simulation$param_list
    folder <- "sim"
  }
  
  output2$time_taken <- time_taken
  output2$file_name <- paste0("fits/", folder,"/", dir, "/fit_", model_name,"_k=", k, "_strat=", strat, "_likelihood=", likelihood_type, "_epi_phase=", epi_phase,"_", tag, ".rds")
  
  if (keep_fit==TRUE) output2$fit <- fit
  
  saveRDS(output2, output2$file_name)
  
  return(output2)
}

k <- 1
file_list <- list()
R_hat_list <- list()
model_name_list <- list()
divergences_list <- list()
time_taken <- list()

for(j in c(list.files("simulation/sim_2/"))[9:28]){
  simulation <- readRDS(paste0("simulation/sim_2/",j))$outbreak
  
  cum_inc <- cumsum(simulation$incidence$incidence)
  tot_inc <- max(cum_inc)
  lower <- which(cum_inc > tot_inc * 0.05)[1]
  upper <- which(cum_inc > tot_inc * 0.95)[1]-28
  times <- round(seq(lower, upper, length.out=10))
  
  tag <- which(j==list.files("simulation/sim_2/")[1:100])
  
  for(i in times){
    for(p in c("expdelayP")){#"lrwP", "R_random5", "R_biased5", "expdelayP", "incdelayP", "powdelayP", "sigmoidP")){
      a5 <- fit_func(simulation, model_name=p, range_limits=c(0), strat="annual", cutoff=i, tag=tag, keep_fit=TRUE, epi_phase=which(times==i), decay_type=NA, dir="fit_sim3", iter=1000, chains=4)
      file_list[[k]] <- a5$file_name
      R_hat_list[[k]] <- a5$max_Rhat
      model_name_list[[k]] <- p
      divergences_list[[k]] <- a5$total_divergences
      time_taken[[k]] <- a5$time_taken
      k <- k+1
      cat("\n")
    }
  }
}

l <- 1

for(j in c(list.files("simulation/sim_homo/")[1])){
  simulation <- readRDS(paste0("simulation/sim_homo/",j))$outbreak
  
  cum_inc <- cumsum(simulation$incidence$incidence)
  tot_inc <- max(cum_inc)
  lower <- which(cum_inc > tot_inc * 0.05)[1]
  upper <- which(cum_inc > tot_inc * 0.95)[1]-28
  times <- round(seq(lower, upper, length.out=10))
  
  tag <- which(j==list.files("simulation/sim_homo/"))
  
  for(i in times){
    for(p in c("R_biased5")){#"lrwP", "R_random5", "R_biased5", "expdelayP", "incdelayP", "powdelayP", "sigmoidPS")){
      a5 <- fit_func(simulation, model_name=p, range_limits=c(0), strat="annual", cutoff=i, tag=tag, keep_fit=TRUE, epi_phase=which(times==i), decay_type=NA, dir="fit_homo") 
      file_list[[l]] <- a5$file_name
      R_hat_list[[l]] <- a5$max_Rhat
      model_name_list[[l]] <- p
      divergences_list[[l]] <- a5$total_divergences
      l <- l+1
      cat("\n")
    }
  }
}

m <- 1
file_list_het <- list()
R_hat_list_het <- list()
model_name_list_het <- list()
divergences_list_het <- list()
time_taken_het <- list()

for(j in c(list.files("simulation/sim_het/"))[8:30]){
  #print(j)
  simulation <- readRDS(paste0("simulation/sim_het/",j))$outbreak
  
  cum_inc <- cumsum(simulation$incidence$incidence)
  tot_inc <- max(cum_inc)
  lower <- which(cum_inc > tot_inc * 0.05)[1]
  upper <- which(cum_inc > tot_inc * 0.95)[1]-28
  times <- round(seq(lower, upper, length.out=10))
  
  tag <- which(j==list.files("simulation/sim_het/"))
  
  for(i in times){
    for(p in c("lrwP", "R_random5", "R_biased5", "expdelayP", "incdelayP", "powdelayP", "sigmoidP")){
      a5 <- fit_func(simulation, model_name=p, range_limits=c(0), strat="annual", cutoff=i, tag=tag, keep_fit=TRUE, epi_phase=which(times==i), decay_type=NA, dir="fit_het")
      file_list_het[[m]] <- a5$file_name
      R_hat_list_het[[m]] <- a5$max_Rhat
      model_name_list_het[[m]] <- p
      divergences_list_het[[m]] <- a5$total_divergences
      m <- m+1
      cat("\n")
    }
  }
}

for(j in c("BCN", "SF", "NYC", "MAD")){
  simulation <- j#readRDS(paste0("simulation/",j))$outbreak
  
  cum_inc_raw <- cumsum(dat[,j])
  first_nonzero <- which(cum_inc_raw != 0)[1]
  cum_inc <- c(rep(0,7), cum_inc_raw[first_nonzero:length(cum_inc_raw)])
  
  tot_inc <- max(cum_inc)
  lower <- which(cum_inc > tot_inc * 0.05)[1]
  upper <- which(cum_inc > tot_inc * 0.95)[1]-28
  times <- round(seq(lower, upper, length.out=10))
  
  tag <- j
  
  for(i in times){
    for(p in c("lrwP", "R_random5", "R_biased5", "expdelayP", "incdelayP", "powdelayP", "sigmoidP")){
      a0 <- fit_func(simulation, model_name=p, range_limits=c(0), strat="annual", cutoff=i, tag=tag, keep_fit=TRUE, epi_phase=which(times==i), decay_type=NA)
    }
  }
}

for(j in list.files("simulation/")[1:100]){
  simulation <- readRDS(paste0("simulation/",j))$outbreak
  tag <- which(j==list.files("simulation/")[1:100])
  
  if(j == list.files("simulation/")[1]) skeleton <- data.frame(tag=tag, simulation$param_list)
  else skeleton <- rbind(skeleton, data.frame(tag=tag, simulation$param_list))
}  

plotting_func <- function(fit_output){
  input <- fit_output$input_list
  fit <- fit_output$fit
  incidence_full <- fit_output$incidence_full
  
  df1 <- cbind(data.frame(x=1:input$N, matrixStats::colQuantiles(apply(rstan::extract(fit)$infections,c(1,2),sum), probs=c(0.025, 0.5, 0.975))[1:input$N,], incidence=c(rowSums(as.matrix(input$incidence_strat))))) %>% magrittr::set_colnames(c("x", "lower", "median", "upper", "incidence"))
  df2 <- cbind(data.frame(x=(input$N+1):(input$N+input$h), matrixStats::colQuantiles(apply(rstan::extract(fit)$forecast,  c(1,2),sum), probs=c(0.025, 0.5, 0.975)), incidence=NA)) %>% magrittr::set_colnames(c("x", "lower", "median", "upper", "incidence"))
  
  df <- rbind(df1, df2) %>% mutate(incidence=rowSums(incidence_full))
  
  p <- ggplot(df, aes(x=x)) + geom_point(aes(y=incidence)) + geom_line(aes(y=median)) + geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) + theme_bw()
  
  inf <- rstan::extract(fit)$infections
  forecast <- rstan::extract(fit)$forecast
  
  q_array <- apply(as.array(inf), c(2, 3), quantile, probs = c(0.025, 0.5, 0.975))
  
  df2 <- melt(q_array, varnames = c("quantile", "time", "stratum")) %>%
    pivot_wider(names_from = quantile, values_from = value) %>%
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
  
  print(p)
  print(q)
  
  r <- ggpubr::ggarrange(p, q, nrow=2)
  print(r)
  
}

### CITIES AGG
files <- list.files("fits/cities/", full.names = TRUE) 

info <- file.info(files)
files <- files[!info$isdir]

sorted_files <- files[order(info$mtime[!info$isdir], decreasing = TRUE)]

for(i in sorted_files[1:8]){
  print(i)
  print(which(sorted_files==i)) 
  
  fit_output <- readRDS(i)
  timesteps <- fit_output$input_list$N
  projection <- fit_output$input_list$h
  
  true_infections <- data.frame(day=1:(fit_output$input_list$N+fit_output$input_list$h)) %>%
    mutate(true_value=fit_output$incidence_full[1:(fit_output$input_list$N+fit_output$input_list$h)]) %>%
    mutate(total_infection=sum(true_value))
  
  inf_samples <- cbind(apply(fit_output$infections, c(1,2), sum)[,1:timesteps], apply(fit_output$forecast, c(1,2), sum)) %>%
    as.data.frame() %>%
    magrittr::set_colnames(1:ncol(.)) %>%
    mutate(sample=row_number()) %>%
    filter(sample %% 10 == 0) %>%
    tidyr::gather("day", "prediction", 1:(ncol(.)-1)) %>%
    mutate(day=as.numeric(day)) %>%
    left_join(., true_infections, by="day") %>%
    mutate(cutoff=fit_output$input_list$N, 
           model=fit_output$model_name,
           epi_phase=fit_output$input_list$epi_phase[1],
           strata=paste0(fit_output$strat,"_",fit_output$input_list$k)) %>% 
    mutate(place=fit_output$tag)
  
  inf_plot <- inf_samples %>% dplyr::select(-sample) %>% 
    dplyr::group_by(day, true_value, place, cutoff, model, epi_phase, total_infection) %>%
    summarise(lower=quantile(prediction, probs=0.025),
              median=quantile(prediction, probs=0.5),
              upper=quantile(prediction, probs=0.975)) %>%
    mutate(cutoff2=factor(as.character(cutoff)))%>%
    mutate(epi_phase=as.factor(epi_phase)) 
  
  inf_last_week <- inf_samples %>% filter((day-cutoff) %in% c(22:28)) %>%
    group_by(sample, model, place, cutoff, epi_phase, total_infection) %>%
    summarise(prediction=sum(prediction), 
              true_value=sum(true_value)) %>%
    mutate(epi_phase=as.factor(epi_phase)) 
  
  
  if(i == sorted_files[1]){
    skeleton_plot_cit <- inf_plot
    skeleton_last_week_cit <- inf_last_week
  }       
  else{
    skeleton_plot_cit <- rbind(skeleton_plot_cit, inf_plot)
    skeleton_last_week_cit <- rbind(skeleton_last_week_cit, inf_last_week)
  } 
  cat("\n")
}

saveRDS(skeleton_plot_cit, "fits/outputs/skeleton_plot_cit.rds")
saveRDS(skeleton_last_week_cit, "fits/outputs/skeleton_last_week_cit.rds")


### SIM AGG


#skeleton_plot_sim2 <- skeleton_plot_sim[!duplicated(skeleton_plot_sim), ]
#skeleton_last_week_sim2 <- skeleton_last_week_sim[!duplicated(skeleton_last_week_sim), ]

skeleton_plot_sim2 <- readRDS("fits/outputs/skeleton_plot_sim_all6.rds")
skeleton_last_week_sim2 <- readRDS("fits/outputs/skeleton_last_week_sim_all6.rds")

skeleton_plot_s <- bind_rows(skeleton_plot_sim, skeleton_plot_sim2)
skeleton_last_week_s <- bind_rows(skeleton_last_week_sim, skeleton_last_week_sim2)

saveRDS(skeleton_plot_sim, "fits/outputs/skeleton_plot_sim_all7.rds")
saveRDS(skeleton_last_week_sim, "fits/outputs/skeleton_last_week_sim_all7.rds")

skeleton_plot_sim <- readRDS("fits/outputs/skeleton_plot_sim.rds") %>% mutate(epi_phase=factor(epi_phase))
skeleton_last_week_sim <- readRDS("fits/outputs/skeleton_last_week_sim.rds") %>% mutate(epi_phase=factor(epi_phase))

plotting_func(readRDS("fits/fit_sat_k=1_strat=annual_likelihood=1_cutoff=315_1.rds"))

skeleton_stats <- skeleton %>% group_by(scenario, num_nodes, infectious_time, transmission_prob, assortativity, ass_v_param, spatial, isolation_time_mean, isolation_time_sd, directionality, incubation_period) %>%
  summarise(true_value=max(true_value))

ggplot(skeleton %>% mutate(trans=round(transmission_prob*infectious_time*num_nodes,2)), aes(x=day, y=true_value)) +
  geom_line(aes(color=trans)) +
  facet_wrap(~trans)

ggplot(skeleton_plot %>% filter(scenario==92), aes(x=day, fill=factor(epi_phase))) +
  geom_point(aes(y=true_value)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
  geom_line(aes(y=median)) +
  theme_bw() +
  #scale_y_log10() +
  facet_grid(model~scenario+cutoff, scales="free") +
  theme(panel.grid=element_blank()) +
  geom_vline(aes(xintercept=cutoff, color=epi_phase), linetype="dashed")

ggplot(skeleton_plot %>% filter(scenario==92, median>1, true_value>0), aes(x=day, fill=factor(epi_phase))) +
  geom_point(aes(y=true_value)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
  geom_line(aes(y=median)) +
  geom_vline(aes(xintercept=cutoff), linetype="dashed") +
  theme_bw() +
  scale_y_log10() +
  facet_grid(model~scenario+cutoff, scales="free_y") +
  theme(panel.grid=element_blank())

a <- transform_forecasts(dplyr::distinct(skeleton_last_week), offset=1) %>% filter(scale=="log")

scoring_results_raw <- scoringutils::score(a, by=c("prediction", "true_value"))  %>%
  mutate(model_type = model_type_dict[model])

comparison_func <- function(scenario_interest, model_list, epi_phase_list, cutoff_graph=F){
  if(cutoff_graph==T) skeleton_plot <- skeleton_plot %>% filter(cutoff<=day+10)
  
  p1 <- ggplot(skeleton_plot %>% filter(scenario %in% scenario_interest, model %in% model_list, epi_phase %in% epi_phase_list, true_value>0), aes(x=day)) +
    geom_point(aes(y=true_value)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
    geom_line(aes(y=median)) +
    geom_vline(aes(xintercept=cutoff), linetype="dashed") +
    theme_bw() +
    #scale_y_log10() +
    facet_grid(model+epi_phase~scenario, scales="free") +
    theme(panel.grid=element_blank())
  
  p2 <- ggplot(skeleton_plot %>% filter(scenario %in% scenario_interest, model %in% model_list, epi_phase %in% epi_phase_list, median > 1, true_value>0), aes(x=day)) +
    geom_point(aes(y=true_value)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
    geom_line(aes(y=median)) +
    geom_vline(aes(xintercept=cutoff), linetype="dashed") +
    theme_bw() +
    scale_y_log10() +
    facet_grid(model+epi_phase~scenario, scales="free") +
    theme(panel.grid=element_blank())
  
  q <- scoring_results_raw %>% filter(scenario %in% scenario_interest, model %in% model_list, epi_phase %in% epi_phase_list) %>% arrange(scenario, crps)
  
  print(ggpubr::ggarrange(p1, p2, nrow=1))
  print(q)
  
}

comparison_func(scenario_interest=c(52:57), model_list=c("expdelay2", "expdelay", "R_biased", "biaslrw"), epi_phase_list=3, cutoff_graph = F)
comparison_func(scenario_interest=c(80:89), model_list=c("expdelay", "lrw", "rw", "biaslrw"), epi_phase_list=2, cutoff_graph = F)
comparison_func(scenario_interest=c(10:19), model_list=c("expdelay", "inc2delay", "biaslrw", "rw"), epi_phase_list=c(2), cutoff_graph = F)

scoring_results <- scoring_results_raw %>%
  group_by(model, epi_phase, cutoff, scenario, assortativity, dd_upper, spatial, dd_param, transmission_prob,
           infectious_time, directionality, ass_v_param, incubation_period, rec_daily, versatile) %>%
  summarise(crps=median(crps), log_score=median(log_score), ae_median=median(ae_median), se_mean=median(se_mean), bias=mean(bias), total_infection=mean(total_infection)) %>%
  group_by("model", "scenario") %>%
  mutate(model_type = model_type_dict[model]) %>%
  mutate(model=factor(model, levels=c("rw", "lrw", "biaslrw", "R_biased", "phenom", 
                                      "exp", "expdelay", "expdelay2", "lin", "lindelay", "inc", "incdelay", "inc2delay", 
                                      "comp"))) %>%
  mutate(epi_phase=factor(epi_phase, levels=c(1,2,3,4,5)))

dat <- scoring_results %>% filter(model=="expdelay") %>% ungroup() %>% 
  select(assortativity, dd_upper, dd_upper, spatial, dd_param, transmission_prob,
         infectious_time, directionality, ass_v_param, incubation_period, rec_daily, versatile, crps) %>% 
  na.omit()

X <- model.matrix(crps ~ ., dat)[, -1]
y <- dat$crps

cvfit <- cv.glmnet(X, y, alpha = 1)
fit_coefs <- coef(cvfit, s = "lambda.min") 



ggplot(scoring_results, aes(y=crps, x=model)) +
  geom_col(position=position_dodge()) +
  facet_grid(paste0("epi_phase=",epi_phase)~scenario, scales="free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))

ggplot(scoring_results, aes(y=log_score, x=model)) +
  geom_col(position=position_dodge()) +
  facet_grid(paste0("epi_phase=",epi_phase)~scenario, scales="free_y") +
  theme_bw() 

ggplot(scoring_results %>% filter(model=="expdelay"), aes(x=assortativity, y=log_score)) +
  geom_point() +
  theme_bw()

boot_mean <- function(data, idx) mean(data[idx], na.rm = TRUE)


ranked_model_func <- function(scoring_results, total_min, filters=list()){
  ranked_models_raw <- scoring_results %>% 
    filter(model != "biaslrw") %>%
    filter(total_infection > total_min) %>%
    group_by(epi_phase, scenario) %>%
    arrange(epi_phase, scenario, .by_group = TRUE) 
  
  for (nm in names(filters)) {
    bounds <- filters[[nm]]
    ranked_models_raw <- ranked_models_raw %>%
      filter(between(.data[[nm]], bounds[1], bounds[2]))
  }
  
  ranked_models <- ranked_models_raw %>%
    mutate(
      rank_log_score = row_number(log_score),
      rank_crps      = row_number(crps),
      rank_ae_median = row_number(ae_median),
      rank_se_mean   = row_number(se_mean)
    ) %>%
    ungroup() 
  
  bias_table <- ranked_models %>% filter(epi_phase != 5) %>%
    group_by(epi_phase, model) %>%
    summarise(mean_bias = mean(bias, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = model, values_from = mean_bias)
  
  o <- bias_table %>%
    tidyr::pivot_longer(-epi_phase, names_to = "model", values_to = "bias") %>%
    ggplot(aes(x = model, y = factor(epi_phase), fill = bias)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(x = "Model", y = "Epidemic Phase", fill = "Mean Bias") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ranked_models_all <- ranked_models %>%
    tidyr::pivot_longer(cols = starts_with("rank_"),
                        names_to = "metric",
                        values_to = "rank") %>%
    group_by(model, epi_phase, metric) %>%
    summarise(boot = list(boot(rank, boot_mean, R = 2000)),
              .groups = "drop") %>%
    mutate(mean_rank = map_dbl(boot, ~ mean(.$t, na.rm = TRUE)),
           lower_ci  = map_dbl(boot, ~ quantile(.$t, 0.025, na.rm = TRUE)),
           upper_ci  = map_dbl(boot, ~ quantile(.$t, 0.975, na.rm = TRUE))) %>%
    select(-boot) %>%
    mutate(metric = sub("^rank_", "", metric),
           model_type = model_type_dict[as.character(model)])
  
  message("Total fits = ", length(unique(ranked_models$scenario)))
  
  p <- ggplot(ranked_models_all, aes(x=model, y=mean_rank, fill=model_type)) + 
    geom_col() +
    geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=0.5) +
    facet_grid(metric~epi_phase, scales = "free_y") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
      legend.position = "none",
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold")
    ) +
    labs(y="Mean rank with bootstrapped 95% CrI", x="Model")
  
  print(p)
  
  ggsave(filename = "figures/ranked_models.png", plot = p, width = 12, height = 12, dpi = 300)
  
  return(ranked_models_all)
}

ranked_model_func(scoring_results, total_min=500, filters=list())

#






ske <- data.frame(simulation=numeric(), assortativity=numeric())
simuls <- list.files("simulation/", full.names=T)
for(i in 1:(length(simuls)-1)){
  ske <- rbind(ske, data.frame(simulation=i, assortativity=readRDS(simuls[i])$g$assortativity_kernel))
}

ranked_models_av <- ranked_models %>% group_by(model, epi_phase) %>%
  mutate(rank_count=rank*count) %>%
  summarise(average=sum(rank_count)/sum(count))

ggplot(ranked_models_av, aes(x=model, y=average)) + 
  geom_col() +
  facet_wrap(~epi_phase)

ggplot(ranked_models %>% filter(rank %in% c(1,2,3)), aes(x=model, y=count)) +
  geom_col() +
  facet_grid(epi_phase~rank)

scoring_results_agg <- scoring_results_raw %>%
  group_by(model, strata, scenario) %>%
  filter((day-cutoff) %in% c(28)) %>%
  summarise(crps=mean(crps), log_score=mean(log_score)) %>%
  mutate(border=ifelse(strata=="none", "1", "0"))

ggplot(scoring_results_agg, aes(y=crps, x=model, fill=strata)) +
  geom_col(position=position_dodge()) +
  facet_grid(~scenario) +
  theme_bw() 

scoring_results_agg2 <- scoring_results_raw %>%
  group_by(model, strata) %>%
  filter((day-cutoff) %in% c(28)) %>%
  summarise(crps=mean(crps), log_score=mean(log_score)) %>%
  mutate(border=ifelse(strata=="none", "1", "0"))

ggplot(scoring_results_agg2, aes(y=crps, x=model, fill=strata)) +
  geom_col(position=position_dodge()) +
  theme_bw() 


a <- readRDS("fits/fit_sat_k=5_strat=recurrent_likelihood=2_cutoff=75.rds")




plotting_func(readRDS("fits/fit_sat_k=4_strat=recurrent_likelihood=2_cutoff=100_ass.rds"))



true_infections <- data.frame(day=1:(fit_output$input_list$N+fit_output$input_list$h)) %>%
  mutate(true_value=simulation$incidence$incidence[1:(fit_output$input_list$N+fit_output$input_list$h)])

inf_samples <- apply(rstan::extract(fit_output$fit)$infections, c(1,2), sum) %>%
  as.data.frame() %>%
  magrittr::set_colnames(1:ncol(.)) %>%
  mutate(sample=row_number()) %>%
  filter(sample %% 10 == 0) %>%
  tidyr::gather("day", "prediction", 1:(ncol(.)-1)) %>%
  mutate(day=as.numeric(day)) %>%
  left_join(., true_infections, by="day") %>%
  mutate(model="strat") %>%
  select(model, sample, prediction, true_value, day) %>%
  filter(day>fit_output$input_list$N)

ggplot(inf_samples, aes(x=day, y=prediction)) +
  geom_point()


a <- readRDS("fits/fit_sat_k=5_strat=annual_likelihood=2_cutoff=75.rds")
model_name <- "sat"
strat <- "recurrent"

plotting_func(fit_output=a)
plotting_func(fit_output=fit_weighted_com2_an)

plotting_func(fit_output=fit_weighted_sat0_an)
plotting_func(fit_output=fit_weighted_sat2_an)

plotting_func(fit_output=fit_ran0_an)
plotting_func(fit_output=fit_ran2_an)

plotting_func(fit_output=fit_lrw0_an)
plotting_func(fit_output=fit_lrw2_an)

#model_sat_strat0

input_list_sat_strat1  <- list(N=100, incidence_strat=incidence_sat_strat0$incidence_strat[1:100,], h=28, S0=incidence_sat_strat0$n_people, likelihood_type=2, k=incidence_sat_strat0$n)
fit_sat_strat1 <- sampling(model_sat, input_list_sat_strat1, iter=1000, chain=1, cores=4, control=list(adapt_delta=0.99), init=list(list(gi_mean=10, gi_sd=2, R0=2, diff=0.1, alpha=rep(0.02,3))), seed=3)
sat_strat1 <- list(fit=fit_sat_strat1, input_list=input_list_sat_strat1)
df_sat_strat1 <- cbind(data.frame(x=1:(input_list_sat_strat1$N+input_list_sat_strat1$h), 
                                  type="sm",
                                  incidence=c(rowSums(input_list_sat_strat1$incidence), rep(NA, input_list_sat_strat1$h))), 
                       matrixStats::colQuantiles(rstan::extract(fit_sat_strat1)$total_infections, probs=c(0.025, 0.5, 0.975))) %>% 
  magrittr::set_colnames(c("x", "type", "incidence", "lower", "median", "upper"))
#ggplot(df_sat_strat1, aes(x=x)) + geom_point(aes(y=incidence)) + geom_line(aes(y=median)) + geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) + theme_bw()
plotting_func(fit_output=sat_strat1)

incidence_sat <- incidence_strat_generator(simulation, range_limits=c(0), strat="annual")
input_list_sat  <- list(N=100, incidence_strat=matrix(incidence_sat$incidence_strat[1:100,], ncol=1), h=28, 
                        population=array(as.integer(incidence_sat$n_people), dim = 1), likelihood_type=2, k=incidence_sat$n)
fit_sat <- sampling(model_sat_strat0, input_list_sat, iter=1000, chain=1, cores=4, control=list(adapt_delta=0.99), init=list(list(gi_mean=10, gi_sd=2, R0=2, diff=0.1, alpha=array(as.numeric(0.02), dim = 1))), seed=3)
sat <- list(fit=fit_sat, input_list=input_list_sat)
df_sat <- cbind(data.frame(x=1:(input_list_sat$N+input_list_sat$h), 
                           type="sm",
                           incidence=c(rowSums(input_list_sat$incidence), rep(NA, input_list_sat$h))), 
                matrixStats::colQuantiles(rstan::extract(fit_sat)$total_infections, probs=c(0.025, 0.5, 0.975))) %>% 
  magrittr::set_colnames(c("x", "type", "incidence", "lower", "median", "upper"))
#ggplot(df_sat_strat1, aes(x=x)) + geom_point(aes(y=incidence)) + geom_line(aes(y=median)) + geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) + theme_bw()
plotting_func(fit_output=sat)

incidence_sat_strat0 <- incidence_strat_generator(simulation, range_limits=c(0, 16, 256), strat="annual")
input_list_sat_strat0  <- list(N=100, incidence_strat=incidence_sat_strat0$incidence_strat[1:100,], h=28, population=incidence_sat_strat0$n_people, likelihood_type=1, k=incidence_sat_strat0$n)
fit_sat_strat0 <- sampling(model_sat_strat0, input_list_sat_strat0, iter=1000, chain=1, cores=4, control=list(adapt_delta=0.99), init=list(list(gi_mean=10, gi_sd=2, R0=2, diff=0.1, alpha=rep(0.02,3))), seed=3)
sat_strat0 <- list(fit=fit_sat_strat0, input_list=input_list_sat_strat0)
df_sat_strat0 <- cbind(data.frame(x=1:(input_list_sat_strat0$N+input_list_sat_strat0$h), 
                                  type="sm",
                                  incidence=c(rowSums(input_list_sat_strat0$incidence), rep(NA, input_list_sat_strat0$h))), 
                       matrixStats::colQuantiles(rstan::extract(fit_sat_strat0)$total_infections, probs=c(0.025, 0.5, 0.975))) %>% 
  magrittr::set_colnames(c("x", "type", "incidence", "lower", "median", "upper"))
#ggplot(df_sat_strat0, aes(x=x)) + geom_point(aes(y=incidence)) + geom_line(aes(y=median)) + geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) + theme_bw()
plotting_func(fit_output=sat_strat0)




model_strat_strat <- stan_model("functions/MPox_strat_forecast_strat4.stan")

rstan::expose_stan_functions(model_strat)
rstan::expose_stan_functions(model_strat_strat)

C <- SIR_model_strat_R(a = c(1, 2, 4), 
                       beta =   0.05,
                       gamma =  15, 
                       S0 =     c(20000, 5000, 5000), 
                       E0 =     c(0,0,0), 
                       I0 =     c(2,2,2), 
                       R0 =     c(0,0,0), 
                       m =      350, 
                       k =      3, 
                       delta =  0.5) 

D <- SIR_model_strat(a = c(1, 2, 4), 
                     beta =   0.05,
                     gamma =  15, 
                     S0 =     c(20000, 5000, 5000), 
                     E0 =     c(0,0,0), 
                     I0 =     c(2,2,2), 
                     R0 =     c(0,0,0), 
                     m =      350, 
                     k =      3, 
                     delta =  0.5) 

#####
B <- SIR_model_strat_strat(a = c(1, 2, 4), 
                           beta =   0.025,
                           gamma =  45, 
                           S0 =     c(20000, 5000, 5000), 
                           E0 =     c(0,0,0), 
                           I0 =     c(2,2,2), 
                           R0 =     c(0,0,0), 
                           m =      350, 
                           k =      3, 
                           delta =  0.25) 

plot(rowSums(simplify2array(B)[,5,]))


input_n1  <- list(N=350, 
                  incidence=round(simplify2array(B)[,5,],0), 
                  h=28, 
                  S0=colSums(simplify2array(B)[1,1:4,]),  
                  R0=rep(0,3),  
                  k=3, 
                  gamma=45)

input_n2  <- list(N=350, 
                  incidence=round(rowSums(simplify2array(B)[,5,],0)), 
                  h=28, 
                  S0=colSums(simplify2array(B)[1,1:4,]),  
                  R0=rep(0,3),  
                  k=3, 
                  gamma=45)

fit_strat_n1 <- sampling(model_strat_strat, input_n1, iter=1000, chain=1, cores=4, control=list(adapt_delta=0.95), init=list(list(beta=c(0.1), phi=0.1)))
fit_strat_n2 <- sampling(model_strat,       input_n2, iter=1000, chain=1, cores=4, control=list(adapt_delta=0.95), init=list(list(beta=c(0.1), phi=0.1)))

fit_strat_n2 <- sampling(model_strat_strat, input_list_strat_n, iter=1000, chain=1, cores=4, control=list(adapt_delta=0.95), init=list(list(beta=c(0.1))))


plot(rowSums(simplify2array(B)[,5,]))
#lines(rowSums(simplify2array(B)[,5,]))
lines(rowSums(rstan::extract(fit_strat_n1)$infections[400,,]))
points(simplify2array(B)[,5,3],col="red")
lines(rstan::extract(fit_strat_n1)$infections[400,,3], col="red")
points(simplify2array(B)[,5,2],col="blue")
lines(rstan::extract(fit_strat_n1)$infections[400,,2], col="blue")
points(simplify2array(B)[,5,1],col="green")
lines(rstan::extract(fit_strat_n1)$infections[400,,1], col="green")

plot(rowSums(simplify2array(B)[,5,]))
lines((rstan::extract(fit_strat_n2)$infections[400,]))


SIR_model_strat_R <- function(a, beta, gamma, S0, E0, I0, R0, m, k, delta) {
  # Initialize a list to store SIR matrices for each stratum
  SIR <- array(0, dim = c(m, 5, k)) # 3D array to store S, E, I, R, and new infections for k strata
  tot <- numeric(k)
  prop_infected <- numeric(k)
  
  # Initial conditions
  for (i in 1:k) {
    SIR[1, 1, i] <- S0[i]
    SIR[1, 2, i] <- E0[i]
    SIR[1, 3, i] <- I0[i]
    SIR[1, 4, i] <- R0[i]
    SIR[1, 5, i] <- E0[i]
    
    tot[i] <- S0[i] + E0[i] + I0[i] + R0[i]
  }
  
  # Generate mobility matrix
  MM <- MM_generator_R(tot, a, delta)/tot*sum(tot)
  
  # Time loop
  for (t in 2:m) {
    for (i in 1:k) {
      # Calculate proportion infected for each stratum
      
      # Compute new infections
      new_infections <- min(SIR[t - 1, 1, i]/tot[i] * beta * sum(MM[i, ] * SIR[t - 1, 3, ]), SIR[t - 1, 1, i])
      new_infections
      
      # Update compartments
      SIR[t, 5, i] <- new_infections
      SIR[t, 1, i] <- SIR[t - 1, 1, i] - new_infections
      SIR[t, 2, i] <- SIR[t - 1, 2, i] + new_infections - ifelse(t > 5, SIR[t - 5, 5, i], 0) 
      SIR[t, 3, i] <- SIR[t - 1, 3, i]                  + ifelse(t > 5, SIR[t - 5, 5, i], 0) - ifelse(t==gamma, I0[i], 0) - ifelse(t > (5 + gamma), SIR[t - 5 - gamma, 5, i], 0)
      SIR[t, 4, i] <- SIR[t - 1, 4, i]                                                       + ifelse(t==gamma, I0[i], 0) + ifelse(t > (5 + gamma), SIR[t - 5 - gamma, 5, i], 0)
    }
  }
  
  return(SIR)
}

MM_generator_R <- function(p, a, delta) {
  n <- length(p)  # Number of elements in p and a
  M <- matrix(0, nrow = n, ncol = n)
  A <- matrix(0, nrow = n, ncol = n)
  denominator <- sum(p * a)  # Compute the dot product of p and a
  
  A <- diag(a * p)
  
  # Compute the matrix elements
  for (i in 1:n) {
    for (j in 1:n) {
      M[i, j] <- (p[i] * a[i] * p[j] * a[j]) / denominator
    }
  }
  
  M <- M + (A - M) * delta
  
  # Ensure M is symmetric and compute its eigenvalues
  eigenvalues <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  
  # Find the dominant eigenvalue (largest eigenvalue)
  dominant_eigenvalue <- max(eigenvalues)
  
  # Normalize the matrix M by dividing it by the dominant eigenvalue
  return(M / dominant_eigenvalue)
}

#####
A <- SIR_model_strat(a =     summary(fit_strat_strat5_annual$fit_strat_n, pars="a")$summary[,"mean"], 
                     beta =  summary(fit_strat_strat5_annual$fit_strat_n, pars="beta")$summary[,"mean"],
                     gamma = 30, 
                     S0 =    fit_strat_strat5_annual$input_list_strat_n$S0, 
                     E0 =    summary(fit_strat_strat5_annual$fit_strat_n, pars="E0")$summary[,"mean"], 
                     I0 =    summary(fit_strat_strat5_annual$fit_strat_n, pars="I0")$summary[,"mean"], 
                     R0 =    fit_strat_strat5_annual$input_list_strat_n$R0, 
                     m =     350, 
                     k =     5, 
                     delta = 1/(1+exp(summary(fit_strat_strat5_annual$fit_strat_n, pars="delta")$summary[,"mean"])))


C <- SIR_model_strat_strat(a = c(1, 2, 4), 
                           beta =   0.105,
                           gamma =  15, 
                           S0 =     c(20000, 5000, 1000), 
                           E0 =     c(0,0,0), 
                           I0 =     c(1,0,1), 
                           R0 =     c(0,0,0), 
                           m =      350, 
                           k =      3, 
                           delta =  0.5) 

D <- SIR_model_strat(a = c(1, 2, 4), 
                     beta =   0.105,
                     gamma =  15, 
                     S0 =     c(20000, 5000, 1000), 
                     E0 =     c(0,0,0), 
                     I0 =     c(1,0,1), 
                     R0 =     c(0,0,0), 
                     m =      350, 
                     k =      3, 
                     delta =  0.5) 

#A[[1]][1:3,]
B[90:100,,1]
C[[1]][90:100,]
D[[1]][90:100,]

MM_generator_R(c(20000, 5000, 5000), c(1, 2, 4), delta=0.5)
MM_generator(c(20000, 5000, 5000), c(1, 2, 4), delta=0.5)*c(20000,5000,5000)

dimnames(B) <- list(
  Time = 1:dim(B)[1],                 # Names for the time dimension
  Compartment = c("Susceptible", "Exposed", "Infectious", "Recovered", "Incidence"), # Names for compartments
  Stratum = 1:dim(B)[3]) 

model_output <- reshape2::melt(B, varnames = c("Time", "Compartment", "Stratum"), value.name = "Value")

ggplot(model_output %>% filter(Compartment=="Incidence") %>% mutate(Stratum=factor(Stratum)), aes(x=Time, color=Stratum)) +
  geom_line(aes(y=Value)) +
  theme_bw() #+
#facet_wrap(~Stratum, scales="free_y")

incidence_input <- fit_strat_strat5_annual$input_list_strat_n$incidence %>% magrittr::set_colnames(1:dim(B)[3]) %>% as.data.frame() %>%
  mutate(Time=row_number()) %>%     # Add a Time column
  tidyr::gather("Stratum", "Incidence", 1:dim(B)[3]) %>%
  mutate(Stratum=as.numeric(Stratum))

combined <- left_join(model_output, incidence_input, by=c("Time", "Stratum")) %>% mutate(Stratum=factor(Stratum))

ggplot(combined %>% filter(Compartment=="Incidence"), aes(x=Time, color=Stratum)) +
  geom_line(aes(y=Value)) +
  geom_point(aes(y=Incidence)) +
  theme_bw() +
  facet_wrap(~Stratum, scales="free_y")

fit_strat_strat5_annual$r

incidence <- A[[1]][,5] + A[[2]][,5] + A[[3]][,5] + A[[4]][,5] + A[[5]][,5]  

#####

stratified_fit_plotter <- function(input_list_strat_strat, fit_strat_strat){
  k <- input_list_strat_strat$k
  
  a <- apply(rstan::extract(fit_strat_strat)$infections, 2, colQuantiles, probs=c(0.025,0.5,0.975))
  
  df_strat_strat_group <-  as.data.frame(a) %>%
    mutate(group=rep(1:(nrow(a)/3),3)) %>%
    mutate(quantile=rep(c("lower", "median", "upper"), rep(k,3))) %>%
    tidyr::pivot_longer(cols = starts_with("V"), names_to = "day", values_to = "value") %>%
    mutate(day=as.numeric(gsub("[^0-9]", "", day))) %>%
    tidyr::pivot_wider(names_from=quantile, values_from=value)
  
  incidence_strat_df <- input_list_strat_strat$incidence %>%
    as.data.frame() %>%
    magrittr::set_colnames(1:k) %>%
    mutate(day=1:nrow(.)) %>%
    tidyr::gather("group", "incidence", 1:k) %>%
    mutate(group=as.numeric(group), day=as.numeric(day))
  
  df_strat_strat_plot <- left_join(df_strat_strat_group, incidence_strat_df, by=c("group", "day")) %>%
    mutate(group=factor(group, levels=unique(sort(group))))
  
  p <- ggplot(df_strat_strat_plot, aes(x=day)) + 
    geom_point(aes(y=incidence)) + 
    geom_line(aes(y=median)) + 
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) + 
    theme_bw() +
    facet_grid(~group, scales="free_y")
  
  #print(p)
  return(p)
}

input_list_generator <- function(simulation, n, strata="annual", strat_strat=F){
  if(strata=="annual"){
    max_degree <- max(simulation$degrees_annual)
    range_limits_n <- floor(exp(seq(log(1), log(max_degree), length.out = n+1)))[1:n]
    ranges_n <- ranges_generator(simulation, range_limits=range_limits_n)
    
    incidence_strat_n <- incidence_strat_generator(simulation, ranges=ranges_n)  
  }
  
  if(strata=="main"){
    incidence_strat_n <- list()
    incidence_strat_n$incidence_strat <- simulation$incidence_matrix_recurrent %>% unclass() %>% as.matrix() 
    n <- ncol(incidence_strat_n$incidence_strat) 
    incidence_strat_n$n_people <- setNames(as.integer(table(degree(simulation$g_main))),c(0:4))
  }  
  
  if(strat_strat==F) input_list_strat_n  <- list(N=timesteps, incidence=simulation$incidence$incidence, h=28, S0=incidence_strat_n$n_people,  E0=rep(0,n),  R0=rep(0,n),  k=n, gamma=30, delta=0)
  
  if(strat_strat==T) input_list_strat_n  <- list(N=timesteps, incidence=incidence_strat_n$incidence_strat, h=28, S0=incidence_strat_n$n_people,  E0=rep(0,n),  R0=rep(0,n),  k=n, gamma=30, delta=0)
  
  return(list(incidence_strat_n=incidence_strat_n, input_list_strat_n=input_list_strat_n))
}

fit_strat_n <- function(simulation, n, strata="annual", strat_strat=F){
  a <- input_list_generator(simulation, n, strata=strata, strat_strat=F)
  incidence_strat_n <- a$incidence_strat_n
  input_iist_strat_n <- a$input_list_strat_n
  if(strata=="main") n=5
  
  if(strat_strat==F){
    input_list_strat_n  <- list(N=timesteps, incidence=simulation$incidence$incidence, h=28, S0=incidence_strat_n$n_people,  E0=rep(0,n),  R0=rep(0,n),  k=n, gamma=20, delta=0)
    fit_strat_n <- sampling(model_strat, input_list_strat_n, iter=1000, chain=1, cores=4, control=list(adapt_delta=0.95))
    df_strat_n <- cbind(data.frame(x=1:(input_list_strat_n$N+input_list_strat_n$h), type="strat"), matrixStats::colQuantiles(rstan::extract(fit_strat_n)$infections, probs=c(0.025, 0.5, 0.975)), incidence=c(input_list_strat_n$incidence, rep(NA, input_list_strat_n$h))) %>% magrittr::set_colnames(c("x", "type", "lower", "median", "upper", "incidence"))
    p <- ggplot(df_strat_n, aes(x=x)) + geom_point(aes(y=incidence)) + geom_line(aes(y=median)) + geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) + theme_bw()
    print(p)
    
    return(list(df_strat_n=df_strat_n, fit_strat_n=fit_strat_n, p=p, input_list_strat_n=input_list_strat_n))
  } 
  
  if(strat_strat==T){
    #timesteps=200
    input_list_strat_n  <- list(N=timesteps, incidence=incidence_strat_n$incidence_strat[1:timesteps,], h=28, S0=incidence_strat_n$n_people,  R0=rep(0,n),  k=n, gamma=45, delta=0)
    
    fit_strat_n <- sampling(model_strat_strat, input_list_strat_n, iter=1000, chain=1, cores=4, control=list(adapt_delta=0.95), init=list(list(beta=c(0.01))))
    df_strat_n <- cbind(data.frame(x=1:(input_list_strat_n$N+input_list_strat_n$h), type="strat"), matrixStats::colQuantiles(apply(rstan::extract(fit_strat_n)$infections, c(1,2), sum), probs=c(0.025, 0.5, 0.975)), incidence=c(rowSums(input_list_strat_n$incidence), rep(NA, input_list_strat_n$h))) %>% magrittr::set_colnames(c("x", "type", "lower", "median", "upper", "incidence"))
    p <- ggplot(df_strat_n, aes(x=x)) + geom_point(aes(y=incidence)) + geom_line(aes(y=median)) + geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) + theme_bw()
    q <- stratified_fit_plotter(input_list_strat_strat=input_list_strat_n, fit_strat_strat=fit_strat_n)
    r <- ggpubr::ggarrange(p, q, nrow=2)
    print(r)
    
    return(list(df_strat_n=df_strat_n, fit_strat_n=fit_strat_n, p=p, q=q, r=r, input_list_strat_n=input_list_strat_n))
  } 
  
  
}

fit_strat2_annual       <- fit_strat_n(simulation, n=3,    strata="annual", strat_strat=F)
fit_strat5_annual       <- fit_strat_n(simulation, n=5,    strata="annual", strat_strat=F)
fit_strat_main          <- fit_strat_n(simulation, n=NULL, strata="main",   strat_strat=F)

fit_strat_strat2_annual <- fit_strat_n(simulation, n=3,    strata="annual", strat_strat=T)
fit_strat_strat5_annual <- fit_strat_n(simulation, n=5,    strata="annual", strat_strat=T)
fit_strat_strat8_annual <- fit_strat_n(simulation, n=8,    strata="annual", strat_strat=T)

fit_strat_strat_main    <- fit_strat_n(simulation, n=NULL, strata="main",   strat_strat=T)

a <- fit_strat_strat_main$input_list_strat_n$incidence %>% as.data.frame() %>% mutate(day=1:nrow(.)) %>% tidyr::gather("n_connections_rec", "value", 1:(ncol(.)-1)) %>% group_by(n_connections_rec) %>% mutate(roll_value=rollmean(value, 7, fill=NA)) %>%
  left_join(., data.frame(table(degree(simulation$g_main))) %>% magrittr::set_colnames(c("n_connections_rec", "total")), by="n_connections_rec") %>%
  mutate(prop_infected=roll_value/total)

ggplot(a, aes(x=day, y=roll_value, color=n_connections_rec)) + geom_line() + theme_bw()
ggplot(a, aes(x=day, y=prop_infected, color=n_connections_rec)) + geom_line() + theme_bw()

b <- simulation$generation_interval_matrix %>% mutate(week=floor(t_E/7)) %>% group_by(week) %>% summarise(n_connections_rec=mean(n_connections_rec), n=n(), n_connections_sd=sd(n_connections_rec)/sqrt(n))

ggplot(b, aes(x=week, y=n_connections_rec)) + geom_line() + theme_bw()

b2 <- simulation$generation_interval_matrix %>% mutate(week=floor(t_E/7), n_connections_rec=as.character(n_connections_rec)) %>% group_by(week) %>% select(person, week, n_connections_rec) %>% 
  group_by(week, n_connections_rec) %>% summarise(count=n()) %>% mutate(n_connections_rec=as.character(n_connections_rec)) %>%
  left_join(., data.frame(table(degree(simulation$g_main))) %>% magrittr::set_colnames(c("n_connections_rec", "total")) %>% mutate(n_connections_rec=as.character(n_connections_rec)), by="n_connections_rec") %>%
  mutate(prop=count/total, n_connections_rec=factor(n_connections_rec, levels=c(0:4)))

ggplot(b2, aes(x=week, y=prop, color=n_connections_rec)) + geom_line() + theme_bw()

c <- simulation$generation_interval_matrix %>% mutate(week=floor(t_E/7)) %>% group_by(week) %>% summarise(n_connections_annual=mean(n_connections_annual), n=n(), n_connections_sd=sd(n_connections_annual)/sqrt(n())) %>% filter(n>20)

ggplot(c, aes(x=week, y=n_connections_annual)) + geom_line() + theme_bw()

stratified_fit_plotter(input_list_strat_strat=fit_strat_strat2_annual$input_list_strat_n, fit_strat_strat=fit_strat_strat2_annual$fit_strat_n)
stratified_fit_plotter(input_list_strat_strat=fit_strat_strat2_main$input_list_strat_n, fit_strat_strat=fit_strat_strat2_main$fit_strat_n)

a <- MM_generator(p=fit_strat_strat5_annual$input_list_strat_n$S0+fit_strat_strat5_annual$input_list_strat_n$E0+fit_strat_strat5_annual$input_list_strat_n$R0+summary(fit_strat_strat5_annual$fit_strat_n, pars="I0")$summary[,"mean"],
                  a=summary(fit_strat_strat5_annual$fit_strat_n, pars="a")$summary[,"mean"], 
                  delta=1/(1+exp(summary(fit_strat_strat5_annual$fit_strat_n, pars="delta")$summary[,"mean"])))


comparison_function <- function(input_lists, N, h){
  output <-  data.frame(x=numeric(), stop=numeric(), type=character(), lower=numeric(), median=numeric(), upper=numeric(), incidence=integer())
  samples <- data.frame(sample=numeric(), day=numeric(), prediction=numeric(), true_value=numeric(), model=character(), cutoff=numeric())
  
  input_list2 <- input_list_generator(simulation, 2, "annual", F)$input_list
  input_list5 <- input_list_generator(simulation, 5, "annual", F)$input_list
  input_list_main <- input_list_generator(simulation, NULL, "main", F)$input_list
  
  input_list_strat2 <- input_list_generator(simulation, 2, "annual", T)$input_list
  input_list_strat5 <- input_list_generator(simulation, 5, "annual", T)$input_list
  input_list_strat_main <- input_list_generator(simulation, NULL, "main", T)$input_list
  
  models     <-  list(rw=model_rw,       sm=model_sm,       model_comp=model_comp, model_strat2=model_strat, model_strat5=model_strat, model_strat_main=model_strat, model_strat_strat2=model_strat_strat, model_strat_strat5=model_strat_strat, model_strat_strat_main=model_strat_strat)
  input_lists <- list(input_list_origin, input_list_origin, input_list_comp,       input_list2,              input_list5,              input_list_main,              input_list_strat2,                    input_list_strat5,                    input_list_strat_main)
  
  #for(k in 1:5){
  for(j in c(150)){
    for(i in 7){#1:length(models)){
      print(c(i, j))
      
      model <- models[[i]]
      input_list <- input_lists[[i]]
      
      input_list$N <- j
      
      if(i<7){
        input_list$incidence <- input_list$incidence[1:j]
        input_list$incidence1 <- input_list$incidence1[1:j]
      }
      if(i>=7) input_list$incidence <- input_list$incidence[1:j,]
      
      input_list$data_days <- 1:j
      
      forecast_days <- (input_list$N+1):(input_list$N+input_list$h)
      type <- names(models)[i]
      
      fit <- sampling(model, input_list, iter=1000, chain=1, cores=4)
      
      if(i>=7)  output <- rbind(output,data.frame(x=1:(j+input_list$h), stop=j, type=type, matrixStats::colQuantiles(apply(rstan::extract(fit)$infections,c(1,2),sum), probs=c(0.025, 0.5, 0.975)), incidence=rowSums(input_lists[[i]]$incidence[1:(j+input_list$h),])) %>% magrittr::set_colnames(colnames(output)))
      else      output <- rbind(output,data.frame(x=1:(j+input_list$h), stop=j, type=type, matrixStats::colQuantiles(apply(rstan::extract(fit)$infections,c(1,2),sum), probs=c(0.025, 0.5, 0.975)), incidence=input_lists[[i]]$incidence[1:(j+input_list$h)]) %>% magrittr::set_colnames(colnames(output)))
      
      joint_fit <- list(input_raw=input_lists[[i]], stan_input=input_list, fit=fit) 
      saveRDS(joint_fit, paste0("fits/output_",type,"_",j,".rds"))
      
      #if(i==2) samples_raw <- (exp(rstan::extract(joint_fit$fit)$log_incidence[,forecast_days])-1) 
      print(dim(samples_raw))
      if(i>=7) samples_raw <- apply(rstan::extract(joint_fit$fit)$infections,c(1,2),sum)[,forecast_days] 
      else     samples_raw <- rstan::extract(joint_fit$fit)$infections[,forecast_days] 
      
      true_value <- joint_fit$input_raw$incidence[forecast_days]
      
      samples_mod <- samples_raw %>%
        magrittr::set_colnames(forecast_days) %>%
        as.data.frame() %>% mutate(sample=row_number()) %>%
        tidyr::pivot_longer(cols=c(1:joint_fit$stan_input$h), names_to="day", values_to="prediction") %>%
        mutate(day=as.numeric(day)) %>%
        left_join(., data.frame(day=forecast_days, 
                                true_value=true_value), by="day") %>%
        mutate(model=type, cutoff=j)
      
      samples <- rbind(samples, samples_mod)
      
    }
  }
  output_list <- list(input=input_list, output=output, samples=samples)
  saveRDS(output_list, "fits/output_combined2.rds")
}

output_list <- readRDS("fits/output_combined.rds")

output_plot <- output_list$output %>% filter(x < stop+26) %>% mutate(stop=factor(stop))  %>%
  mutate(type=factor(type, levels=c("rw", "lrw", "sm", "model_comp", "model_strat2", "model_strat5", "model_strat_main", "model_strat_strat2", "model_strat_strat5", "model_strat_strat_main"))) 
#group_by(stop, seed, type, x) %>%
#summarise(lower=mean(lower), median=mean(median), upper=mean(upper), incidence=mean(incidence)) 

fits <- ggplot(output_plot %>% filter(), aes(x=x, fill=type)) +
  geom_line(aes(y=median, color=type)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.5) +
  geom_point(aes(y=incidence, color=type)) +
  theme_bw() +
  facet_grid(stop ~ type, scales="free") +
  geom_vline(aes(xintercept = as.numeric(as.character(stop))), linetype="dashed") +
  theme(panel.grid=element_blank()) #+
#lims(y=c(0,200))

fits
ggsave("figures/fits.png", fits)

fits2 <- ggplot(data=output_plot %>% mutate(stop=as.numeric(as.character(stop))) %>% filter(x>=stop, x<stop+28), aes(x=x, fill=type)) +
  geom_line(aes(y=median, color=type)) + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.5) +
  geom_point(data=output_plot, aes(y=incidence, color=type)) +
  theme_bw() +
  facet_wrap(~ type, scales="free", nrow=3) +
  geom_vline(aes(xintercept = as.numeric(as.character(stop))), linetype="dashed") +
  theme(panel.grid=element_blank(),
        legend.position="none") +
  labs(y="Incidence", x="Day")

fits2
ggsave("figures/fits2.png", fits2, width=10, height=3)

samples <- output_list$samples

sc_scores <- score(samples)
sc_summary <- summarise_scores(sc_scores, by=c("cutoff","model")) %>% arrange(cutoff) #%>%
#mutate(model=factor(model, levels=c("rw", "lrw", "sm", "comp", "homo", "homo_dd", "ass", "ass_dd", "ass_strat", "ass_strat_dd")))
sc_summary

fits3 <- ggplot(sc_summary, aes(x = cutoff, y = crps, fill = model)) +
  geom_col(position = "dodge") +  # "dodge" puts bars next to each other
  labs(title = "Continuous Ranked Probability Score by Cutoff and Model", x = "Cutoff", y = "CRPS") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  scale_y_log10() +
  theme()#legend.position="none")

fits3
ggsave("figures/fits3.png", fits3, width=10, height=3)

















###
model_homo <- stan_model("functions/MPox_strat_homo_forecast.stan")
k=10
input_list_comp <- list(N=timesteps, incidence=incidence, h=28, V=rep(0,148), S0=rep(sum(n_people)/2,2), R0=rep(0, 2), k=2, edges_group=rep(sum(round(n_people*n_partners))/2,2), n_partners_k=rep(sum(round(n_people*n_partners))/2,2)/rep(sum(n_people)/2,2))
input_list_homo <- list(N=timesteps, incidence=incidence, h=28, V=rep(0,148), S0=n_people, R0=rep(0, k), k=length(n_people), edges_group=round(n_people*n_partners), n_partners_k=n_partners)
fit_homo <- sampling(model_homo, input_list_homo, iter=1000, chain=1, cores=4, init=list(list(beta=0.01)))#rep(0.01, 163))))
df_homo <- cbind(data.frame(x=1:(input_list_homo$N+input_list_homo$h), type="strat2"), matrixStats::colQuantiles(rstan::extract(fit_homo)$infections, probs=c(0.025, 0.5, 0.975)), incidence=c(input_list_homo$incidence, rep(NA, input_list_homo$h))) %>% magrittr::set_colnames(c("x", "type", "lower", "median", "upper", "incidence"))
ggplot(df_homo, aes(x=x)) + geom_line(aes(y=median)) + geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) + geom_point(aes(y=incidence)) + theme_bw()#

model_ass <- stan_model("functions/MPox_strat_ass_forecast.stan")
input_list_ass <- list(N=timesteps, incidence=incidence[1:timesteps], h=28, V=rep(0,75), S0=n_people, R0=rep(0, k), k=k, edges_group=round(n_people*n_partners), n_partners_k=n_partners)
fit_ass <- sampling(model_ass, input_list_ass, iter=1000, chain=1, cores=4, init=list(list(beta=0.01)))#rep(0.01, 163))))
df_ass <- cbind(data.frame(x=1:(input_list_ass$N+input_list_ass$h), type="strat2"), matrixStats::colQuantiles(rstan::extract(fit_ass)$infections, probs=c(0.025, 0.5, 0.975)), incidence=c(input_list_ass$incidence, rep(NA, input_list_ass$h))) %>% magrittr::set_colnames(c("x", "type", "lower", "median", "upper", "incidence"))
ggplot(df_ass, aes(x=x)) + geom_line(aes(y=median)) + geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) + geom_point(aes(y=incidence)) + theme_bw()

model_homo_dd_upper <- stan_model("functions/MPox_strat_homo_d_upper_forecast.stan")
input_list_homo_dd <- list(N=timesteps, incidence=incidence[1:timesteps], h=28, V=rep(0,148), S0=n_people, R0=rep(0, k), k=k, dd_param=1.81, delta=0)
fit_homo_dd <- sampling(model_homo_dd_upper, input_list_homo_dd, iter=1000, chain=1, cores=4, init=list(list(beta=0.01)))#rep(0.01, 163))))
df_homo_dd <- cbind(data.frame(x=1:(input_list_homo_dd$N+input_list_homo_dd$h), type="strat2"), matrixStats::colQuantiles(rstan::extract(fit_homo_dd)$infections, probs=c(0.025, 0.5, 0.975)), incidence=c(input_list_homo_dd$incidence, rep(NA, input_list_homo_dd$h))) %>% magrittr::set_colnames(c("x", "type", "lower", "median", "upper", "incidence"))
ggplot(df_homo_dd, aes(x=x)) + geom_line(aes(y=median)) + geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) + geom_point(aes(y=incidence)) + theme_bw()

model_ass_dd_upper <- stan_model("functions/MPox_strat_ass_d_upper_forecast.stan")
input_list_ass_dd <- list(N=timesteps, incidence=incidence[1:timesteps], h=28, V=rep(0,148), S0=n_people, R0=rep(0, k), k=k, dd_param=1.81)
fit_ass_dd <- sampling(model_ass_dd_upper, input_list_ass_dd, iter=1000, chain=1, cores=4, init=list(list(beta=0.01)))#rep(0.01, 163))))
df_ass_dd <- cbind(data.frame(x=1:(input_list_ass_dd$N+input_list_ass_dd$h), type="strat2"), matrixStats::colQuantiles(rstan::extract(fit_ass_dd)$infections, probs=c(0.025, 0.5, 0.975)), incidence=c(input_list_ass_dd$incidence, rep(NA, input_list_ass_dd$h))) %>% magrittr::set_colnames(c("x", "type", "lower", "median", "upper", "incidence"))
ggplot(df_ass_dd, aes(x=x)) + geom_line(aes(y=median)) + geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) + geom_point(aes(y=incidence)) + theme_bw()

model_ass_strat <- stan_model("functions/MPox_strat_ass_strat_forecast.stan")
input_list_ass_strat <- list(N=timesteps, incidence=incidence_strat[1:timesteps,], h=28, V=rep(0,148), S0=n_people, R0=rep(0, k), k=k, edges_group=round(n_people*n_partners), n_partners_k=n_partners)
fit_ass_strat <- sampling(model_ass_strat, input_list_ass_strat, iter=1000, chain=1, cores=4, init=list(list(beta=0.01)))#rep(0.01, 163))))
df_ass_strat <- cbind(data.frame(x=1:(input_list_ass_strat$N+input_list_ass_strat$h), type="strat2"), matrixStats::colQuantiles(apply(rstan::extract(fit_ass_strat)$infections,c(1,2),sum), probs=c(0.025, 0.5, 0.975)), incidence=c(rowSums(input_list_ass_strat$incidence), rep(NA, input_list_ass_strat$h))) %>% magrittr::set_colnames(c("x", "type", "lower", "median", "upper", "incidence"))
ggplot(df_ass_strat, aes(x=x)) + geom_line(aes(y=median)) + geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) + geom_point(aes(y=incidence)) + theme_bw()

model_ass_strat_dd_upper <- stan_model("functions/MPox_strat_ass_strat_d_upper_forecast.stan")
input_list_ass_strat_dd <- list(N=timesteps, incidence=incidence_strat[1:timesteps,], h=28, V=rep(0,75), S0=n_people, R0=rep(0, k), k=k, dd_param=1.81)
fit_ass_strat_dd <- sampling(model_ass_strat_dd_upper, input_list_ass_strat_dd, iter=1000, chain=1, cores=4, init=list(list(beta=0.01)))#rep(0.01, 163))))
df_ass_strat_dd <- cbind(data.frame(x=1:(input_list_ass_strat_dd$N+input_list_ass_strat_dd$h), type="strat2"), matrixStats::colQuantiles(apply(rstan::extract(fit_ass_strat_dd)$infections,c(1,2),sum), probs=c(0.025, 0.5, 0.975)), incidence=c(rowSums(input_list_ass_strat_dd$incidence), rep(NA, input_list_ass_strat_dd$h))) %>% magrittr::set_colnames(c("x", "type", "lower", "median", "upper", "incidence"))
ggplot(df_ass_strat_dd, aes(x=x)) + geom_line(aes(y=median)) + geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) + geom_point(aes(y=incidence)) + theme_bw()

#####

M <- matrix(NA, nrow=k, ncol=k)
for(i in 1:nrow(M)) for(j in 1:ncol(M)) M[i,j] <- n_people[i]*n_partners[i]*n_people[j]*n_partners[j]/((n_people[1]*n_partners[1])^2) * sum(n_people*n_partners)/((sum(n_people*n_partners)/(n_people[1]*n_partners[1]))^2)/2
#M <- M/(max(eigen(M)$values))

n_connections_df <- simulation$gi_mat %>% select(person, n_connections, log_connections_decile)

pairs_df <- igraph::as_data_frame(simulation$G, what="edges") %>% left_join(., n_connections_df, by=c("from"="person")) %>%
  left_join(., n_connections_df, by=c("to"="person")) 

bucketed_ass_mat <- pairs_df %>% group_by(log_connections_decile.x, log_connections_decile.y) %>%
  summarise(count=n(), .groups="drop") 

M_ass <- matrix(bucketed_ass_mat$count, nrow=k, ncol=k)
M_ass

ggplot(bucketed_ass_mat, aes(x=log_connections_decile.x, y=log_connections_decile.y, fill=count)) +
  geom_tile() +
  theme_bw()

assortativity_degree(simulation$G)
cor(pairs_df$log_connections_decile.x, pairs_df$log_connections_decile.y)

ggplot(pairs_df, aes(x=n_connections.x, y=n_connections.y)) +
  geom_point()

assortative_mixing_function <- function(mixing_matrix, delta){
  output_matrix <- matrix(nrow=nrow(mixing_matrix), ncol=ncol(mixing_matrix))
  k <- dim(mixing_matrix)[1]
  for(i in 1:nrow(mixing_matrix)) for(j in 1:ncol(mixing_matrix)) output_matrix[i,j] <- ifelse(i==j, mixing_matrix[i,j]+delta*(k-1)*mixing_matrix[i,j], mixing_matrix[i,j]-delta*mixing_matrix[i,j])
  return(output_matrix)
}

trials <- seq(0, 0.5, 0.01)
storage <- vector(length=51)

for(i in 1:51){
  delta=trials[i]
  print(delta)
  M2 <- assortative_mixing_function(M, delta)
  storage[i] <- sqrt(sum((M2-M_ass)^2))
}

M2 <- assortative_mixing_function(M, trials[which(storage==min(storage))])

M3 <- assortative_mixing_function(M, delta=0.5)

input_list_strat_ass <- list(N=timesteps, incidence=incidence, h=28, V=rep(0,148), S0=n_people, R0=rep(0, k), mu=n_partners/100, k=k, assortativity=M)
input_list_strat_ass2 <- input_list_strat_ass
input_list_strat_ass2$assortativity <- assortative_mixing_function(M, trials[which(storage==min(storage))])
input_list_strat_ass3 <- input_list_strat_ass
input_list_strat_ass3$assortativity <- assortative_mixing_function(M, delta=0.5)
input_list_strat_ass_strat <- list(N=timesteps, incidence=incidence_strat, h=28, V=rep(0,148), S0=n_people, R0=rep(0, k), mu=n_partners/100, k=k, assortativity=M)

expose_stan_functions(model_strat_ass_strat_forecast)
SIR_strat <- SIR_model


expose_stan_functions(model_strat_ass_forecast)
SIR_ass <- SIR_model

SIR_ass_R <- function(beta, gamma, S0, I0, R0, m, k, assortativity) {
  # Create an array to store S, I, R compartments
  SIR <- array(0, dim = c(m, 3, k))  # 3 compartments (S, I, R) for each of k populations
  
  tot <- numeric(k)  # Vector to store total population for each group
  
  # Initial conditions
  for (i in 1:k) {
    SIR[1, 1, i] <- S0[i]  # Initial susceptible
    SIR[1, 2, i] <- I0[i]  # Initial infected
    SIR[1, 3, i] <- R0[i]  # Initial recovered
    
    tot[i] <- S0[i] + I0[i] + R0[i]  # Total population per group
  }
  
  # Time loop (start from 2 because 1st time step is initial condition)
  for (t in 2:m) {
    for (i in 1:k) {
      effective_contact <- sum(assortativity[i,]*SIR[t-1,2,]/tot)
      
      # Update compartments for S, I, R
      SIR[t, 1, i] <- SIR[t - 1, 1, i] - beta[t] * SIR[t - 1, 1, i] * effective_contact / tot[i]  # Susceptible
      SIR[t, 2, i] <- SIR[t - 1, 2, i] + beta[t] * SIR[t - 1, 1, i] * effective_contact / tot[i] - gamma * SIR[t - 1, 2, i]  # Infected
      SIR[t, 3, i] <- SIR[t - 1, 3, i] + gamma * SIR[t - 1, 2, i]  # Recovered
    }
  }
  
  return(SIR)  # Return array with all compartments
}


ass_I <- SIR_ass(beta=rep(400, 100), gamma=0, S0=c(5000,2500,1000), I0=c(1,2,3), R0=rep(0,3), 
                 m=100, k=3, assortativity=matrix(1, ncol=3, nrow=3))

ass_R <- SIR_ass_R(beta=rep(400, 100), gamma=0, S0=c(5000,2500,1000), I0=c(1,2,3), R0=rep(0,3), 
                   m=100, k=3, assortativity=matrix(1, ncol=3, nrow=3))

strat_I <- SIR_strat(beta=rep(0.08, 100), gamma=0, S0=c(5000,2500,1000), I0=c(1,2,3), R0=rep(0,3), 
                     m=100, k=3, mu=c(1,2,5))

strat_R <- SIR_strat_R(beta=rep(0.08, 100), gamma=0, S0=c(5000,2500,1000), I0=c(1,2,3), R0=rep(0,3), 
                       m=100, k=3, mu=c(1,2,5))

plot(ass_R[,,1][,2], col="red")
lines(ass_I[[1]][,2], col="red")
points(strat_R[,2,1], col="blue")
lines(strat_I[[1]][,2], col="blue")


fit_strat_strat2 <- sampling(model_strat_strat, input_list_strat_strat2, iter=1000, chain=1, cores=4, control=list(adapt_delta=0.99), init=list(list(beta=c(0.1), E0=rep(1,2))))
df_strat_strat2 <- cbind(data.frame(x=1:(input_list_strat_strat2$N+input_list_strat_strat2$h), type="strat"), matrixStats::colQuantiles(apply(rstan::extract(fit_strat_strat2)$infections, c(1,2), sum), probs=c(0.025, 0.5, 0.975)), incidence=c(rowSums(input_list_strat_strat2$incidence), rep(NA, input_list_strat_strat2$h))) %>% magrittr::set_colnames(c("x", "type", "lower", "median", "upper", "incidence"))
ggplot(df_strat_strat2, aes(x=x)) + geom_point(aes(y=incidence)) + geom_line(aes(y=median)) + geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) + theme_bw()

fit_strat_strat3 <- sampling(model_strat_strat, input_list_strat_strat3, iter=1000, chain=1, cores=4, control=list(adapt_delta=0.99), init=list(list(beta=c(0.1), E0=rep(1,3))))
df_strat_strat3 <- cbind(data.frame(x=1:(input_list_strat_strat3$N+input_list_strat_strat3$h), type="strat"), matrixStats::colQuantiles(apply(rstan::extract(fit_strat_strat3)$infections, c(1,2), sum), probs=c(0.025, 0.5, 0.975)), incidence=c(rowSums(input_list_strat_strat3$incidence), rep(NA, input_list_strat_strat3$h))) %>% magrittr::set_colnames(c("x", "type", "lower", "median", "upper", "incidence"))
ggplot(df_strat_strat3, aes(x=x)) + geom_point(aes(y=incidence)) + geom_line(aes(y=median)) + geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) + theme_bw()

fit_strat_strat5 <- sampling(model_strat_strat, input_list_strat_strat5, iter=1000, chain=1, cores=4, control=list(adapt_delta=0.99), init=list(list(beta=c(0.1))))
df_strat_strat5 <- cbind(data.frame(x=1:(input_list_strat_strat5$N+input_list_strat_strat5$h), type="strat"), matrixStats::colQuantiles(apply(rstan::extract(fit_strat_strat5)$infections, c(1,2), sum), probs=c(0.025, 0.5, 0.975)), incidence=c(rowSums(input_list_strat_strat5$incidence), rep(NA, input_list_strat_strat5$h))) %>% magrittr::set_colnames(c("x", "type", "lower", "median", "upper", "incidence"))
ggplot(df_strat_strat5, aes(x=x)) + geom_point(aes(y=incidence)) + geom_line(aes(y=median)) + geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) + theme_bw()

fit_strat_strat10 <- sampling(model_strat_strat, input_list_strat_strat10, iter=1000, chain=1, cores=4, control=list(adapt_delta=0.99), init=list(list(beta=c(0.1))))
df_strat_strat10 <- cbind(data.frame(x=1:(input_list_strat_strat10$N+input_list_strat_strat10$h), type="strat"), matrixStats::colQuantiles(apply(rstan::extract(fit_strat_strat10)$infections, c(1,2), sum), probs=c(0.025, 0.5, 0.975)), incidence=c(rowSums(input_list_strat_strat10$incidence), rep(NA, input_list_strat_strat10$h))) %>% magrittr::set_colnames(c("x", "type", "lower", "median", "upper", "incidence"))
ggplot(df_strat_strat10, aes(x=x)) + geom_point(aes(y=incidence)) + geom_line(aes(y=median)) + geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) + theme_bw()

a <- SIR_model_strat(a=summary(fit_strat2, pars="a")$summary[,"mean"], 
                     beta=summary(fit_strat2, pars="beta")$summary[,"mean"],
                     gamma=1/30, 
                     S0=input_list_strat2$S0, 
                     E0=input_list_strat2$E0, 
                     I0=summary(fit_strat2, pars="I0")$summary[,"mean"], 
                     R0=input_list_strat2$R0, 
                     m=350, k=2, 
                     delta=1/(1+exp(summary(fit_strat2, pars="delta")$summary[,"mean"])))



rstan::expose_stan_functions(model_strat_strat)

model_sat <- stan_model("functions/MPox_saturating.stan")
incidence_sat <- incidence_strat_generator(simulation, range_limits=c(0), strat="annual")
input_list_sat  <- list(N=100, incidence1=simulation$incidence$incidence[1:100], h=28, population=incidence_sat$n_people)
fit_sat <- sampling(model_sat, input_list_sat, iter=1000, chain=1, cores=4, control=list(adapt_delta=0.99), init=list(list(gi_mean=10, gi_sd=2, R0=2, alpha=0.02)), seed=3)
df_sat <- cbind(data.frame(x=1:(input_list_sat$N+input_list_sat$h), type="sm"), matrixStats::colQuantiles(rstan::extract(fit_sat)$incidence, probs=c(0.025, 0.5, 0.975)), incidence=c(input_list_sat$incidence, rep(NA, input_list_sat$h))) %>% magrittr::set_colnames(c("x", "type", "lower", "median", "upper", "incidence"))
ggplot(df_sat, aes(x=x)) + geom_point(aes(y=incidence)) + geom_line(aes(y=median)) + geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) + theme_bw()


model_sat_strat <- stan_model("functions/MPox_saturating_strata.stan")
incidence_sat_strat <- incidence_strat_generator(simulation, range_limits=c(0, 4, 16, 64, 256), strat="annual")
input_list_sat_strat  <- list(N=100, incidence=incidence_sat_strat$incidence_strat[1:100,], h=28, population=incidence_sat_strat$n_people, K=incidence_sat_strat$n)
fit_sat_strat <- sampling(model_sat_strat, input_list_sat_strat, iter=1000, chain=1, cores=4, control=list(adapt_delta=0.99), init=list(list(gi_mean=10, gi_sd=2, R0=rep(2,5), alpha=rep(0.02,5))), seed=3)
sat_strat <- list(fit=fit_sat_strat, input_list=input_list_sat_strat)
#df_sat_strat <- cbind(data.frame(x=1:(input_list_sat$N+input_list_sat$h), type="sm"), matrixStats::colQuantiles(rstan::extract(fit_sat)$incidence, probs=c(0.025, 0.5, 0.975)), incidence=c(input_list_sat$incidence, rep(NA, input_list_sat$h))) %>% magrittr::set_colnames(c("x", "type", "lower", "median", "upper", "incidence"))
#ggplot(df_sat, aes(x=x)) + geom_point(aes(y=incidence)) + geom_line(aes(y=median)) + geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) + theme_bw()
plotting_func(fit_output=sat_strat)


time_vec <- vector(length=5)
divs <- vector(length=5)

trial_deltas <- c(0.8, 0.85, 0.9, 0.95, 0.99)

for(ad in trial_deltas){
  p <- proc.time()
  i <- which(trial_deltas==ad)
  fit1 <- sampling(R_random3, input_list, iter=iter, chains=4, cores=8, seed=(3), init=init, control = list(adapt_delta = ad)) 
  time_vec[i] <- (proc.time() - p)["elapsed"]
  divs[i] <- sum(sapply(rstan::get_sampler_params(fit1, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))
}
print(time_vec)
print(divs)

time_vec2 <- vector(length=5)
divs2 <- vector(length=5)

trial_deltas <- c(0.8, 0.85, 0.9, 0.95, 0.99)

for(ad in trial_deltas){
  p <- proc.time()
  i <- which(trial_deltas==ad)
  fit1 <- sampling(R_random4, input_list, iter=iter, chains=4, cores=8, seed=(3), init=init, control = list(adapt_delta = ad)) 
  time_vec2[i] <- (proc.time() - p)["elapsed"]
  divs2[i] <- sum(sapply(rstan::get_sampler_params(fit1, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))
  print(time_vec2)
  print(divs2)
}

time_vec3 <- vector(length=5)
divs3 <- vector(length=5)

trial_deltas <- c(0.8, 0.85, 0.9, 0.95, 0.99)

for(ad in trial_deltas){
  p <- proc.time()
  i <- which(trial_deltas==ad)
  fit1 <- sampling(R_random5, input_list, iter=iter, chains=4, cores=8, seed=(3), init=init, control = list(adapt_delta = ad)) 
  time_vec3[i] <- (proc.time() - p)["elapsed"]
  divs3[i] <- sum(sapply(rstan::get_sampler_params(fit1, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))
  print(time_vec3)
  print(divs3)
}

time_vec4 <- vector(length=5)
divs4 <- vector(length=5)

trial_deltas <- c(0.8, 0.85, 0.9, 0.95, 0.99)

for(ad in trial_deltas){
  p <- proc.time()
  i <- which(trial_deltas==ad)
  fit1 <- sampling(R_random6, input_list, iter=iter, chains=4, cores=8, seed=(3), init=init, control = list(adapt_delta = ad)) 
  time_vec4[i] <- (proc.time() - p)["elapsed"]
  divs4[i] <- sum(sapply(rstan::get_sampler_params(fit1, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))
  print(time_vec4)
  print(divs4)
}

###
###
###
###


