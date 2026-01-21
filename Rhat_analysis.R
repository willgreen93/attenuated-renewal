library(stringr)

file_list <- list()
R_hat_list <- list()
model_name_list <- list()
divergences_list <- list()

df <- data.frame(file_name=character(), Rhat=numeric(), divergences=numeric())
j <- 1

for(i in files_list$file_name){#c(list.files("fits/sim/fit_sim3", full.names = TRUE), list.files("fits/sim/fit_homo", full.names = TRUE), list.files("fits/sim/fit_het", full.names = TRUE), list.files("fits/cities", full.names = TRUE))){
  
  file <- readRDS(i)
  
  divs <- sapply(rstan::get_sampler_params(file$fit, inc_warmup = FALSE), function(x) sum(x[, "divergent__"]))
  if(sum(divs)>0) print(divs)
  
  df[j,] <- c(i, file$max_Rhat, file$total_divergences) 
  print(df[j,], quote=F)
  
  j <- j+1
} 

files <- unlist(file_list)
Rhats <- unlist(R_hat_list)
models <- unlist(model_name_list)
divergences <- unlist(divergences_list)

df_RW5 <- data.frame(files=files, Rhats=Rhats, models=models, divergences=divergences) %>%
  filter(!models %in% c("R_biased", "R_random", "sigmoid_inf", "lrw2", "expdelay2")) %>%
  mutate(source = case_when(str_detect(files, "/cities/") ~ "cities", str_detect(files, "/sim/") ~ str_match(files, "/sim/([^/]+)/")[, 2], TRUE ~ NA_character_))

df_RW5

saveRDS(df, "fits/outputs/Rhats_df2.rds")
df <- readRDS("fits/outputs/Rhats_df.rds")

df %>% group_by(models, source) %>%
  summarise(n_div_0 = sum(divergences <= 3,  na.rm = TRUE), n_div_gt0 = sum(divergences > 3,  na.rm = TRUE), .groups = "drop") %>%
  mutate(pct_div = n_div_gt0/(n_div_0+n_div_gt0)) %>%
  select(-c(n_div_0,n_div_gt0)) %>%
  pivot_wider(names_from = source, values_from = pct_div)

df %>% group_by(models, source) %>%
  summarise(n=n(), n_div_gt0 = sum(divergences > 3,  na.rm = TRUE), .groups = "drop") %>%
  mutate(pct=n_div_gt0/n) %>%
  print(n=32)

df %>% group_by(models, source) %>%
  summarise(n=n(), n_div_gt0 = sum(divergences > 50,  na.rm = TRUE), .groups = "drop") %>%
  mutate(pct=n_div_gt0/n) %>%
  print(n=32)

df %>%
  filter(!models %in% c("R_biased", "R_random", "sigmoid_inf", "lrw2", "expdelay2")) %>%
  summarise(n=n(), n_div_gt0 = sum(divergences > 3,  na.rm = TRUE), .groups = "drop") %>%
  mutate(pct=n_div_gt0/n) 

df %>%
  filter(!models %in% c("R_biased", "R_random", "sigmoid_inf", "lrw2", "expdelay2")) %>%
  summarise(n=n(), n_div_gt0 = sum(divergences > 50,  na.rm = TRUE), .groups = "drop") %>%
  mutate(pct=n_div_gt0/n) 

files_list <- df2[df2$divergences >= 1,]

for(j in 1:nrow(files_list)){
  filename <- files_list$file_name[j]
  print(filename)
  a <- readRDS(filename)
  
  if (grepl("cities", filename)) dir <- ""  # No need to extract anything
  else dir <- sub(".*sim/([^/]+)/.*", "\\1", filename)
  print(dir)
  
  simulation <- a$simulation
  if(length(simulation) < 10) simulation <- a$tag
  model_name <- a$model_name
  cutoff <- a$cutoff
  epi_phase <- a$input_list$epi_phase
  
  cum_inc <- cumsum(a$input_list$incidence_strat)

  tot_inc <- max(cum_inc)
  lower <- which(cum_inc > tot_inc * 0.05)[1]
  upper <- which(cum_inc > tot_inc * 0.95)[1]-28
  times <- round(seq(lower, upper, length.out=10))
  
  tag <- a$tag
  #new_model_name <- paste0(model_name,"P")
  
  print(c("Divergences before = ", files_list$divergences[j]), quote=F)
  
  a12 <- fit_func(simulation, model_name=model_name, range_limits=c(0), strat="annual", cutoff=cutoff, tag=tag, keep_fit=TRUE, epi_phase=epi_phase, decay_type=NA, dir=dir, iter=1000, chains=4, adapt_delta=0.8)
  print(sapply(rstan::get_sampler_params(a12$fit, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))
  cat("\n")
  
}

for(j in df$files){
  print(which(df$files[df$divergences>0]==j))
  
  source <- case_when(str_detect(j, "/cities/") ~ "cities", str_detect(j, "/sim/") ~ str_match(j, "/sim/([^/]+)/")[, 2], TRUE ~ NA_character_)
  
  a <- readRDS(j)
  
  divs <- sapply(rstan::get_sampler_params(a$fit, inc_warmup = FALSE), function(x) sum(x[, "divergent__"]))
  if(sum(divs<=3)) next
  
  model_name <- a$model_name
  cutoff <- a$cutoff
  epi_phase <- a$input_list$epi_phase
  
  if(source=="cities"){
    simulation <- a$tag
    cum_inc <- cumsum(a$input_list$incidence_strat)
    dir <- ""
  }   
  if(source != "cities"){
    simulation <- a$simulation
    cum_inc <- cumsum(simulation$incidence$incidence)
    dir <- source
  } 
  
  tot_inc <- max(cum_inc)
  lower <- which(cum_inc > tot_inc * 0.05)[1]
  upper <- which(cum_inc > tot_inc * 0.95)[1]-28
  times <- round(seq(lower, upper, length.out=10))
  
  tag <- a$tag
  iter <- a$fit@sim$iter
  
  print(paste0("Source = ", source, "     Model = ", model_name, "     Iterations = ", iter, "     Divergent iterations = ", sum(divs)), quote=F)
  
  a12 <- fit_func(simulation, model_name="R_random5", range_limits=c(0), strat="annual", cutoff=cutoff, tag=tag, keep_fit=TRUE, epi_phase=epi_phase, decay_type=NA, dir=dir, iter=2000, chains=4, adapt_delta=0.99)
  a13 <- fit_func(simulation, model_name="R_random5", range_limits=c(0), strat="annual", cutoff=cutoff, tag=tag, keep_fit=TRUE, epi_phase=epi_phase, decay_type=NA, dir=dir, iter=2000, chains=4, adapt_delta=0.99)
  a14 <- fit_func(simulation, model_name="R_random5", range_limits=c(0), strat="annual", cutoff=cutoff, tag=tag, keep_fit=TRUE, epi_phase=epi_phase, decay_type=NA, dir=dir, iter=2000, chains=4, adapt_delta=0.99)
  a15 <- fit_func(simulation, model_name="R_random5", range_limits=c(0), strat="annual", cutoff=cutoff, tag=tag, keep_fit=TRUE, epi_phase=epi_phase, decay_type=NA, dir=dir, iter=2000, chains=4, adapt_delta=0.99)
  a16 <- fit_func(simulation, model_name="R_random5", range_limits=c(0), strat="annual", cutoff=cutoff, tag=tag, keep_fit=TRUE, epi_phase=epi_phase, decay_type=NA, dir=dir, iter=2000, chains=4, adapt_delta=0.99)
  #a13 <- fit_func(simulation, model_name="R_random6", range_limits=c(0), strat="annual", cutoff=cutoff, tag=tag, keep_fit=TRUE, epi_phase=epi_phase, decay_type=NA, dir=dir, iter=iter, chains=4, adapt_delta=0.99)
  
  print(sapply(rstan::get_sampler_params(a12$fit, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))
}

df <- data.frame(simul=numeric(), R0_m=numeric(), R0_2sd=numeric(), rw_sigma_m=numeric(), rw_sigma_2sd=numeric(), divergences = numeric())

simuls <- 1:50

for(i in simuls){
  print(i)
  fit <- readRDS(paste0("fits/sim/fit_sim2/fit_R_random5_k=1_strat=annual_likelihood=1_epi_phase=1_",i,".rds"))$fit
  
  divergences <- sum(sapply(rstan::get_sampler_params(fit, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))
  print(divergences)
  df[i,1] <- i
  df[i, c(2,3)] <- round(summary(fit)$summary["R0[1]",c("50%", "97.5%")],2)
  df[i, c(4,5)] <- round(summary(fit)$summary["rw_sigma",c("50%", "97.5%")],2)
  df[i, 6] <- divergences
} 

df2 <- data.frame(simul=numeric(), R0_m=numeric(), R0_2sd=numeric(), rw_sigma_m=numeric(), rw_sigma_2sd=numeric(), divergences = numeric())

for(i in simuls){
  print(i)
  fit <- readRDS(paste0("fits/sim/fit_sim2/fit_R_random5_k=1_strat=annual_likelihood=1_epi_phase=10_",i,".rds"))$fit
  
  divergences <- sum(sapply(rstan::get_sampler_params(fit, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))
  print(divergences)
  df2[i,1] <- i
  df2[i, c(2,3)] <- round(summary(fit)$summary["R0[1]",c("50%", "97.5%")],2)
  df2[i, c(4,5)] <- round(summary(fit)$summary["rw_sigma",c("50%", "97.5%")],2)
  df2[i, 6] <- divergences
} 

exp(mean(log(df2$R0_m)))
sd(log(df2$R0_m))

ggplot(df %>% tidyr::gather(param, value, 2:5), aes(x=simul, y=value, color=divergences>10)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~param, scales="free_y")

new_Rhats <- list()
new_files <- list()

for(j in files[Rhats>1.05 & !models %in% c("R_random", "R_biased")]){
  i <- which(files[Rhats>1.05 & !models %in% c("R_random", "R_biased")] == j)
  
  a <- readRDS(j)
  new_Rhats[[i]] <- (summary(a$fit)$summary[,"Rhat"] %>% max())
  new_files[[i]] <- j
  print(j)
  print(new_Rhats[[i]])
  
}

unlist(new_Rhats)

mod_files <- list.files("fits/sim/fit_sim2/", full.names=T, pattern = "powdelay")
base <- readRDS(mod_files[1])$fit@stanmodel@model_code

k <- 1

for(i in mod_files){
  print(k)
  file <- readRDS(i)
  if(file$fit@stanmodel@model_code==base) print(TRUE)
  else print(i)
  k <- k+1
} 

###
###
###
###
p <- proc.time()
a1 <- fit_func(simulation, model_name="R_random3", range_limits=c(0), strat="annual", cutoff=cutoff, tag=tag, keep_fit=TRUE, epi_phase=epi_phase, decay_type=NA, dir="old", iter=iter, chains=4, adapt_delta=0.99)
print((proc.time()-p)["elapsed"])

a2 <- fit_func(simulation, model_name="R_random4", range_limits=c(0), strat="annual", cutoff=cutoff, tag=tag, keep_fit=TRUE, epi_phase=epi_phase, decay_type=NA, dir="old", iter=iter, chains=4, adapt_delta=0.99)

p <- proc.time()
a3 <- fit_func(simulation, model_name="R_random5", range_limits=c(0), strat="annual", cutoff=cutoff, tag=tag, keep_fit=TRUE, epi_phase=epi_phase, decay_type=NA, dir="old", iter=iter, chains=4, adapt_delta=0.99)
print((proc.time()-p)["elapsed"])
###
###
###
###

