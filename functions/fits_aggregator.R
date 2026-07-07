library(stringr)

files <- c(list.files("fits/sim/MSM_like/", full.names=T), list.files("fits/sim/homo/", full.names=T), list.files("fits/sim/hetero/", full.names=T), list.files("fits/cities/", full.names=T))

files2 <- files[grepl("sigmoidP9", files) | grepl("expI4", files) | grepl("lrwP", files) | grepl("R_random5", files) | grepl("R_biased5", files)]
model     <- sub(".*fit_(.*)_k=.*", "\\1", files2)
epi_phase <- as.numeric(sub(".*epi_phase=(\\d+)_.*", "\\1", files2))
tags       <- idx <- str_extract(files2, "(?<=_)[0-9]+(?=T?\\.rds$)") |> as.numeric()
network <- sub(".*sim/([^/]+)/fit.*", "\\1", files2)

table(model)
table(epi_phase)
table(tags)
table(network)

for(i in files2){
  print(which(files2==i), quote=F)
  fit <- readRDS(i)
  network <- ifelse(grepl("sim", i), sub(".*sim/([^/]+)/fit.*", "\\1", i), sub(".*_([^_]+)\\.rds$", "\\1", i))
  tag <- ifelse(grepl("sim", i), as.numeric(gsub("T$", "", fit$tag)), 0)
  
  Rhat_df_raw <- data.frame(file=i,
                       model=fit$model_name,
                       tag=tag, 
                       epi_phase=fit$input_list$epi_phase,
                       cutoff=fit$cutoff,
                       max_Rhat=fit$max_Rhat,
                       iter=nrow(fit$forecast),
                       divergences=fit$total_divergences,
                       network=network)
  
  idx <- round(seq(1, nrow(fit$infections), length.out = 100))
  
  plotter_raw <- data.frame(day=1:(fit$input_list$N+fit$input_list$h),
             true_infections=fit$incidence_full) %>%
    mutate(lower=colQuantiles(fit$infections, prob=c(0.025)), 
           median=colQuantiles(fit$infections, prob=c(0.5)),
           upper=colQuantiles(fit$infections, prob=c(0.975))) %>%
    mutate(model=fit$model_name, tag=tag, epi_phase=fit$input_list$epi_phase, cutoff=fit$cutoff, network=network)
  
  main_params_raw <- as.data.frame(fit$summary) %>%
    mutate(model=fit$model_name, tag=tag, epi_phase=fit$input_list$epi_phase, cutoff=fit$cutoff, network=network)
    
  CRPS_df_raw <- as.data.frame(t(fit$infections[idx,])) %>%
    mutate(day=1:(fit$input_list$N+fit$input_list$h),
           true_infections=fit$incidence_full) %>% 
    filter(day > fit$cutoff) %>% 
    mutate(f_week=floor((day-fit$cutoff-1)/7)) %>%
    tidyr::gather("sample", "value", 1:(ncol(.)-3)) %>%
    mutate(sample=as.numeric(sub("^V", "", sample))) %>%
    group_by(sample, f_week) %>%
    summarise(true_infections=sum(true_infections), value=sum(value), .groups="drop") %>%
    mutate(model=fit$model_name, tag=tag, epi_phase=fit$input_list$epi_phase, cutoff=fit$cutoff, network=network) 
  
  if(i == files2[1]){
    Rhat_df <- Rhat_df_raw
    CRPS_df <- CRPS_df_raw
    plotter <- plotter_raw
    main_params <- main_params_raw
  } 
  else{
    Rhat_df <- rbind(Rhat_df, Rhat_df_raw)
    CRPS_df <- rbind(CRPS_df, CRPS_df_raw)
    plotter <- rbind(plotter, plotter_raw)
    main_params <- rbind(main_params, main_params_raw)
  } 
}

df <- readRDS("new_outputs_hh.rds")
saveRDS(object=list(Rhat_df=Rhat_df, CRPS_df=CRPS_df, plotter=plotter, main_params=main_params), file="fit_cit2.rds")


ggplot(plotter %>% filter(network=="homo", epi_phase==10, model=="sigmoidP4"), aes(x=day)) +
  geom_point(aes(y=true_infections)) +
  geom_line(aes(y=median)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
  theme_bw() +
  facet_wrap(~tag, scales="free") +
  theme(panel.grid = element_blank())

out <- readRDS("new_outputs2.rds")
df_filtered <- Rhat_df %>% filter(divergences >0 | max_Rhat > 1.05 | is.nan(max_Rhat)) %>% arrange(desc(max_Rhat)) %>% filter(!network %in% c("MSM_like", "homo", "hetero"))


for(k in 3:4){
  print(df_filtered[k,])
  fit <- readRDS(df_filtered$file[k])
  tag <- readr::parse_number(fit$tag)
  
  adapt_delta <- 0.9
  iter <- 3000
  
  if(grepl("cities", df_filtered$file[k])){
    simulation <- gsub("[^a-zA-Z ]", "", fit$tag)
    
    cum_inc_raw <- cumsum(dat[,j])
    first_nonzero <- which(cum_inc_raw != 0)[1]
    cum_inc <- c(rep(0,7), cum_inc_raw[first_nonzero:length(cum_inc_raw)])
    
    tot_inc <- max(cum_inc)
    lower <- which(cum_inc > tot_inc * 0.05)[1]
    upper <- which(cum_inc > tot_inc * 0.95)[1]-28
    times <- round(seq(lower, upper, length.out=10))
    
    tag <- fit$tag
    cutoff <- times[fit$epi_phase]
    
    a0 <- fit_func2(simulation, model_name=p, cutoff=cutoff, iter=2000, epi_phase=which(times==i), dir="CIT", tag=tag, adapt_delta=0.98)
    
  }
  
  else{
    dir <- sub(".*sim/([^/]+)/fit.*", "\\1", df_filtered$file[k])
    
    if(dir != "MSM_like") simulation <- readRDS(paste0("simulation2/",dir,"/epidemic_", dir, "_", ifelse(dir=="homo", 30+as.numeric(tag), as.numeric(substr(tag, start=1, stop=nchar(tag)-1))),".rds"))
    if(dir == "MSM_like") simulation <- readRDS(paste0("simulation2/",dir,"/epidemic_", ifelse(dir=="homo", 30+as.numeric(tag), as.numeric(substr(tag, start=1, stop=nchar(tag)-1))),".rds"))
    
    inc <- simulation$outbreak$incidence_vector
    cum_inc <- cumsum(inc)
    tot_inc <- max(cum_inc)
    lower <- which(cum_inc > tot_inc * 0.05)[1]
    
    peak_inc <- max(inc)
    peak_inc_t <- which(inc==max(inc))
    
    
    if(df_filtered$network[k] != "hetero") upper <- which(cum_inc > tot_inc * 0.95)[1]-28
    if(df_filtered$network[k] == "hetero") upper <- which(inc<(peak_inc*0.03) & (1:length(inc)) > peak_inc_t)[1]
    
    times <- round(seq(lower, upper, length.out=10))
    
    cutoff <- times[fit$epi_phase]
    
    #if(fit$total_divergences == 0 & fit$max_Rhat < 1.05 & !is.nan(fit$max_Rhat) ) next
    
    if(fit$total_divergences > 0) adapt_delta <- 0.995
    if(fit$max_Rhat > 1.05 | is.nan(fit$max_Rhat)) iter <- 5000
    
    a5 <- fit_func2(simulation=fit$simulation, model_name=fit$model_name, 
                    cutoff=cutoff, epi_phase=fit$epi_phase, 
                    dir=dir, iter=iter, chains=4, tag=fit$tag, adapt_delta=adapt_delta, keep_fit=F)
    
    plotting_func(a5)
  }  
  
}
