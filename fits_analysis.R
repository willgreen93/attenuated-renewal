fit_files <- list.files("fits/sim/fit_sim4/", full.names=T)

for(i in fit_files){
  
  fit <- readRDS(i)
  
  print(fit$model_name, quote=F)
  
  plotting_func(fit)
  
  df_raw <- data.frame(file=i, Rhat=fit$max_Rhat, divergences=fit$total_divergences, 
                   simul=fit$tag, epi_phase=fit$input_list$epi_phase)
  
  if(i==fit_files[1]) df <- df_raw 
  else df <- rbind(df, df_raw)
}

i <- fit_files[1]

df %>% filter(divergences > 1)
