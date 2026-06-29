outbreak_files <-  c(list.files("simulation2/MSM_like/", full.names=T)) # c(list.files("simulation2/homo/", full.names=T), list.files("simulation2/hetero/", full.names=T))

for(i in outbreak_files){
  simul <- readRDS(i)
  inc <- simul$outbreak$incidence_vector
  cum_inc <- cumsum(inc)
  tot_inc <- sum(inc)
  peak_inc <- max(inc)
  peak_inc_t <- which(inc==max(inc))
  first_lower_5 <- which(inc<(peak_inc*0.03) & (1:length(inc)) > peak_inc_t)[1]
  
  df_incidence_raw <- data.frame(file=i, n=simul$graph$graph_params$n, tag2=simul$graph$tag, tot_inf = sum(simul$outbreak$incidence_vector), 
                                 time=1:length(simul$outbreak$incidence_vector), incidence=simul$outbreak$incidence_vector,
                                 dd_param=simul$graph$graph_params$alpha, dd_upper=simul$graph$graph_params$kmax,
                                 start=which(cum_inc > (tot_inc*0.05))[1],
                                 end=which(cum_inc > (tot_inc*0.95))[1]-28,
                                 end2 = first_lower_5)
  
  if(i==outbreak_files[1]) df_incidence <- df_incidence_raw
  else df_incidence <- rbind(df_incidence, df_incidence_raw)
  
}

df_incidence


ggplot(df_incidence, aes(x=time, y=incidence)) +
  geom_line() +
  facet_wrap(~tag2, scales="free") +
  theme_bw() +
  geom_vline(aes(xintercept=end+28)) +
  geom_vline(aes(xintercept=start)) +
  geom_vline(aes(xintercept=end2), linetype="dashed", color="red") +
  lims(x=c(0,365))

ggplot(df_incidence %>% filter(n>10000), aes(x=time, y=incidence)) +
  geom_line() +
  facet_wrap(~n, scales="free") +
  theme_bw() +
  lims(x=c(0,400))

ggplot(df_incidence %>% filter(time==365, n>10000), aes(x=dd_upper, y=tot_inf)) +
  geom_point() +
  theme_bw() 

