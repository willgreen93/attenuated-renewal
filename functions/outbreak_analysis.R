library(ggplot2)
library(zoo)
library(dplyr)

simul_list <- list.files("simulation2/MSM_like/", full.names=T)

for(i in simul_list[1:100]){
  simulation <- readRDS(i)
  
  print(sum(simulation$outbreak$incidence_vector))
  
  prop_infected <- data.frame(day=1:simulation$outbreak$outbreak_params$timesteps,
                              incidence=simulation$outbreak$incidence_vector,
                              prop_susc=1-(cumsum(simulation$outbreak$incidence_vector)/simulation$graph$graph_params$n))
  
  R_case <- simulation$outbreak$generation_interval_matrix %>% left_join(., prop_infected, by=c("exposed_time"="day")) %>% 
    group_by(index, infector, infector_infectious, exposed_time, prop_susc) %>% summarise(count=n()) %>%
    mutate(adjusted_count=count/prop_susc)  %>% 
    group_by(infector, infector_infectious) %>% summarise(count=sum(adjusted_count)) %>%
    group_by(infector_infectious) %>% summarise(R_adj=mean(count)) 
  
  plotting <- left_join(prop_infected, R_case, by=c("day"="infector_infectious")) %>% mutate(R_adj=ifelse(is.na(R_adj), 0, R_adj)) %>% mutate(R_adj_m=rollmean(R_adj, 7, fill=NA)) %>% filter(incidence>10) %>% mutate(tag=factor(which(simul_list==i)))
  
  if(i==simul_list[1]) plotting_df <- plotting
  else plotting_df <- rbind(plotting, plotting_df)
  
  a1 <- ggplot(plotting, aes(x=day, y=R_adj_m)) + geom_line() + theme_bw() + theme(panel.grid=element_blank()) + scale_x_continuous(limits = c(0, NA))
  a2 <- ggplot(plotting, aes(x=1-prop_susc, y=R_adj_m)) + geom_line() + theme_bw() + theme(panel.grid=element_blank())
  
  print(ggpubr::ggarrange(a1, a2))
  
} 

a1 <- ggplot(plotting_df, aes(x=day, y=R_adj_m, color=tag)) + geom_line() + theme_bw() + theme(panel.grid=element_blank()) + scale_x_continuous(limits = c(0, NA))
a2 <- ggplot(plotting_df, aes(x=1-prop_susc, y=R_adj_m, color=tag)) + geom_line() + theme_bw() + theme(panel.grid=element_blank())

plotting_df %>% group_by()
