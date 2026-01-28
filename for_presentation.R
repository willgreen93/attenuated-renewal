library(ggplot2)

q1 <- ggplot(skeleton_plot_fil_10 %>% filter(true_value>0 | day > cutoff), aes(x=day, fill=factor(epi_phase))) +
  geom_point(aes(y=true_value)) +
  #geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
  #geom_line(aes(y=median)) +
  #geom_vline(aes(xintercept=cutoff), linetype="dashed") +
  theme_minimal() +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm"),
        panel.grid=element_blank(),
        legend.position="none") +
  #scale_y_log10() +
  labs(x="Day", y="Cases") +
  #facet_grid(model~place, scales="free") +
  theme(panel.grid=element_blank(),
        panel.background = element_rect(fill = 'white')) 

q1

q2 <- ggplot(skeleton_plot_fil_10 %>% filter(true_value>0 | day > cutoff), aes(x=day, fill=factor(epi_phase))) +
  geom_point(aes(y=true_value)) +
  #geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
  #geom_line(aes(y=median)) +
  #geom_vline(aes(xintercept=cutoff), linetype="dashed") +
  theme_minimal() +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm"),
        panel.grid=element_blank(),
        legend.position="none") +
  geom_vline(aes(xintercept=70), linetype="dashed") +
  labs(x="Day", y="Cases") +
  #facet_grid(model~place, scales="free") +
  theme(panel.grid=element_blank(),
        panel.background = element_rect(fill = 'white')) 

q2

q3 <- ggplot(skeleton_plot_fil_5 %>% filter(true_value>0 | day > cutoff) %>% mutate(true_value=ifelse(day>cutoff, NA, true_value)), aes(x=day, fill=factor(epi_phase))) +
  geom_point(aes(y=true_value)) +
  #geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
  #geom_line(aes(y=median)) +
  geom_vline(aes(xintercept=cutoff), linetype="dashed") +
  theme_minimal() +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm"),
        panel.grid=element_blank(),
        legend.position="none") +
  #scale_y_log10() +
  labs(x="Day", y="Cases") +
  xlim(c(0,139)) +
  #facet_grid(model~place, scales="free") +
  theme(panel.grid=element_blank(),
        panel.background = element_rect(fill = 'white')) 

q3

q4 <- ggplot(skeleton_plot_fil_5 %>% filter(true_value>0 | day > cutoff) %>% mutate(true_value=ifelse(day>cutoff, NA, true_value), median=ifelse(day<cutoff, NA, median), lower=ifelse(day<cutoff, NA, lower), upper=ifelse(day<cutoff, NA, upper)), aes(x=day, fill=factor(epi_phase))) +
  geom_point(aes(y=true_value)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
  #geom_line(aes(y=median)) +
  geom_vline(aes(xintercept=cutoff), linetype="dashed") +
  theme_minimal() +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm"),
        panel.grid=element_blank(),
        legend.position="none") +
  #scale_y_log10() +
  labs(x="Day", y="Cases") +
  xlim(c(0,139)) +
  #facet_grid(model~place, scales="free") +
  theme(panel.grid=element_blank(),
        panel.background = element_rect(fill = 'white')) 

q4

q5 <- ggplot(skeleton_plot_fil_5 %>% filter(true_value>0 | day > cutoff) %>% mutate(lower=ifelse(day<cutoff, NA, lower), upper=ifelse(day<cutoff, NA, upper)), aes(x=day, fill=factor(epi_phase))) +
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
  labs(x="Day", y="Cases") +
  xlim(c(0,139)) +
  #facet_grid(model~place, scales="free") +
  theme(panel.grid=element_blank(),
        panel.background = element_rect(fill = 'white')) 

q5

ggsave("figures/q1.png", q1, height=4, width=4)
ggsave("figures/q2.png", q2, height=4, width=4)
ggsave("figures/q3.png", q3, height=4, width=4)
ggsave("figures/q4.png", q4, height=4, width=4)
ggsave("figures/q5.png", q5, height=4, width=4)

df <- data.frame(x=0:30, y=dgamma(0:30, shape=12, rate=0.8))


q6 <- ggplot(df, aes(x=x, y=y)) +
  geom_col(fill="blue") +
  theme_bw() +
  theme(panel.grid=element_blank())

q6
ggsave("figures/q6.png", q6, height=4, width=4)

df2 <- data.frame(x=seq(0,50,1)) %>% mutate(y=exp(-((x-35)/15)^2))

q7 <- ggplot(df2, aes(x=x, y=y)) +
  geom_col(fill="orange") +
  theme_bw() +
  theme(panel.grid=element_blank())

q7
ggsave("figures/q7.png", q7, height=4, width=4)

r1 <- ggplot(data.frame(x=c(1,2,3), y=c(0.25, 0.5, 0.25)), aes(x=x, y=y)) +
  geom_col(fill="blue") +
  theme_bw() +
  labs(x="Days into infection", y="Proportion of infectiousness") +
  theme(panel.grid=element_blank())

ggsave("figures/r1.png", r1, height=3, width=5)


R_types <- data.frame(x=1:100) %>% 
  mutate(cuminc=cumsum(50*exp(-0.005*(x-50)^2))/5000) %>%
  mutate(EXPT=ifelse(x>30, 4*exp(-0.04*(x-30)), 4),
         SIGT=0.1+(4-0.1)/(1+exp(-0.15*(40-x))), 
         POWI=4*(1+cuminc)^-7,
         EXPI=4*exp(-10*cuminc))

R_types$RRAND <- NA
R_types$RBIAS <- NA
R_types$RRAND[1] <- 4
R_types$RBIAS[1] <- 4

set.seed(3)

for(i in 2:nrow(R_types)){
  R_types$RRAND[i] <- ifelse(i%%5 ==0, exp(rnorm(n=1, log(R_types$RRAND[i-4]), sd=0.1)), R_types$RRAND[i-1])
  R_types$RBIAS[i] <- ifelse(i%%5 == 0, exp(rnorm(n=1, log(R_types$RBIAS[i-4]-0.1), sd=0.1)), R_types$RBIAS[i-1])
}

R_plot <- R_types %>% tidyr::gather(model, value, 3:8) %>% mutate(model=factor(model, levels=c("RRAND", "RBIAS", "EXPT", "SIGT")))

Rrand_plot <- ggplot(R_plot %>% filter(model %in% c("RRAND")), aes(x=x, y=value)) +
  geom_line() +
  theme_bw() +
  labs(x="Day", y="Effective R") +
  theme(panel.grid=element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_rect(fill = "grey90", colour = NA),
        strip.placement = "outside") +
  facet_wrap(~model, nrow=4)

ggsave("figures/Rrand_plot.png", Rrand_plot, height=5*1/3, width=3)

Rt_plot <- ggplot(R_plot %>% filter(model %in% c("RBIAS", "EXPT", "SIGT")), aes(x=x, y=value)) +
  geom_line() +
  theme_bw() +
  labs(x="Day", y="Effective R") +
  theme(panel.grid=element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_rect(fill = "grey90", colour = NA),
        strip.placement = "outside") +
  facet_wrap(~model, nrow=4)

ggsave("figures/Rt_plot.png", Rt_plot, height=5*3/4, width=3)

Rinc_plot <- ggplot(R_plot %>% filter(!model %in% c("RRAND", "RBIAS", "EXPT", "SIGT")), aes(x=cuminc, y=value)) +
  geom_line() +
  theme_bw() +
  labs(x="Attack rate", y="Effective R") +
  theme(panel.grid=element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_rect(fill = "grey90", colour = NA),
        strip.placement = "outside") +
  facet_wrap(~model, nrow=4)

ggsave("figures/Rinc_plot.png", Rinc_plot, height=2.5, width=3)
