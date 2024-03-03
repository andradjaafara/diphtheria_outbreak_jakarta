### This file contains analyses evaluating contributions of:
### 1. infection type to Rt
### 2. susceptibility to transmission dynamics

### 1. Start of evaluation of the contribution of each infection type to Rt!

beta = as.numeric(mod1_summary[1,3])
S0 = as.numeric(mod1_summary[2,3])
theta1 = as.numeric(mod1_summary[3,3])
theta2 = as.numeric(mod1_summary[4,3])
theta3 = as.numeric(mod1_summary[5,3])
p1 = as.numeric(mod1_summary[6,3])
p2 = as.numeric(mod1_summary[7,3])
kappa = as.numeric(mod1_summary[8,3])
sigmaE_1 = as.numeric(mod1_summary[9,3])
sigmaE_2 = as.numeric(mod1_summary[10,3])
rho =  40000 * S0 * 1.974 # 40000 roughly average per day vaccinated; times susceptible times 1.974531 ratio S in children to overall S
tau = 1/3
gammaI_1 = 1/(3.88+2)
gammaI_2 = 1/(1.12+2)
gammaC1 = 1/18
gammaC3 = 1/18
delta = 0.70
eta = 0.05
rho1 = 0
kappa1 = 1
epsilon = 0.95
sigmaE = 0
sigmaC = 0

s0 <- seq(1,0,by=-0.01)
v0 <- seq(0,1,by=0.01)
rt_E <- vector()
rt_EV <- vector()
rt_I <- vector()
rt_C1 <- vector()
rt_C2 <- vector()
rt_C3 <- vector()
rt <- vector()
rt_exposed <- vector()
rt_carriers <- vector()
contrib_exposed <- vector()
contrib_carriers <- vector()
contrib_clin <- vector()

for (i in 1:length(s0)){
  npop = 10400000
  I  = 1
  S  = as.numeric(npop * s0[i] - I)
  E  = 0
  C1 = 0
  C2 = 0
  C3 = 0
  R  = 0     
  V  = as.numeric((npop - S - I))
  Ev = 0
  
  rt_E[i] <- S/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE_2)
  rt_EV[i] <- V/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE_2)
  rt_I[i] <- delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*beta*kappa1/(p2*gammaI_1+(1-p2)*gammaC1)
  rt_C1[i] <- ((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(theta1*beta*kappa1/(gammaC1))
  rt_C2[i] <- (delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*((1-p2-eta)*gammaC1/(p2*gammaI_1+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa1/(gammaC1)) 
  rt_C3[i] <- (delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*((1-p2-eta)*gammaC1/(p2*gammaI_1+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa1/(gammaC3))*(1-epsilon)
  rt[i] <- rt_E[i] + rt_EV[i] + rt_I[i] + rt_C1[i] + rt_C2[i] + rt_C3[i]
  rt_exposed[i] <- rt_E[i] + rt_EV[i]
  rt_carriers[i] <- rt_C1[i] + rt_C2[i] + rt_C3[i]
  contrib_exposed[i] <- rt_exposed[i]/rt[i]
  contrib_carriers[i] <- rt_carriers[i]/rt[i]
  contrib_clin[i] <- rt_I[i]/rt[i]
}

contribution.df <- data.frame(S0=s0,E=rt_E,EV=rt_EV,I=rt_I,C1=rt_C1,C2=rt_C2,C3=rt_C3,Rt_all=rt,
                              Exposed=rt_exposed,Carriers=rt_carriers,Clinical=rt_I,
                              Exposed.contrib=contrib_exposed,Carriers.contrib=contrib_carriers,
                              Clinical.contrib=contrib_clin,
                              immune=v0)

### ggplot friendly format
contribution.df2 <- contribution.df %>% 
  dplyr::select(immune,`Carriers/Asymptomatic`=Carriers.contrib,
                Exposed=Exposed.contrib,`Clinical/Symptomatic`=Clinical.contrib) %>% 
  pivot_longer(-immune,names_to="component",values_to="proportion") %>% 
  mutate(immune=immune*100,
         proportion=proportion*100,
         component=factor(component,
                          levels=c("Carriers/Asymptomatic","Exposed","Clinical/Symptomatic")))

p4_A <- contribution.df2 %>% 
  ggplot(aes(x=immune,y=proportion,fill=component)) +
  geom_area() + 
  scale_fill_manual(values=safe_colorblind_palette[c(7:9)]) +
  theme_classic() +
  labs(x="Proportion of the immune population (%)",
       y="Contribution of each infectious component (%)",
       fill=NULL) +
  theme(legend.position = "bottom") +
  # geom_vline(xintercept=87,linewidth=1.5,linetype=2,col="black") +
  geom_segment(aes(x=87,xend=87,y=0,yend=100),linewidth=1.25,linetype=2,col="white") +
  annotate("text",x=92,y=75,label="50.5%",col="white") +
  annotate("text",x=92,y=33,label="41.0%",col="white") +
  annotate("text",x=91,y=4,label="8.5%",col="white") + labs(tag="A")
  
# p4_A <- ggplot(contribution.df,aes(x=immune*100)) +
#   geom_line(aes(y=Carriers.contrib*100,col="Carriers/Asymptomatic"),linewidth=1.5) +
#   geom_line(aes(y=Exposed.contrib*100,col="Exposed"),linewidth=1.5) +
#   geom_line(aes(y=Clinical.contrib*100,col="Clinical/Symptomatic"),linewidth=1.5) +
#   theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   labs(x="Proportion of the immune population (%)",
#        y="Contribution of each infectious component (%)",
#        col=NULL) +
#   scale_colour_manual(breaks=c("Carriers/Asymptomatic", "Exposed", "Clinical/Symptomatic"),
#                       labels=c("Carriers/Asymptomatic", "Exposed", "Clinical/Symptomatic"),
#                       values=c("blue","purple","red")) +
#   theme(legend.position=c(0.175,0.875)) + ylim(0,80) +
#   geom_vline(xintercept=87,linewidth=1.5,linetype=2,col="green") +
#   annotate("text",x=87,y=50.5,label="50.5%") +
#   annotate("text",x=87,y=41,label="41.0%") +
#   annotate("text",x=87,y=8.5,label="8.5%") + labs(tag="A") +
#   theme(legend.background = element_rect(fill='transparent'))

### 2. Start of susceptibility evaluation
s0 <- seq(1,0,by=-0.01)
v0 <- seq(0,1,by=0.01)
diff_rt <- vector()
rti_save <- vector()
rt_save <- vector()
prop <- vector()

# R0 <- 4.649771

for (i in 1:length(s0)){
  npop = 10400000
  I  = 1
  S  = as.numeric(npop * s0[i] - I)
  E  = 0
  C1 = 0
  C2 = 0
  C3 = 0
  R  = 0     
  V  = as.numeric((npop - S - I))
  Ev = 0
  
  rt_save[i] <- S/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE_2) +
    V/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE_2) +
    delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*beta*kappa1/(p2*gammaI_1+(1-p2)*gammaC1) +
    ((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(theta1*beta*kappa1/(gammaC1)) +
    (delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*((1-p2-eta)*gammaC1/(p2*gammaI_1+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa1/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*((1-p2-eta)*gammaC1/(p2*gammaI_1+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa1/(gammaC3))*(1-epsilon)
  
  rti <- S/npop*delta*(S/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE_2) +
                         V/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE_2) +
                         delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*beta*kappa1/(p2*gammaI_1+(1-p2)*gammaC1) +
                         ((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(theta1*beta*kappa1/(gammaC1)) +
                         (delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*((1-p2-eta)*gammaC1/(p2*gammaI_1+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa1/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*((1-p2-eta)*gammaC1/(p2*gammaI_1+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa1/(gammaC3))*(1-epsilon))
  
  prop[i] <- S/npop*delta
  
  rti_save[i] <- rti
  diff_rt[i] <- rt_save[i] - rt_save[1]
}

pop_immunity_impact <- data.frame(susceptible=s0,
                                  immune=v0,
                                  diff_rt,diff_rt,
                                  prop=prop)

p4_B <- pop_immunity_impact %>% 
  ggplot(aes(x=immune*100,y=diff_rt)) +
  geom_line(col="red",linewidth=1.5) +
  theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="Proportion of immune population\n before the outbreak (%)",
       y="Reff value differences\n to 100% susceptible population") + labs(tag="B")

p4_C <- pop_immunity_impact %>% 
  ggplot(aes(x=immune*100,y=prop*100)) +
  geom_line(col="red",linewidth=1.5) +
  theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="Proportion of immune population\n before the outbreak (%)",
       y="Proportion of infections\n to be symptomatic cases (%)") + labs(tag="C")

  
p4_BC <- plot_grid(p4_B,p4_C,ncol=2)
p4 <- plot_grid(p4_A,p4_BC,nrow=2,rel_heights = c(1,0.75))

ggsave("output/figures/p4.png", p4,
       width = 16, height = 20, units = "cm")