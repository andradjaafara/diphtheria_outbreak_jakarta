library(raster)

### This file contains analyses evaluating contributions of:
### 1. infection type to Rt
### 2. susceptibility to transmission dynamics
### 3. contact tracing to Rt

### 1. Start of evaluation of the contribution of each infection type to Rt!

beta = as.numeric(mod_optimistic_summary[1,3])
S0 = as.numeric(mod_optimistic_summary[2,3])
theta1 = as.numeric(mod_optimistic_summary[3,3])
theta2 = as.numeric(mod_optimistic_summary[4,3])
theta3 = as.numeric(mod_optimistic_summary[5,3])
p1 = as.numeric(mod_optimistic_summary[6,3])
p2 = as.numeric(mod_optimistic_summary[7,3])
kappa = as.numeric(mod_optimistic_summary[8,3])
sigmaE_1 = as.numeric(mod_optimistic_summary[9,3])
sigmaE_2 = as.numeric(mod_optimistic_summary[10,3])
sigmaC_1 = as.numeric(mod_optimistic_summary[11,3])
sigmaC_2 = as.numeric(mod_optimistic_summary[12,3])
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
sigmaC_2 <- 0
sigmaE_2 <- 0

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
  rt_C1[i] <- ((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(theta1*beta*kappa1/(gammaC1+sigmaC_2))
  rt_C2[i] <- (delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*((1-p2-eta)*gammaC1/(p2*gammaI_1+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(gammaC1/(gammaC1+sigmaC_2)))*(theta2*beta*kappa1/(gammaC1)) 
  rt_C3[i] <- (delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*((1-p2-eta)*gammaC1/(p2*gammaI_1+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(gammaC1/(gammaC1+sigmaC_2)))*(theta2*beta*kappa1/(gammaC3))*(1-epsilon)
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

p4 <- ggplot(contribution.df,aes(x=immune*100)) +
  geom_line(aes(y=Carriers.contrib*100,col="Carriers/Asymptomatic"),linewidth=1.5) +
  geom_line(aes(y=Exposed.contrib*100,col="Exposed"),linewidth=1.5) +
  geom_line(aes(y=Clinical.contrib*100,col="Clinical/Symptomatic"),linewidth=1.5) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="Proportion of the immune population (%)",
       y="Contribution of each infectious component (%)",
       col=NULL) +
  scale_colour_manual(breaks=c("Carriers/Asymptomatic", "Exposed", "Clinical/Symptomatic"),
                      labels=c("Carriers/Asymptomatic", "Exposed", "Clinical/Symptomatic"),
                      values=c("blue","purple","red")) +
  theme(legend.position=c(0.175,0.875)) + ylim(0,80) +
  geom_vline(xintercept=84,linewidth=1.5,linetype=2,col="green") +
  annotate("text",x=84,y=25,label="27%") +
  annotate("text",x=84,y=69,label="66%") +
  annotate("text",x=84,y=10,label="7%")

ggsave("figures/new/p4.png", p4,
       width = 16, height = 12, units = "cm")

### 2. Start of susceptibility evaluation

sigmaC_2 <- 0
sigmaE_2 <- 0

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
    ((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(theta1*beta*kappa1/(gammaC1+sigmaC_2)) +
    (delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*((1-p2-eta)*gammaC1/(p2*gammaI_1+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(gammaC1/(gammaC1+sigmaC_2)))*(theta2*beta*kappa1/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*((1-p2-eta)*gammaC1/(p2*gammaI_1+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(gammaC1/(gammaC1+sigmaC_2)))*(theta2*beta*kappa1/(gammaC3))*(1-epsilon)
  
  rti <- S/npop*delta*(S/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE_2) +
                         V/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE_2) +
                         delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*beta*kappa1/(p2*gammaI_1+(1-p2)*gammaC1) +
                         ((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(theta1*beta*kappa1/(gammaC1+sigmaC_2)) +
                         (delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*((1-p2-eta)*gammaC1/(p2*gammaI_1+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(gammaC1/(gammaC1+sigmaC_2)))*(theta2*beta*kappa1/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*((1-p2-eta)*gammaC1/(p2*gammaI_1+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(gammaC1/(gammaC1+sigmaC_2)))*(theta2*beta*kappa1/(gammaC3))*(1-epsilon))
  
  prop[i] <- S/npop*delta
  
  rti_save[i] <- rti
  diff_rt[i] <- rt_save[i] - rt_save[1]
}

pop_immunity_impact <- data.frame(susceptible=s0,
                                  immune=v0,
                                  diff_rt,diff_rt,
                                  prop=prop)

p_immune_rt <- pop_immunity_impact %>% 
  ggplot(aes(x=immune*100,y=diff_rt)) +
  geom_line(col="red",linewidth=1.5) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="Proportion of immune population\n before the outbreak (%)",
       y="Reff value differences\n to 100% susceptible population")

p_immune_symptomatic <- pop_immunity_impact %>% 
  ggplot(aes(x=immune*100,y=prop*100)) +
  geom_line(col="red",linewidth=1.5) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="Proportion of immune population\n before the outbreak (%)",
       y="Proportion of infections\n to be symptomatic cases (%)")
  
p_immune_impact <- plot_grid(p_immune_rt,p_immune_symptomatic,ncol=2,
                             labels=c("A","B"))
ggsave("figures/new/p_immune_impact.png", p_immune_impact,
       width = 16, height = 8, units = "cm")

### 3. Start of CT evaluation

coverage <- seq(0,0.99,by = 0.01)

sigmaC_param <- gammaC1*coverage/(1-coverage)
sigmaE_param <- tau*coverage/(1-coverage)

param_grid <- expand.grid(sigmaC_param,rev(sigmaE_param))

Rt_estim_ct_1 <- data.frame(sigmaC=NA,sigmaE=NA,covC=NA,covE=NA,Rt=NA,lt1=NA) # data for ct coverage map
Rt_estim_ct_2 <- data.frame(sigmaC=NA,sigmaE=NA,covC=NA,covE=NA,Rt=NA,lt1=NA) # data for ct coverage map
Rt_estim_ct_3 <- data.frame(sigmaC=NA,sigmaE=NA,covC=NA,covE=NA,Rt=NA,lt1=NA) # data for ct coverage map
Rt_estim_ct_4 <- data.frame(sigmaC=NA,sigmaE=NA,covC=NA,covE=NA,Rt=NA,lt1=NA) # data for ct coverage map

# 1. 100%
npop = 10400000
I  = 1
S  = as.numeric(npop * 1 - I)
E  = 0
C1 = 0
C2 = 0
C3 = 0
R  = 0     
V  = as.numeric((npop - S - I))
Ev = 0

for (i in 1:nrow(param_grid)){
  sigmaC_2 <- param_grid[i,1]
  sigmaE_2 <- param_grid[i,2]
  covC <- sigmaC_2/(gammaC1+sigmaC_2)
  covE <- sigmaE_2/(tau+sigmaE_2)
  rt <- S/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE_2) + 
    V/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE_2) + 
    delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*beta*kappa1/(p2*gammaI_1+(1-p2)*gammaC1) + 
    ((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(theta1*beta*kappa1/(gammaC1+sigmaC_2)) + 
    (delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*((1-p2-eta)*gammaC1/(p2*gammaI_1+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(gammaC1/(gammaC1+sigmaC_2)))*(theta2*beta*kappa1/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*((1-p2-eta)*gammaC1/(p2*gammaI_1+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(gammaC1/(gammaC1+sigmaC_2)))*(theta2*beta*kappa1/(gammaC3))*(1-epsilon)
  lt1 <- as.numeric(rt<1)
  Rt_estim_ct_1[i,] <- c(sigmaC_2,sigmaE_2,covC,covE,rt,lt1)
}

# 2. 50%
npop = 10400000
I  = 1
S  = as.numeric(npop * 0.5 - I)
E  = 0
C1 = 0
C2 = 0
C3 = 0
R  = 0     
V  = as.numeric((npop - S - I))
Ev = 0

for (i in 1:nrow(param_grid)){
  sigmaC_2 <- param_grid[i,1]
  sigmaE_2 <- param_grid[i,2]
  covC <- sigmaC_2/(gammaC1+sigmaC_2)
  covE <- sigmaE_2/(tau+sigmaE_2)
  rt <- S/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE_2) + 
    V/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE_2) + 
    delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*beta*kappa1/(p2*gammaI_1+(1-p2)*gammaC1) + 
    ((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(theta1*beta*kappa1/(gammaC1+sigmaC_2)) + 
    (delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*((1-p2-eta)*gammaC1/(p2*gammaI_1+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(gammaC1/(gammaC1+sigmaC_2)))*(theta2*beta*kappa1/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*((1-p2-eta)*gammaC1/(p2*gammaI_1+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(gammaC1/(gammaC1+sigmaC_2)))*(theta2*beta*kappa1/(gammaC3))*(1-epsilon)
  lt1 <- as.numeric(rt<1)
  Rt_estim_ct_2[i,] <- c(sigmaC_2,sigmaE_2,covC,covE,rt,lt1)
}

# 3. 25%
npop = 10400000
I  = 1
S  = as.numeric(npop * 0.25 - I)
E  = 0
C1 = 0
C2 = 0
C3 = 0
R  = 0     
V  = as.numeric((npop - S - I))
Ev = 0

for (i in 1:nrow(param_grid)){
  sigmaC_2 <- param_grid[i,1]
  sigmaE_2 <- param_grid[i,2]
  covC <- sigmaC_2/(gammaC1+sigmaC_2)
  covE <- sigmaE_2/(tau+sigmaE_2)
  rt <- S/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE_2) + 
    V/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE_2) + 
    delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*beta*kappa1/(p2*gammaI_1+(1-p2)*gammaC1) + 
    ((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(theta1*beta*kappa1/(gammaC1+sigmaC_2)) + 
    (delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*((1-p2-eta)*gammaC1/(p2*gammaI_1+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(gammaC1/(gammaC1+sigmaC_2)))*(theta2*beta*kappa1/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*((1-p2-eta)*gammaC1/(p2*gammaI_1+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(gammaC1/(gammaC1+sigmaC_2)))*(theta2*beta*kappa1/(gammaC3))*(1-epsilon)
  lt1 <- as.numeric(rt<1)
  Rt_estim_ct_3[i,] <- c(sigmaC_2,sigmaE_2,covC,covE,rt,lt1)
}

# 4. 10%
npop = 10400000
I  = 1
S  = as.numeric(npop * 0.1 - I)
E  = 0
C1 = 0
C2 = 0
C3 = 0
R  = 0     
V  = as.numeric((npop - S - I))
Ev = 0

for (i in 1:nrow(param_grid)){
  sigmaC_2 <- param_grid[i,1]
  sigmaE_2 <- param_grid[i,2]
  covC <- sigmaC_2/(gammaC1+sigmaC_2)
  covE <- sigmaE_2/(tau+sigmaE_2)
  rt <- S/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE_2) + 
    V/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE_2) + 
    delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*beta*kappa1/(p2*gammaI_1+(1-p2)*gammaC1) + 
    ((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(theta1*beta*kappa1/(gammaC1+sigmaC_2)) + 
    (delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*((1-p2-eta)*gammaC1/(p2*gammaI_1+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(gammaC1/(gammaC1+sigmaC_2)))*(theta2*beta*kappa1/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*((1-p2-eta)*gammaC1/(p2*gammaI_1+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(gammaC1/(gammaC1+sigmaC_2)))*(theta2*beta*kappa1/(gammaC3))*(1-epsilon)
  lt1 <- as.numeric(rt<1)
  Rt_estim_ct_4[i,] <- c(sigmaC_2,sigmaE_2,covC,covE,rt,lt1)
}

r_1 <- raster(xmn = 0, xmx = 0.99, ymn = 0, ymx = 0.99, nrows = 100, ncols = 100)
r_1[] <- Rt_estim_ct_1$lt1
r_2 <- raster(xmn = 0, xmx = 0.99, ymn = 0, ymx = 0.99, nrows = 100, ncols = 100)
r_2[] <- Rt_estim_ct_2$lt1
r_3 <- raster(xmn = 0, xmx = 0.99, ymn = 0, ymx = 0.99, nrows = 100, ncols = 100)
r_3[] <- Rt_estim_ct_3$lt1
r_4 <- raster(xmn = 0, xmx = 0.99, ymn = 0, ymx = 0.99, nrows = 100, ncols = 100)
r_4[] <- Rt_estim_ct_4$lt1

p_ct_100 <- ggplot(Rt_estim_ct_1, aes(covC, covE)) +
  geom_raster(aes(fill = as.character(lt1))) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="Proportion of carriers\n removed by contact tracing (%)",
       y="Proportion of exposed individuals\n removed by contact tracing (%)",
       fill=NULL) +
  theme(legend.position = "top") +
  scale_fill_manual(breaks=c("0","1"),
                    labels=c("Rt >= 1","Rt < 1"),
                    values=rev(terrain.colors(2)))
p_ct_50 <- ggplot(Rt_estim_ct_2, aes(covC, covE)) +
  geom_raster(aes(fill = as.character(lt1))) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="Proportion of carriers\n removed by contact tracing (%)",
       y="Proportion of exposed individuals\n removed by contact tracing (%)",
       fill=NULL) +
  theme(legend.position = "top") +
  scale_fill_manual(breaks=c("0","1"),
                    labels=c("Rt >= 1","Rt < 1"),
                    values=rev(terrain.colors(2)))
p_ct_25 <- ggplot(Rt_estim_ct_3, aes(covC, covE)) +
  geom_raster(aes(fill = as.character(lt1))) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="Proportion of carriers\n removed by contact tracing (%)",
       y="Proportion of exposed individuals\n removed by contact tracing (%)",
       fill=NULL) +
  theme(legend.position = "top") +
  scale_fill_manual(breaks=c("0","1"),
                    labels=c("Rt >= 1","Rt < 1"),
                    values=rev(terrain.colors(2)))
p_ct_10 <- ggplot(Rt_estim_ct_4, aes(covC, covE)) +
  geom_raster(aes(fill = as.character(lt1))) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="Proportion of carriers\n removed by contact tracing (%)",
       y="Proportion of exposed individuals\n removed by contact tracing (%)",
       fill=NULL) +
  theme(legend.position = "top") +
  scale_fill_manual(breaks=c("0","1"),
                    labels=c("Rt >= 1","Rt < 1"),
                    values=rev(terrain.colors(2)))
p_ct_legend <- get_legend(p_ct_100)

p_ct_all <- plot_grid(p_ct_100 + theme(legend.position = "none"),
                      p_ct_50 + theme(legend.position = "none"),
                      p_ct_25 + theme(legend.position = "none"),
                      p_ct_10 + theme(legend.position = "none"),
                      labels=c("A","B","C","D"),ncol=2)
p_ct_impact <- plot_grid(p_ct_all,ggdraw(p_ct_legend),
                         rel_heights = c(1,0.1),ncol=1)
ggsave("figures/new/p_ct_impact.png", p_ct_impact,
       width = 16, height = 18, units = "cm")