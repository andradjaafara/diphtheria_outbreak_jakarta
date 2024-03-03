## Load some required packages
##############
library(deSolve)
library(dplyr)
library(ggplot2)
library(rstan)
library(chron)
library(bayesplot)
library(loo)
library(robustbase)
library(matrixStats)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## Model with stan

##################################################################################
# set inits vector
##################################################################################
npop = 10400000
inits <- list(S = 0,
              E = 0,
              I = 0,
              C1 = 0,
              C2 = 0,
              C3 = 0,
              R = 0,
              V = 0,
              Ev = 0,
              In = 0,
              As = 0,
              Rt = 0)

##################################################################################
# set vector of parameters
##################################################################################
params <- list(beta = 0,
               theta1 = 0,
               theta2 = 0,
               theta3 = 0,
               p1 = 0,
               p2 = 0,
               kappa = 0,
               sigmaE_1 = 0,
               sigmaE_2 = 0)

param_init <- list(beta = 0.4,
                   theta1 = 0.1,
                   theta2 = 0.05,
                   theta3 = 0.9,
                   p1 = 0.5,
                   p2 = 0.8,
                   kappa = 0.7,
                   sigmaE_1 = 0.15,
                   sigmaE_2 = 0.15,
                   S0=0.4)

##################################################################################
# model fitted for all Jakarta for using weekly data
# first vaccination campaign started at 10 Dec 2017 = day 344
##################################################################################

cases_jakarta_week <- cases_admin2_week %>% 
  group_by(date) %>% 
  summarise(cases=sum(cases)) %>% 
  ungroup() %>% 
  mutate(date_of_the_year=yday(date)) %>% # next line convert timeline of 2018 data to cont. from 2017
  mutate(date_of_the_year=ifelse(date_of_the_year<300,date_of_the_year+365,date_of_the_year))

##################################################################################
# set samples timing
##################################################################################
fake_ts = 303:428

##################################################################################
# set rstan input
##################################################################################
stan_d = list(n_obs = length(cases_jakarta_week$date_of_the_year),
              n_params = 23,
              n_difeq = length(inits),
              n_fake = length(fake_ts),
              y = cases_jakarta_week$cases,
              N = npop,
              t0 = 300,
              ts = cases_jakarta_week$date_of_the_year,
              fake_ts = fake_ts,
              n_samples = length(cases_jakarta_week$date_of_the_year),
              init=param_init)

##################################################################################
# set parameters to monitor
##################################################################################
params_monitor = c("beta","theta1","theta2","theta3","p1","p2","kappa",
                   "sigmaE_1","sigmaE_2","S0","R0","inci_samp","log_lik")

##################################################################################
# debug + run the model
##################################################################################
comp_model1 <- stan_model('codes/stan/diphtheria-model1.stan')
comp_model2 <- stan_model('codes/stan/diphtheria-model2.stan')

# runif(2,1,9999999); 6025123 8398911
seed_model1 <- 6025123
seed_model2 <- 8398911

# mod1_all = sampling(comp_model1,
#                     data = stan_d,
#                     pars = params_monitor,
#                     chains = 2,
#                     seed = seed_model1,
#                     warmup = 300,
#                     init = "random",
#                     iter = 2000)
# saveRDS(mod1_all,file="output/model/mod1_all.rds")

# mod2_all = sampling(comp_model2,
#                     data = stan_d,
#                     pars = params_monitor,
#                     chains = 2,
#                     seed = seed_model2,
#                     warmup = 300,
#                     init = "random",
#                     iter = 2000)
# saveRDS(mod2_all,file="output/model/mod2_all.rds")

mod1_all <- readRDS(file="output/model/mod1_all.rds")
mod2_all <- readRDS(file="output/model/mod2_all.rds")

##################################################################################
# Compare models
##################################################################################
log_lik_mod1 <- extract_log_lik(mod1_all, merge_chains = FALSE)
log_lik_mod2 <- extract_log_lik(mod2_all, merge_chains = FALSE)

r_eff_mod1 <- relative_eff(exp(log_lik_mod1))
r_eff_mod2 <- relative_eff(exp(log_lik_mod2))

loo_mod1 <- loo(log_lik_mod1, r_eff = r_eff_mod1, cores = 4)
loo_mod2 <- loo(log_lik_mod2, r_eff = r_eff_mod2, cores = 4)

print(loo_mod1)
print(loo_mod2)

comp <- loo_compare(loo_mod1, loo_mod2)
print(comp)

##################################################################################
# extract simulated values
##################################################################################
# fake_I.post <- extract(mod_6_new)$fake_I

mod1_posterior <- data.frame(beta=rstan::extract(mod1_all,
                                                 inc_warmup=FALSE)$beta,
                             S0=rstan::extract(mod1_all, 
                                               inc_warmup=FALSE)$S0,
                             theta1=rstan::extract(mod1_all,
                                                   inc_warmup=FALSE)$theta1,
                             theta2=rstan::extract(mod1_all,
                                                   inc_warmup=FALSE)$theta2,
                             theta3=rstan::extract(mod1_all,
                                                   inc_warmup=FALSE)$theta3,
                             p1=rstan::extract(mod1_all, 
                                               inc_warmup=FALSE)$p1,
                             p2=rstan::extract(mod1_all, 
                                               inc_warmup=FALSE)$p2,
                             kappa=rstan::extract(mod1_all, 
                                                  inc_warmup=FALSE)$kappa,
                             sigmaE_1=rstan::extract(mod1_all,
                                                     inc_warmup=FALSE)$sigmaE_1,
                             sigmaE_2=rstan::extract(mod1_all,
                                                     inc_warmup=FALSE)$sigmaE_2,
                             R0=rstan::extract(mod1_all,inc_warmup=FALSE)$R0)

mod1_summary <- data.frame(param=colnames(mod1_posterior),
                           mean=as.numeric(apply(mod1_posterior,2,mean)),
                           median=as.numeric(apply(mod1_posterior,2,quantile,
                                                   probs=0.5)),
                           lower=as.numeric(apply(mod1_posterior,2,quantile,
                                                  probs=0.025)),
                           upper=as.numeric(apply(mod1_posterior,2,quantile,
                                                  probs=0.975)))

mod2_posterior <- data.frame(beta=rstan::extract(mod2_all,
                                                 inc_warmup=FALSE)$beta,
                             S0=rstan::extract(mod2_all, 
                                               inc_warmup=FALSE)$S0,
                             theta1=rstan::extract(mod2_all,
                                                   inc_warmup=FALSE)$theta1,
                             theta2=rstan::extract(mod2_all,
                                                   inc_warmup=FALSE)$theta2,
                             theta3=rstan::extract(mod2_all,
                                                   inc_warmup=FALSE)$theta3,
                             p1=rstan::extract(mod2_all, 
                                               inc_warmup=FALSE)$p1,
                             p2=rstan::extract(mod2_all, 
                                               inc_warmup=FALSE)$p2,
                             kappa=rstan::extract(mod2_all, 
                                                  inc_warmup=FALSE)$kappa,
                             sigmaE_1=rstan::extract(mod2_all,
                                                     inc_warmup=FALSE)$sigmaE_1,
                             sigmaE_2=rstan::extract(mod2_all,
                                                     inc_warmup=FALSE)$sigmaE_2,
                             R0=rstan::extract(mod2_all,inc_warmup=FALSE)$R0)

mod2_summary <- data.frame(param=colnames(mod2_posterior),
                           mean=as.numeric(apply(mod2_posterior,2,mean)),
                           median=as.numeric(apply(mod2_posterior,2,quantile,
                                                   probs=0.5)),
                           lower=as.numeric(apply(mod2_posterior,2,quantile,
                                                  probs=0.025)),
                           upper=as.numeric(apply(mod2_posterior,2,quantile,
                                                  probs=0.975)))

mod1_samp <- rstan::extract(mod1_all, inc_warmup=FALSE)$inci_samp
mod2_samp <- rstan::extract(mod2_all, inc_warmup=FALSE)$inci_samp

mod1_median_posterior <- colMedians(as.matrix(mod1_posterior))
names(mod1_median_posterior) <- colnames(mod1_posterior)

mod2_median_posterior <- colMedians(as.matrix(mod2_posterior))
names(mod2_median_posterior) <- colnames(mod2_posterior)

mod1_quantiles_inci_samp <- colQuantiles(as.matrix(mod1_samp),
                                         probs=c(0.025,0.5,0.975))
mod2_quantiles_inci_samp <- colQuantiles(as.matrix(mod2_samp),
                                         probs=c(0.025,0.5,0.975))

# plot estim and true data
# model 1
plot(1:18,cases_jakarta_week$cases,col="black",cex=1.5,
     ylim=c(0,max(mod1_quantiles_inci_samp)),pch=16)
lines(1:18,mod1_quantiles_inci_samp[,2],lwd=2,col="red")
lines(1:18,mod1_quantiles_inci_samp[,1],lwd=2,lty=2,col="red")
lines(1:18,mod1_quantiles_inci_samp[,3],lwd=2,lty=2,col="red")
# model 2
plot(1:18,cases_jakarta_week$cases,col="black",cex=1.5,
     ylim=c(0,max(mod2_quantiles_inci_samp)),pch=16)
lines(1:18,mod2_quantiles_inci_samp[,2],lwd=2,col="red")
lines(1:18,mod2_quantiles_inci_samp[,1],lwd=2,lty=2,col="red")
lines(1:18,mod2_quantiles_inci_samp[,3],lwd=2,lty=2,col="red")

##################################################################################
# get traceplot and parameter estimates from model
##################################################################################
## traceplots
p_trace_mod1 <- stan_trace(mod1_all, pars = c("beta","theta1","theta2",
                                              "theta3","p1","p2","kappa",
                                              "sigmaE_1","sigmaE_2","S0","R0"),
                           inc_warmup = TRUE,ncol=3) +
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank())

p_trace_mod2 <- stan_trace(mod2_all, pars = c("beta","theta1","theta2",
                                              "theta3","p1","p2","kappa",
                                              "sigmaE_1","sigmaE_2","S0","R0"),
                           inc_warmup = TRUE,ncol=3) +
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank())

p_trace_all <- plot_grid(p_trace_mod1,p_trace_mod2,
                         labels=c("A","B","C"),ncol=1)

ggsave("output/figures/p_trace.png", p_trace_all,
       width = 16, height = 26, units = "cm")

## parameter estimates
# simulate having posteriors for two different models each with parameters beta[1],..., beta[4]
# posterior_1 <- matrix(rnorm(4000), 1000, 4)
# posterior_2 <- matrix(rnorm(4000), 1000, 4)
# posterior_3 <- matrix(rnorm(4000), 1000, 4)
# colnames(posterior_1) <- colnames(posterior_2) <- paste0("beta[", 1:4, "]")

# use bayesplot::mcmc_intervals_data() function to get intervals data in format easy to pass to ggplot
combined_posteriors <- rbind(mcmc_intervals_data(mod1_all,prob_outer = 0.95),
                             mcmc_intervals_data(mod2_all,prob_outer = 0.95)) %>% 
  filter(parameter %in% c("beta","theta1","theta2",
                          "theta3","p1","p2","kappa",
                          "sigmaE_1","sigmaE_2","S0","R0")) %>% 
  mutate(parameter=factor(parameter,levels=c("theta1","theta2",
                                             "theta3","p1","p2","kappa",
                                             "sigmaE_1","sigmaE_2","beta","S0","R0")))

combined_posteriors$model <- rep(c("Model 1", "Model 2"), 
                                 each = length(c("beta","theta1","theta2",
                                                 "theta3","p1","p2","kappa",
                                                 "sigmaE_1","sigmaE_2","S0","R0")))

# make the plot using ggplot 
pos <- position_nudge(y = ifelse(combined_posteriors$model == "Model 1", 0.125,
                          ifelse(combined_posteriors$model == "Model 2",-0.125,0)))

p_estimates <- ggplot(combined_posteriors, aes(x = m, y = parameter, color = model)) + 
  geom_point(position = pos, size=2) +
  geom_linerange(aes(xmin = ll, xmax = hh), position = pos) +
  geom_linerange(aes(xmin = l, xmax = h), position = pos, linewidth=1.25) +
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank()) +
  labs(col="Model",x="Parameter estimates",y="Model parameters") +
  theme(legend.position = "top") +
  scale_colour_manual(breaks=c("Model 1","Model 2"),
                      labels=c("Model 1","Model 2"),
                      values=safe_colorblind_palette[c(7,2)])

ggsave("output/figures/p_estimates.png", p_estimates,
       width = 16, height = 12, units = "cm")

# tables of parameter estimates

params_summary_mod1 <- rownames_to_column(data.frame(summary(mod1_all,
                                                             params_monitor[1:11])$summary))
params_summary_mod2 <- rownames_to_column(data.frame(summary(mod2_all,
                                                             params_monitor[1:11])$summary))

write_csv(params_summary_mod1,"output/model/parameter_estimates_mod1.csv")
write_csv(params_summary_mod2,"output/model/parameter_estimates_mod2.csv")





