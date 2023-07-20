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
comp_model_optimistic <- stan_model('codes/stan/diphtheria-model-optimistic.stan')
comp_model_pessimistic <- stan_model('codes/stan/diphtheria-model-pessimistic.stan')

# runif(2,1,9999999); 6813636 6592353
seed_optimistic <- 6813636
seed_pessimistic <- 6592353

# mod_optimistic_all = sampling(comp_model_optimistic,
#                               data = stan_d,
#                               pars = params_monitor,
#                               chains = 2,
#                               seed = seed_optimistic,
#                               warmup = 300,
#                               init = "random",
#                               iter = 2000)
# saveRDS(mod_optimistic_all,file="output/model/mod_optimistic_all.rds")

# mod_pessimistic_all = sampling(comp_model_pessimistic,
#                                data = stan_d,
#                                pars = params_monitor,
#                                chains = 2,
#                                seed = seed_pessimistic,
#                                warmup = 300,
#                                init = "random",
#                                iter = 2000)
# saveRDS(mod_pessimistic_all,file="output/model/mod_pessimistic_all.rds")

mod_optimistic_all <- readRDS(file="output/model/mod_optimistic_all.rds")
mod_pessimistic_all <- readRDS(file="output/model/mod_pessimistic_all.rds")

##################################################################################
# Compare models
##################################################################################
log_lik_optimistic <- extract_log_lik(mod_optimistic_all, merge_chains = FALSE)
log_lik_pessimistic <- extract_log_lik(mod_pessimistic_all, merge_chains = FALSE)

r_eff_optimistic <- relative_eff(exp(log_lik_optimistic))
r_eff_pessimistic <- relative_eff(exp(log_lik_pessimistic)) 

loo_optimistic <- loo(log_lik_optimistic, r_eff = r_eff_optimistic, cores = 4)
loo_pessimistic <- loo(log_lik_pessimistic, r_eff = r_eff_pessimistic, cores = 4)
print(loo_optimistic)
print(loo_pessimistic)

comp <- loo_compare(loo_optimistic, loo_pessimistic)
print(comp)

##################################################################################
# extract simulated values
##################################################################################
# fake_I.post <- extract(mod_6_new)$fake_I

mod_optimistic_posterior <- data.frame(beta=rstan::extract(mod_optimistic_all,inc_warmup=FALSE)$beta,
                                       S0=rstan::extract(mod_optimistic_all, 
                                                         inc_warmup=FALSE)$S0,
                                       theta1=rstan::extract(mod_optimistic_all,
                                                             inc_warmup=FALSE)$theta1,
                                       theta2=rstan::extract(mod_optimistic_all,
                                                             inc_warmup=FALSE)$theta2,
                                       theta3=rstan::extract(mod_optimistic_all,
                                                             inc_warmup=FALSE)$theta3,
                                       p1=rstan::extract(mod_optimistic_all, inc_warmup=FALSE)$p1,
                                       p2=rstan::extract(mod_optimistic_all, inc_warmup=FALSE)$p2,
                                       kappa=rstan::extract(mod_optimistic_all, inc_warmup=FALSE)$kappa,
                                       sigmaE_1=rstan::extract(mod_optimistic_all,
                                                               inc_warmup=FALSE)$sigmaE_1,
                                       sigmaE_2=rstan::extract(mod_optimistic_all,
                                                               inc_warmup=FALSE)$sigmaE_2,
                                       R0=rstan::extract(mod_optimistic_all,inc_warmup=FALSE)$R0)

mod_optimistic_summary <- data.frame(param=colnames(mod_optimistic_posterior),
                                     mean=as.numeric(apply(mod_optimistic_posterior,2,mean)),
                                     median=as.numeric(apply(mod_optimistic_posterior,2,quantile,
                                                             probs=0.5)),
                                     lower=as.numeric(apply(mod_optimistic_posterior,2,quantile,
                                                            probs=0.025)),
                                     upper=as.numeric(apply(mod_optimistic_posterior,2,quantile,
                                                            probs=0.975)))

mod_pessimistic_posterior <- data.frame(beta=rstan::extract(mod_pessimistic_all,inc_warmup=FALSE)$beta,
                                        S0=rstan::extract(mod_pessimistic_all, 
                                                          inc_warmup=FALSE)$S0,
                                        theta1=rstan::extract(mod_pessimistic_all,
                                                              inc_warmup=FALSE)$theta1,
                                        theta2=rstan::extract(mod_pessimistic_all,
                                                              inc_warmup=FALSE)$theta2,
                                        theta3=rstan::extract(mod_pessimistic_all,
                                                              inc_warmup=FALSE)$theta3,
                                        p1=rstan::extract(mod_pessimistic_all,inc_warmup=FALSE)$p1,
                                        p2=rstan::extract(mod_pessimistic_all, inc_warmup=FALSE)$p2,
                                        kappa=rstan::extract(mod_pessimistic_all,
                                                             inc_warmup=FALSE)$kappa,
                                        sigmaE_1=rstan::extract(mod_pessimistic_all,
                                                                inc_warmup=FALSE)$sigmaE_1,
                                        sigmaE_2=rstan::extract(mod_pessimistic_all,
                                                                inc_warmup=FALSE)$sigmaE_2,
                                        R0=rstan::extract(mod_pessimistic_all,inc_warmup=FALSE)$R0)

mod_pessimistic_summary <- data.frame(param=colnames(mod_pessimistic_posterior),
                                      mean=as.numeric(apply(mod_pessimistic_posterior,2,mean)),
                                      median=as.numeric(apply(mod_pessimistic_posterior,2,quantile,
                                                              probs=0.5)),
                                      lower=as.numeric(apply(mod_pessimistic_posterior,2,quantile,
                                                             probs=0.025)),
                                      upper=as.numeric(apply(mod_pessimistic_posterior,2,quantile,
                                                             probs=0.975)))

mod_optimistic_inci_samp <- rstan::extract(mod_optimistic_all, inc_warmup=FALSE)$inci_samp
mod_pessimistic_inci_samp <- rstan::extract(mod_pessimistic_all, inc_warmup=FALSE)$inci_samp

mod_optimistic_median_posterior <- colMedians(as.matrix(mod_optimistic_posterior))
names(mod_optimistic_median_posterior) <- colnames(mod_optimistic_posterior)
mod_pessimistic_median_posterior <- colMedians(as.matrix(mod_pessimistic_posterior))
names(mod_pessimistic_median_posterior) <- colnames(mod_pessimistic_posterior)
mod_optimistic_quantiles_inci_samp <- colQuantiles(as.matrix(mod_optimistic_inci_samp),
                                                   probs=c(0.025,0.5,0.975))
mod_pessimistic_quantiles_inci_samp <- colQuantiles(as.matrix(mod_pessimistic_inci_samp),
                                                    probs=c(0.025,0.5,0.975))

# plot estim and true data
# optimistic model
plot(1:18,cases_jakarta_week$cases,col="black",cex=1.5,
     ylim=c(0,max(mod_optimistic_quantiles_inci_samp)),pch=16)
lines(1:18,mod_optimistic_quantiles_inci_samp[,2],lwd=2,col="red")
lines(1:18,mod_optimistic_quantiles_inci_samp[,1],lwd=2,lty=2,col="red")
lines(1:18,mod_optimistic_quantiles_inci_samp[,3],lwd=2,lty=2,col="red")
# pessimistic model
plot(1:18,cases_jakarta_week$cases,col="black",cex=1.5,
     ylim=c(0,max(mod_pessimistic_quantiles_inci_samp)),pch=16)
lines(1:18,mod_pessimistic_quantiles_inci_samp[,2],lwd=2,col="red")
lines(1:18,mod_pessimistic_quantiles_inci_samp[,1],lwd=2,lty=2,col="red")
lines(1:18,mod_pessimistic_quantiles_inci_samp[,3],lwd=2,lty=2,col="red")

##################################################################################
# get traceplot and parameter estimates from model
##################################################################################
## traceplots
p_trace_optimistic <- stan_trace(mod_optimistic_all, pars = c("beta","theta1","theta2",
                                                              "theta3","p1","p2","kappa",
                                                              "sigmaE_1","sigmaE_2","S0","R0"),
                                 inc_warmup = TRUE,ncol=3) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p_trace_pessimistic <- stan_trace(mod_pessimistic_all, pars = c("beta","theta1","theta2",
                                                                "theta3","p1","p2","kappa",
                                                                "sigmaE_1","sigmaE_2","S0","R0"),
                                  inc_warmup = TRUE,ncol=3) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p_trace_all <- plot_grid(p_trace_optimistic,p_trace_pessimistic,
                         labels=c("A","B"),ncol=1)

ggsave("output/figures/p_trace.png", p_trace_all,
       width = 16, height = 26, units = "cm")

## parameter estimates
# simulate having posteriors for two different models each with parameters beta[1],..., beta[4]
posterior_1 <- matrix(rnorm(4000), 1000, 4)
posterior_2 <- matrix(rnorm(4000), 1000, 4)
colnames(posterior_1) <- colnames(posterior_2) <- paste0("beta[", 1:4, "]")

# use bayesplot::mcmc_intervals_data() function to get intervals data in format easy to pass to ggplot
combined_posteriors <- rbind(mcmc_intervals_data(mod_optimistic_all,prob_outer = 0.95),
                             mcmc_intervals_data(mod_pessimistic_all,prob_outer = 0.95)) %>% 
  filter(parameter %in% c("beta","theta1","theta2",
                          "theta3","p1","p2","kappa",
                          "sigmaE_1","sigmaE_2","S0","R0")) %>% 
  mutate(parameter=factor(parameter,levels=c("theta1","theta2",
                                             "theta3","p1","p2","kappa",
                                             "sigmaE_1","sigmaE_2","beta","S0","R0")))
combined_posteriors$model <- rep(c("Optimistic", "Pessimistic"), 
                                 each = length(c("beta","theta1","theta2",
                                                 "theta3","p1","p2","kappa",
                                                 "sigmaE_1","sigmaE_2","S0","R0")))

# make the plot using ggplot 
pos <- position_nudge(y = ifelse(combined_posteriors$model == "Pessimistic", -0.2, 0.2))
p_estimates <- ggplot(combined_posteriors, aes(x = m, y = parameter, color = model)) + 
  geom_point(position = pos, size=2) +
  geom_linerange(aes(xmin = ll, xmax = hh), position = pos) +
  geom_linerange(aes(xmin = l, xmax = h), position = pos, linewidth=1.25) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(col="Model",x="Parameter estimates",y="Model parameters") +
  theme(legend.position = "top") +
  scale_colour_manual(breaks=c("Optimistic","Pessimistic"),
                      labels=c("Optimistic","Pessimistic"),
                      values=safe_colorblind_palette[c(7,9)])

ggsave("output/figures/p_estimates.png", p_estimates,
       width = 16, height = 12, units = "cm")

# tables of parameter estimates

params_summary_optimistic <- rownames_to_column(data.frame(summary(mod_optimistic_all,
                                                       params_monitor[1:11])$summary))
params_summary_pessimistic <- rownames_to_column(data.frame(summary(mod_pessimistic_all,
                                                        params_monitor[1:11])$summary))

write_csv(params_summary_optimistic,"output/model/parameter_estimates_optimistic.csv")
write_csv(params_summary_pessimistic,"output/model/parameter_estimates_pessimistic.csv")







