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
# params <- list(gammaI_1 = 5.88,
#                gammaI_2 = 3.12,
#                gammaC1 = 0,
#                gammaC2 = 0,
#                theta1 = 0,
#                theta2 = 0,
#                beta = 0.6,
#                p1 = 0,
#                tau = 0,
#                p = 0,
#                px = 0,
#                rvac = 0,
#                S0 = 0.25,
#                t0 = 0)
# 
params <- list(beta = 0,
               theta1 = 0,
               theta2 = 0,
               theta3 = 0,
               p1 = 0,
               p2 = 0,
               kappa = 0,
               sigmaE_1 = 0,
               sigmaE_2 = 0,
               sigmaC_1 = 0,
               sigmaC_2 = 0)

param_init <- list(beta = 0.4,
                   theta1 = 0.1,
                   theta2 = 0.05,
                   theta3 = 0.9,
                   p1 = 0.5,
                   p2 = 0.8,
                   kappa = 0.7,
                   sigmaE_1 = 0.15,
                   sigmaE_2 = 0.15,
                   sigmaC_1 = 0.15,
                   sigmaC_2 = 0.15,
                   S0=0.4)

##################################################################################
# data from individual-data-viz.R code
# look at the final section of the code 'Read Individual'
# the source is the individual data collected by Jakarta Health Agency
# during the outbreak
# first case Jakarta Selatan 1-Nov = week zero
# week 1 = week 45 - 5-11 November 2017
# we fit model for all Jakarta for now; data.weekly.indiv$All
# set time relative to 1 Jan 2017; first week ended at 7 Jan 2017
# first vaccination campaign started at 10 Dec 2017 = day 343
#
#
# Specifically for this data:
# make sure we are far enough into the season that there is at least one case per week
# no cases at week 3 -> assumed to have 0 case, as it's not possible to model 0 cases
##################################################################################

data_district_week <- read.csv("data/difteri-individual-district-calendarweek-final.csv",header = T)
data_weekly_indiv <- data.frame(Jakarta_Pusat=table(data_district_week$CalendarWeek,
                                                    data_district_week$District)[,2],
                                Jakarta_Utara=table(data_district_week$CalendarWeek,
                                                    data_district_week$District)[,5],
                                Jakarta_Barat=table(data_district_week$CalendarWeek,
                                                    data_district_week$District)[,1],
                                Jakarta_Selatan=table(data_district_week$CalendarWeek,
                                                      data_district_week$District)[,3],
                                Jakarta_Timur=table(data_district_week$CalendarWeek,
                                                    data_district_week$District)[,4])

data_weekly_indiv$All <- rowSums(data_weekly_indiv)

incidence_observed = data_weekly_indiv$All
incidence_observed = c(incidence_observed[1:3],0,incidence_observed[4:17])
week44_rel_to_jan_1 = julian(1,7,2017)+(44-1)*7-julian(1,1,2017)+1
times_of_observed = seq(from=week44_rel_to_jan_1, by=7, length.out=18)
times_of_observed2 = seq(from=week44_rel_to_jan_1, by=7, length.out=35)
times_of_observed3 = seq(from=week44_rel_to_jan_1, by=7, length.out=20)
time_binning = min(diff(times_of_observed))

diphdata = data.frame(time=times_of_observed, obs=incidence_observed)

##################################################################################
# set samples timing
##################################################################################
fake_ts = 301:426

##################################################################################
# set rstan input
##################################################################################
stan_d = list(n_obs = length(times_of_observed),
              n_params = 25,
              n_difeq = length(inits),
              n_fake = length(fake_ts),
              y = incidence_observed,
              N = npop,
              t0 = 300,
              ts = times_of_observed,
              fake_ts = fake_ts,
              n_samples = length(times_of_observed),
              init=param_init)

##################################################################################
# set parameters to monitor
##################################################################################
params_monitor = c("beta","theta1","theta2","theta3","p1","p2","kappa","sigmaE_1","sigmaE_2","sigmaC_1","sigmaC_2","S0","R0","inci_samp","log_lik")

##################################################################################
# debug + run the model
##################################################################################
comp_model_optimistic <- stan_model('codes/diphtheria-model-optimistic.stan')
comp_model_pessimistic <- stan_model('codes/diphtheria-model-pessimistic.stan')

# runif(2,1,9999999); 6813636 6592353
seed_optimistic <- 6813636
seed_pessimistic <- 6592353

# mod_optimistic_c1 = sampling(comp_model_optimistic,
#                              data = stan_d,
#                              pars = params_monitor,
#                              chains = 1,
#                              seed = seed_optimistic,
#                              chain_id = 1,
#                              warmup = 300,
#                              init = "random",
#                              iter = 2000)
# mod_optimistic_c2 = sampling(comp_model_optimistic,
#                              data = stan_d,
#                              pars = params_monitor,
#                              chains = 1,
#                              seed = seed_optimistic,
#                              chain_id = 2,
#                              warmup = 300,
#                              init = "random",
#                              iter = 2000)
# mod_optimistic_all <- sflist2stanfit(list(mod_optimistic_c1, mod_optimistic_c2))
# saveRDS(mod_optimistic_all,file="output/mod_optimistic_all.rds")
# 
# mod_pessimistic_c1 = sampling(comp_model_pessimistic,
#                               data = stan_d,
#                               pars = params_monitor,
#                               chains = 1,
#                               seed = seed_pessimistic,
#                               chain_id = 1,
#                               warmup = 300,
#                               init = "random",
#                               iter = 2000)
# mod_pessimistic_c2 = sampling(comp_model_pessimistic,
#                               data = stan_d,
#                               pars = params_monitor,
#                               chains = 1,
#                               seed = seed_pessimistic,
#                               chain_id = 2,
#                               warmup = 300,
#                               init = "random",
#                               iter = 2000)
# mod_pessimistic_all <- sflist2stanfit(list(mod_pessimistic_c1, mod_pessimistic_c2))
# saveRDS(mod_pessimistic_all,file="output/mod_pessimistic_all.rds")

mod_optimistic_all <- readRDS(file="output/mod_optimistic_all.rds")
mod_pessimistic_all <- readRDS(file="output/mod_pessimistic_all.rds")

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
                                       sigmaC_1=rstan::extract(mod_optimistic_all,
                                                               inc_warmup=FALSE)$sigmaC_1,
                                       sigmaC_2=rstan::extract(mod_optimistic_all,
                                                               inc_warmup=FALSE)$sigmaC_2,
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
                                        sigmaC_1=rstan::extract(mod_pessimistic_all,
                                                                inc_warmup=FALSE)$sigmaC_1,
                                        sigmaC_2=rstan::extract(mod_pessimistic_all,
                                                                inc_warmup=FALSE)$sigmaC_2,
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
plot(1:18,incidence_observed,col="black",cex=1.5,
     ylim=c(0,max(mod_optimistic_quantiles_inci_samp)),pch=16)
lines(1:18,mod_optimistic_quantiles_inci_samp[,2],lwd=2,col="red")
lines(1:18,mod_optimistic_quantiles_inci_samp[,1],lwd=2,lty=2,col="red")
lines(1:18,mod_optimistic_quantiles_inci_samp[,3],lwd=2,lty=2,col="red")
# pessimistic model
plot(1:18,incidence_observed,col="black",cex=1.5,
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
                                                              "sigmaE_1","sigmaE_2",
                                                              "sigmaC_1","sigmaC_2","S0","R0"),
                                 inc_warmup = TRUE,ncol=3) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p_trace_pessimistic <- stan_trace(mod_pessimistic_all, pars = c("beta","theta1","theta2",
                                                                "theta3","p1","p2","kappa",
                                                                "sigmaE_1","sigmaE_2",
                                                                "sigmaC_1","sigmaC_2","S0","R0"),
                                  inc_warmup = TRUE,ncol=3) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p_trace_all <- plot_grid(p_trace_optimistic,p_trace_pessimistic,
                         labels=c("A","B"),ncol=1)

ggsave("figures/new/p_trace.png", p_trace_all,
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
                          "sigmaE_1","sigmaE_2",
                          "sigmaC_1","sigmaC_2","S0","R0")) %>% 
  mutate(parameter=factor(parameter,levels=c("theta1","theta2",
                                             "theta3","p1","p2","kappa",
                                             "sigmaE_1","sigmaE_2",
                                             "sigmaC_1","sigmaC_2","beta","S0","R0")))
combined_posteriors$model <- rep(c("Optimistic", "Pessimistic"), 
                                 each = length(c("beta","theta1","theta2",
                                                 "theta3","p1","p2","kappa",
                                                 "sigmaE_1","sigmaE_2",
                                                 "sigmaC_1","sigmaC_2","S0","R0")))

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

ggsave("figures/new/p_estimates.png", p_estimates,
       width = 16, height = 12, units = "cm")







