library(fitR)
library(gridGraphics)
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
source("codes/diphtheria_model_functions.R")

### load data
mod1_all <- readRDS(file="output/model/mod1_all.rds")
mod2_all <- readRDS(file="output/model/mod2_all.rds")

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

cases_jakarta_week <- cases_admin2_week %>% 
  group_by(date) %>% 
  summarise(cases=sum(cases)) %>% 
  ungroup() %>% 
  mutate(date_of_the_year=yday(date)) %>% # next line convert timeline of 2018 data to cont. from 2017
  mutate(date_of_the_year=ifelse(date_of_the_year<300,date_of_the_year+365,date_of_the_year))

### run deterministic model based on parameter estimates
# seed: 981054 1009482
# process_model1 <- process_model_parameters_deterministic(
#   mod1_median_posterior,mod1_posterior,
#   cases_jakarta_week$date_of_the_year,981054)
# saveRDS(process_model1,"output/model/process_model1.rds")
# 
# process_model2 <- process_model_parameters_deterministic(
#   mod2_median_posterior,mod2_posterior,
#   cases_jakarta_week$date_of_the_year,1009482)
# saveRDS(process_model2,"output/model/process_model2.rds")

process_model1 <- readRDS("output/model/process_model1.rds")
process_model2 <- readRDS("output/model/process_model2.rds")

times_of_observed <- cases_jakarta_week$date_of_the_year
incidence_observed <- cases_jakarta_week$cases

p3 <- process_model_deterministic_plot(process_model1,
                                       times_of_observed,
                                       incidence_observed)
ggsave("output/figures/p3.jpg", ggdraw(p3),
       width = 16, height = 11, units = "cm")

p_mod2 <- process_model_deterministic_plot(process_model2,
                                           times_of_observed,
                                           incidence_observed)

# p_all_model <- plot_grid(ggdraw(p3),ggdraw(p_sens_low),ncol=1,
#                          labels=c("A","B"), align = "v")
ggsave("output/figures/p_mod2.jpg", ggdraw(p_mod2),
       width = 16, height = 11, units = "cm")

### run scenarios
output_model_without_ORI <- calculate_diff(mod1_posterior,model_without_ORI,325663)
output_model_without_CT <- calculate_diff(mod1_posterior,model_without_CT,173401)
output_model_delayed_ORI <- calculate_diff(mod1_posterior,model_delayed_ORI,173674)

saveRDS(output_model_without_ORI,"output/model/output_model_without_ORI.rds")
saveRDS(output_model_without_CT,"output/model/output_model_without_CT.rds")
saveRDS(output_model_delayed_ORI,"output/model/output_model_delayed_ORI.rds")

output_model_without_ORI <- readRDS("output/model/output_model_without_ORI.rds")
output_model_without_CT <- readRDS("output/model/output_model_without_CT.rds")
output_model_delayed_ORI <- readRDS("output/model/output_model_delayed_ORI.rds")

### estimate different intervention impacts
diff.long.all <- data.frame(Scenario=rep(c("Scenario 1","Scenario 2","Scenario 3"),
                                         each=5000),
                            Difference=c(output_model_without_ORI$diff.intv.all,
                                         output_model_without_CT$diff.intv.all,
                                         output_model_delayed_ORI$diff.intv.all))
diff.long <- data.frame(Scenario=rep(c("Scenario 1","Scenario 2","Scenario 3"),each=5000),
                        Difference=c(output_model_without_ORI$diff.intv,
                                     output_model_without_CT$diff.intv,
                                     output_model_delayed_ORI$diff.intv))
diff.total.all <- data.frame(Scenario=c("Scenario 1","Scenario 2","Scenario 3"),
                             Difference=rbind(output_model_without_ORI$diff.intv.all.summary[[1]],
                                              output_model_without_CT$diff.intv.all.summary[[1]],
                                              output_model_delayed_ORI$diff.intv.all.summary[[1]]))
diff.total <- data.frame(Scenario=c("Scenario 1","Scenario 2","Scenario 3"),
                         Difference=rbind(output_model_without_ORI$diff.intv.summary[[1]],
                                          output_model_without_CT$diff.intv.summary[[1]],
                                          output_model_delayed_ORI$diff.intv.summary[[1]]))

p5 <- diff.long.all %>% ggplot(aes(x=Scenario,y=Difference,fill=Scenario)) +
  geom_violin() + scale_y_continuous(trans="log10",
                                     breaks=c(0,10,100,1000,10000,100000,1000000),
                                     labels=c("0","10","100","1,000","10,000","100,000","1,000,000")) +
  theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_point(data=diff.total.all,aes(x=Scenario,y=Difference),col="red") +
  geom_text(data=diff.total.all,aes(x=Scenario,y=0.75*Difference,
                                    label=scales::label_comma(accuracy = NULL)(round(Difference,0)))) +
  labs(x=NULL,y="Estimated total clinical cases averted\nfrom performing the actual interventions\n(ORI & contact tracing)") +
  scale_fill_manual(breaks=c("Scenario 1","Scenario 2","Scenario 3"),
                    labels=c("Scenario 1","Scenario 2","Scenario 3"),
                    values=safe_colorblind_palette[c(1,7,3)]) +
  theme(legend.position = "none")
ggsave("output/figures/p5.png", p5,
       width = 16, height = 11, units = "cm")


