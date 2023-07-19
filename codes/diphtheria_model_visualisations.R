library(fitR)
library(gridGraphics)
source("codes/diphtheria_model_functions.R")

### load data
mod_optimistic_all <- readRDS(file="output/model/mod_optimistic_all.rds")
mod_pessimistic_all <- readRDS(file="output/model/mod_pessimistic_all.rds")

### run deterministic model based on parameter estimates
# seed: 981054 1009482
# process_model_optimistic <- process_model_parameters_deterministic(mod_optimistic_median_posterior,
#                                                                    mod_optimistic_posterior,
#                                                                    cases_jakarta_week$date_of_the_year,
#                                                                    981054)
# process_model_pessimistic <- process_model_parameters_deterministic(mod_pessimistic_median_posterior,
#                                                                     mod_pessimistic_posterior,
#                                                                     cases_jakarta_week$date_of_the_year,
#                                                                     1009482)
# 
# saveRDS(process_model_optimistic,"output/model/process_model_optimistic.rds")
# saveRDS(process_model_pessimistic,"output/model/process_model_pessimistic.rds")

process_model_optimistic <- readRDS("output/model/process_model_optimistic.rds")
process_model_pessimistic <- readRDS("output/model/process_model_pessimistic.rds")

times_of_observed <- cases_jakarta_week$date_of_the_year
incidence_observed <- cases_jakarta_week$cases

p3 <- process_model_deterministic_plot(process_model_optimistic,
                                       times_of_observed,
                                       incidence_observed)
ggsave("output/figures/p3.jpg", ggdraw(p3),
       width = 16, height = 11, units = "cm")

p_pessimistic <- process_model_deterministic_plot(process_model_pessimistic,
                                                  times_of_observed,
                                                  incidence_observed)
# p_all_model <- plot_grid(ggdraw(p3),ggdraw(p_pessimistic),ncol=1,
#                          labels=c("A","B"), align = "v")
ggsave("output/figures/p_pessimistic.jpg", ggdraw(p_pessimistic),
       width = 16, height = 11, units = "cm")

### run scenarios
# output_model_without_ORI <- calculate_diff(mod_optimistic_posterior,model_without_ORI,325663)
# output_model_without_CT <- calculate_diff(mod_optimistic_posterior,model_without_CT,173401)
# output_model_delayed_ORI <- calculate_diff(mod_optimistic_posterior,model_delayed_ORI,173674)
# 
# saveRDS(output_model_without_ORI,"output/model/output_model_without_ORI.rds")
# saveRDS(output_model_without_CT,"output/model/output_model_without_CT.rds")
# saveRDS(output_model_delayed_ORI,"output/model/output_model_delayed_ORI.rds")

output_model_without_ORI <- readRDS("output/model/output_model_without_ORI.rds")
output_model_without_CT <- readRDS("output/model/output_model_without_CT.rds")
output_model_delayed_ORI <- readRDS("output/model/output_model_delayed_ORI.rds")

### estimate different intervention impacts
diff.long.all <- data.frame(Scenario=rep(c("Scenario 1","Scenario 2","Scenario 3"),each=5000),
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
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
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


