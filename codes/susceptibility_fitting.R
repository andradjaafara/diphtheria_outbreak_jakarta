### Load rstan library
library(rstan)

### Data preparations
contact_matrix_long <- read_csv('data/contact_matrix.csv') # IDN urban from prem et al.
age_groups_fitting <- read_csv('data/age_groups.csv')
age_groups_vec <- unique(age_groups_fitting$age_new)
phase_vec <- c("Phase 1","Phase 2","Phase 3","Phase 4")

### 20 to 59 contains 8 groups
### 60+ contains 4 groups

### Create new contact matrix (wide) with new age groups
contact_matrix <- contact_matrix_long %>% dplyr::select(age_contactor,age_cotactee,
                                                        mean_number_of_contacts) %>% 
  left_join(age_groups_fitting,by=join_by(age_contactor == age_old)) %>% 
  rename(age_contactor_new=age_new) %>% 
  left_join(age_groups_fitting,by=join_by(age_cotactee == age_old)) %>% 
  rename(age_contactee_new=age_new) %>% 
  mutate(age_contactor_new=factor(age_contactor_new,levels=age_groups_vec),
         age_contactee_new=factor(age_contactee_new,levels=age_groups_vec)) %>% 
  group_by(age_contactor_new,age_contactee_new) %>% 
  summarise(mean_number_of_contacts=sum(mean_number_of_contacts)) %>% 
  ungroup %>% 
  pivot_wider(names_from=age_contactee_new,values_from=mean_number_of_contacts)

contact_matrix_count <- contact_matrix_long %>% dplyr::select(age_contactor,age_cotactee,
                                                        mean_number_of_contacts) %>% 
  left_join(age_groups_fitting,by=join_by(age_contactor == age_old)) %>% 
  rename(age_contactor_new=age_new) %>% 
  left_join(age_groups_fitting,by=join_by(age_cotactee == age_old)) %>% 
  rename(age_contactee_new=age_new) %>% 
  mutate(age_contactor_new=factor(age_contactor_new,levels=age_groups_vec),
         age_contactee_new=factor(age_contactee_new,levels=age_groups_vec)) %>% 
  group_by(age_contactor_new,age_contactee_new) %>% 
  summarise(count_number_of_contacts=n()) %>% 
  ungroup %>% 
  pivot_wider(names_from=age_contactee_new,values_from=count_number_of_contacts)

### need element wise averaging for 20 to 59 and 60+ rows and columns
### averaging is needed as 20 to 59  and 60+ are comprised multiple age groups initially
contact_matrix[5,2] <- contact_matrix[5,2]/8
contact_matrix[6,2] <- contact_matrix[6,2]/4
contact_matrix[5,3] <- contact_matrix[5,3]/8
contact_matrix[6,3] <- contact_matrix[6,3]/4
contact_matrix[5,4] <- contact_matrix[5,4]/8
contact_matrix[6,4] <- contact_matrix[6,4]/4
contact_matrix[5,5] <- contact_matrix[5,5]/8
contact_matrix[6,5] <- contact_matrix[6,5]/4
contact_matrix[5,6] <- contact_matrix[5,6]/64
contact_matrix[6,6] <- contact_matrix[6,6]/32
contact_matrix[5,7] <- contact_matrix[5,7]/32
contact_matrix[6,7] <- contact_matrix[6,7]/16
contact_matrix[1,6] <- contact_matrix[1,6]/8
contact_matrix[1,7] <- contact_matrix[1,7]/4
contact_matrix[2,6] <- contact_matrix[2,6]/8
contact_matrix[2,7] <- contact_matrix[2,7]/4
contact_matrix[3,6] <- contact_matrix[3,6]/8
contact_matrix[3,7] <- contact_matrix[3,7]/4
contact_matrix[4,6] <- contact_matrix[4,6]/8
contact_matrix[4,7] <- contact_matrix[4,7]/4

### now as a matrix
cm <- as.matrix(contact_matrix[,-1])

### Set the number of age groups, phases
N <- 6 # 6 age groups: 0-4, 5-9, 10-14, 15-19, 20-59, 60+
P <- 4 # 4 phases of the outbreak in Jakarta

### Can actually use different contact matrices for different phases

### Set contact matrix data for fitting
X <- array(numeric(),c(P,N,N))
X[1,,] <- cm
X[2,,] <- cm
X[3,,] <- cm
X[4,,] <- cm

### Set cases data for fitting
cases_age_phase <- read_csv('data/age_phase.csv') %>% 
  mutate(age_group=factor(age_group,levels=age_groups_vec))

y <- matrix(numeric(), P, N)
y[1,] <- cases_age_phase %>% filter(phase=="Phase 1") %>% pull(count)
y[2,] <- cases_age_phase %>% filter(phase=="Phase 2") %>% pull(count)
y[3,] <- cases_age_phase %>% filter(phase=="Phase 3") %>% pull(count)
y[4,] <- cases_age_phase %>% filter(phase=="Phase 4") %>% pull(count)

z <- y

### Set population data for fitting
pop <- pop_age2$pop * 0.01
pop <- round(pop)

### Package the data
data_list <- list(N = N, P = P, X = X, y = y, z = z, pop = pop)

### Specify the Stan model file
### Model 1 - only youngest prior
### Model 2 - with adult prior
stan_file1 <- 'codes/stan/susceptibility_fitting_new.stan' # no prior for adults
stan_file2 <- 'codes/stan/susceptibility_fitting_adult_prior_new.stan' # with prior for adults

### Compile the model:
### Model 1 - only youngest prior
### Model 2 - with adult prior
susceptibility_model1 <- stan_model(stan_file1)
susceptibility_model2 <- stan_model(stan_file2)

### Fit the model
fit1 <- sampling(susceptibility_model1, data = data_list, iter = 2000, 
                 chains = 4, warmup = 500, seed = 2134654)
fit2 <- sampling(susceptibility_model2, data = data_list, iter = 2000, 
                 chains = 4, warmup = 500, seed = 2134654)

### Print the summary of the fit
summary(fit1)
summary(fit2)

### Extract the generated quantities
generated_quantities1 <- rstan::extract(fit1, "expected_z")
generated_quantities2 <- rstan::extract(fit2, "expected_z")

### Generated quantities model 1
generated_quantities1_list <- list()
for(i in seq_len(length(age_groups_vec))){
  gen_quant_list <- list()
  for (j in seq_len(length(phase_vec))){
    gen_quant_list[[j]] <- tibble(phase=phase_vec[j],
                                  age_group=age_groups_vec[i],
                                  cases_sim=generated_quantities1$expected_z[,,i][,j],
                                  iter=1:6000)
  }
  generated_quantities1_list[[i]] <- bind_rows(gen_quant_list)
}
generated_quantities1_df <- bind_rows(generated_quantities1_list)

generated_quantities1_summary <- generated_quantities1_df %>% 
  group_by(phase,age_group) %>% 
  summarise(median=quantile(cases_sim,probs=0.5),
            lower=quantile(cases_sim,probs=0.025),
            upper=quantile(cases_sim,probs=0.975)) %>% 
  ungroup() %>% 
  mutate(age_group=factor(age_group,levels=age_groups_vec)) %>% 
  mutate(type="Model 1")

### Generated quantities model 2
generated_quantities2_list <- list()
for(i in seq_len(length(age_groups_vec))){
  gen_quant_list <- list()
  for (j in seq_len(length(phase_vec))){
    gen_quant_list[[j]] <- tibble(phase=phase_vec[j],
                                  age_group=age_groups_vec[i],
                                  cases_sim=generated_quantities2$expected_z[,,i][,j],
                                  iter=1:6000)
  }
  generated_quantities2_list[[i]] <- bind_rows(gen_quant_list)
}
generated_quantities2_df <- bind_rows(generated_quantities2_list)

generated_quantities2_summary <- generated_quantities2_df %>% 
  group_by(phase,age_group) %>% 
  summarise(median=quantile(cases_sim,probs=0.5),
            lower=quantile(cases_sim,probs=0.025),
            upper=quantile(cases_sim,probs=0.975)) %>% 
  ungroup() %>% 
  mutate(age_group=factor(age_group,levels=age_groups_vec)) %>% 
  mutate(type="Model 2")

### Combine generated quantities of both models and data
cases_age_phase_plotting <- cases_age_phase %>% 
  mutate(median=count,
         lower=count,
         upper=count,
         type="Data")

generated_quantities_summary <- bind_rows(cases_age_phase_plotting,
                                          generated_quantities1_summary,
                                          generated_quantities2_summary)

p2_A <- generated_quantities_summary %>% 
  ggplot(aes(x=age_group,col=type)) + 
  geom_pointrange(aes(ymin=lower, ymax=upper, y=median),
                  position = position_dodge(width=0.75),size=0.1) +
  facet_wrap(.~phase,nrow=2) + theme_classic() +
  scale_colour_manual(breaks=c("Data","Model 1","Model 2"),
                      values=safe_colorblind_palette[c(1,2,3)]) +
  labs(colour=NULL,y="Cases",x="Age groups") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(tag="A")

p2_legend <- generated_quantities_summary %>% 
  ggplot(aes(x=age_group,col=type)) + 
  geom_pointrange(aes(ymin=lower, ymax=upper, y=median),
                  position = position_dodge(width=0.75),size=0.1) +
  facet_wrap(.~phase,nrow=2) + theme_classic() +
  scale_colour_manual(breaks=c("Data","Model 1","Model 2"),
                      values=safe_colorblind_palette[c(1,2,3)]) +
  labs(colour=NULL,y="Cases",x="Age groups") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(tag="A")
p2_legend <- get_legend(p2_legend)

### Extract parameter 's'
s_samples1 <- rstan::extract(fit1, "s")$s
s_samples2 <- rstan::extract(fit2, "s")$s

### Extract susceptibility estimates - Model 1
susceptibility1_list <- list()
for(i in seq_len(length(age_groups_vec))){
  susceptibility1_list[[i]] <- tibble(age_group=age_groups_vec[i],
                                     susceptibility=s_samples1[,i],
                                     iter=1:6000)
}
susceptibility1_df <- bind_rows(susceptibility1_list)

susceptibility_overall1_df <- susceptibility1_df %>% 
  left_join(pop_age2) %>% 
  group_by(iter) %>% 
  summarise(susceptibility=weighted.mean(susceptibility,pop)) %>% 
  ungroup() %>% 
  mutate(age_group="Overall") %>% 
  dplyr::select(age_group,susceptibility,iter)

susceptibility1_summary <- 
  bind_rows(susceptibility1_df,susceptibility_overall1_df) %>% 
  mutate(age_group=factor(age_group,levels=c(age_groups_vec,"Overall"))) %>% 
  group_by(age_group) %>% 
  summarise(median=quantile(susceptibility,probs=0.5),
            lower=quantile(susceptibility,probs=0.025),
            upper=quantile(susceptibility,probs=0.975)) %>% 
  ungroup() %>% 
  mutate(type="Model 1")

### Extract susceptibility estimates - Model 2
susceptibility2_list <- list()
for(i in seq_len(length(age_groups_vec))){
  susceptibility2_list[[i]] <- tibble(age_group=age_groups_vec[i],
                                      susceptibility=s_samples2[,i],
                                      iter=1:6000)
}
susceptibility2_df <- bind_rows(susceptibility2_list)

susceptibility_overall2_df <- susceptibility2_df %>% 
  left_join(pop_age2) %>% 
  group_by(iter) %>% 
  summarise(susceptibility=weighted.mean(susceptibility,pop)) %>% 
  ungroup() %>% 
  mutate(age_group="Overall") %>% 
  dplyr::select(age_group,susceptibility,iter)

susceptibility2_summary <- 
  bind_rows(susceptibility2_df,susceptibility_overall2_df) %>% 
  mutate(age_group=factor(age_group,levels=c(age_groups_vec,"Overall"))) %>% 
  group_by(age_group) %>% 
  summarise(median=quantile(susceptibility,probs=0.5),
            lower=quantile(susceptibility,probs=0.025),
            upper=quantile(susceptibility,probs=0.975)) %>% 
  ungroup() %>% 
  mutate(type="Model 2")

### Combine susceptibility from Model 1 and Model 2
susceptibility_summary <- bind_rows(susceptibility1_summary,
                                    susceptibility2_summary)

p2_B <- susceptibility_summary %>% 
  ggplot(aes(x=age_group,col=type)) + 
  geom_pointrange(aes(ymin=lower*100, ymax=upper*100, y=median*100),
                  position = position_dodge(width=0.5),size=0.25) +
  labs(y="Susceptibility (%)",x="Age groups",col=NULL) +
  scale_colour_manual(breaks=c("Model 1","Model 2"),
                      values=safe_colorblind_palette[c(2,3)]) +
  theme_classic() +
  geom_hline(yintercept=25,linetype=2,col="grey") +
  geom_hline(yintercept=50,linetype=2,col="grey") +
  geom_hline(yintercept=75,linetype=2,col="grey") + labs(tag="B") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  ylim(0,100)

p2_plot <- plot_grid(p2_A,p2_B,nrow=1)
p2 <- plot_grid(p2_plot,p2_legend,nrow=2,rel_heights = c(1,0.08))

ggsave("output/figures/p2.png", p2,
       width = 16, height = 12, units = "cm")

### Fit overall susceptibility to beta distribution
fit_susceptibility1_beta <- 
  fitdistrplus::fitdist(susceptibility_overall1_df$susceptibility,"beta") ## beta(21.66427,143.11779)
fit_susceptibility2_beta <- 
  fitdistrplus::fitdist(susceptibility_overall2_df$susceptibility,"beta") ## beta(63.27874,106.12358)

## 75% - beta(48,16)
## 25% - beta(6,16)

plot(fit_susceptibility1_beta, las = 1)
plot(fit_susceptibility2_beta, las = 1)

stan_trace(fit1, pars = c("s"),
           inc_warmup = TRUE,ncol=3)
stan_trace(fit2, pars = c("s"),
           inc_warmup = TRUE,ncol=3)
