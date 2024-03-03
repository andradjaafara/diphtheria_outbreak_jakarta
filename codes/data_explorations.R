### library
library(sf)
library(maptools)
library(incidence)
library(scales)
library(reshape2)
library(gridExtra)
library(ggplot2)
library(gtable)
library(grid)
library(tidyverse)
library(cowplot)
library(lubridate)
library(fitdistrplus)

### colourblind safe palette
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

admin2_jakarta <- c("Central Jakarta","East Jakarta",
                    "North Jakarta","South Jakarta","West Jakarta")

### read and prepare data
### data from the first week of Nov 17 to first week of Mar 17
### a total of 250 clinical diphtheria cases recorded in Jakarta within this period
### the week start on Monday - and the weekly data represent total cases at the end of the week (Sunday)

cases_admin2_week <- read_csv("data/diphteria-cases-weekly-admin2.csv") %>% 
  mutate(date=ymd(date))
pop_admin2 <- read_csv("data/pop-admin2.csv")

### summary of incidence rate per 100,000 by admin2
cases_admin2_summary <- cases_admin2_week %>% 
  group_by(admin2) %>% 
  summarise(cases=sum(cases)) %>% ungroup() %>% 
  left_join(pop_admin2) %>% 
  mutate(incidence_rate=cases/pop*100000)

### plot of incidence rate per 100,000 by admin2
p1_A <- ggplot(data=cases_admin2_summary, aes(x=admin2,y=incidence_rate,fill=admin2)) + 
  geom_bar(stat="identity") +
  geom_text(aes(label=round(incidence_rate,3)), vjust=-0.3, size=3.5) + 
  scale_fill_manual(breaks=admin2_jakarta,
                    labels=admin2_jakarta,
                    values=safe_colorblind_palette[1:5]) +
  xlab("District") + ylab("Incidence rate per 100,000") +
  theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     legend.position = "none") +
  ylim(0,4.15) + labs(tag="A")

### weekly case lines by admin2
### plot of weekly case lines by admin2
p1_B <- ggplot(data=cases_admin2_week, aes(x=date, y=cases, group=admin2)) +
  geom_line(aes(colour=admin2), linewidth=1.5) +
  geom_point(aes(colour=admin2), size=2.5) + 
  scale_colour_manual(breaks=admin2_jakarta,
                      labels=admin2_jakarta,
                      values=safe_colorblind_palette[1:5]) +
  xlab("") + ylab("Weekly incidence") +
  scale_x_date(date_labels="%b %y") +
  theme_classic() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     legend.title=element_blank()) +
  scale_y_continuous(breaks=c(0,5,10,15,20),limits=c(0,20)) +
  theme(legend.position = c(0.875, 0.65),
        legend.background = element_rect(fill='transparent')) + labs(tag="B")

### cases by age
### age groups
age_groups_fitting <- read_csv('data/age_groups.csv')
age_groups_vec <- unique(age_groups_fitting$age_new)

pop_age <- read_csv("data/pop-age.csv") %>% 
  mutate(age_group=factor(age_group,levels=c("0-4","5-9","10-14",
                                             "15-19","20-24","25-29",
                                             "30-34","35-39","40-44", 
                                             "45-49","50-54","55-59",">=60")))

pop_age2 <- pop_age %>% 
  mutate(age_group_new=age_groups_fitting$age_new[1:13]) %>% 
  mutate(age_group_new=factor(age_group_new,levels=age_groups_vec)) %>% 
  group_by(age_group_new) %>% 
  summarise(pop=sum(pop)) %>% ungroup() %>% 
  dplyr::select(age_group=age_group_new,pop)

cases_age <- read_csv("data/diphtheria-cases-age.csv") %>% 
  mutate(age_group=factor(age_group,levels=c("0-4","5-9","10-14",
                                             "15-19","20-24","25-29",
                                             "30-34","35-39","40-44", 
                                             "45-49","50-54","55-59",">=60"))) %>% 
  left_join(pop_age) %>% 
  mutate(age_group_new=age_groups_fitting$age_new[1:13]) %>% 
  mutate(age_group_new=factor(age_group_new,levels=age_groups_vec)) %>% 
  group_by(age_group_new) %>% 
  summarise(pop=sum(pop),cases=sum(cases)) %>% ungroup() %>% 
  mutate(incidence_rate=cases/pop*100000,
         prop_cases=cases/sum(cases)*100) %>% 
  dplyr::select(age_group=age_group_new,pop,cases,incidence_rate,prop_cases)
  
## create cases by age based on the new age groups
## then for supplementary create for each phases

# write_csv(cases_age,"output/csv/cases_by_age.csv")

p1_C <- ggplot(data=cases_age, aes(x=age_group,y=incidence_rate)) + 
  geom_bar(stat="identity",fill=safe_colorblind_palette[11]) +
  geom_text(aes(label=round(incidence_rate,3)), vjust=-0.3, size=3.5) + 
  xlab("Age group") + ylab("Incidence rate per 100,000") + 
  labs(fill=NULL) +
  theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks=c(0,2,4,6,8),limits=c(0,8)) +
  theme(legend.position = c(0.9, 0.80)) + labs(tag="C")

p1 <- plot_grid(p1_A,p1_B,p1_C,nrow=3)

ggsave("output/figures/p1.png", p1,
       width = 16, height = 20, units = "cm")

### susceptibility estimates
# susceptibility_estimates <- read.csv("data/susceptibility-estimates.csv",header = TRUE,
#                                      stringsAsFactors = FALSE) %>% 
#   mutate(age_group=factor(age_group,levels=c("0-4","5-9","10-14","15-19","20-24",
#                                              "25-29","30-34","35-39","40-44",
#                                              "45-49","50-54","55-59",">=60","Overall")))
# 
# p2 <- susceptibility_estimates %>% 
#   ggplot(aes(x=age_group,y=susceptible_proportion*100,fill=scenario)) +
#   geom_bar(stat="identity", position=position_dodge()) +
#   scale_fill_manual(breaks=c("Optimistic","Pessimistic"),
#                     labels=c("Optimistic","Pessimistic"),
#                     values=safe_colorblind_palette[c(7,9)]) + 
#   labs(fill=NULL,y="Population susceptibility (%)",x="Age group") +
#   theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   scale_y_continuous(breaks=seq(0,50,10),limits=c(0,50)) +
#   theme(legend.position = c(0.9, 0.825))
# 
# ggsave("output/figures/p2.png", p2,
#        width = 16, height = 7, units = "cm")