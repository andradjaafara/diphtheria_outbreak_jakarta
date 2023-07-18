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

### colourblind safe palette
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

admin2_jakarta <- c("Central Jakarta","East Jakarta",
                    "North Jakarta","South Jakarta","West Jakarta")

### read and prepare data
### data from the first week of Nov 17 to first week of Mar 17
### a total of 250 clinical diphtheria cases recorded in Jakarta within this period
data_district_week <- read.csv("data/difteri-individual-district-calendarweek-final-en.csv",header = T)
table_district <- table(data_district_week$District)
population_district <- c(921344,2892783,1781316,2226830,2528065)
incidence_rate_district <- table_district/population_district*100000

jakarta_summary <- data.frame(admin2 = c("Central Jakarta","East Jakarta",
                                         "North Jakarta","South Jakarta","West Jakarta"),
                              incidence=as.vector(table_district), 
                              pop=population_district,
                              incidence_rate=as.vector(incidence_rate_district))

## weekly incidence
inci_date <- as.Date(gsub("/", "-", as.character(data_district_week[,1])),format = "%d-%m-%Y")
inci_district <- as.character(data_district_week$District )

## by district incidence ratio
i_1 <- incidence(inci_date)
i_7 <- incidence(inci_date, interval=7)
i_1_district <- incidence(inci_date, group=inci_district)
i_7_district <- incidence(inci_date, interval=7, group=inci_district)

p1_A <- ggplot(data=jakarta_summary, aes(x=admin2,y=incidence_rate,fill=admin2)) + 
  geom_bar(stat="identity") +
  geom_text(aes(label=round(incidence_rate,3)), vjust=-0.3, size=3.5) + 
  scale_fill_manual(breaks=admin2_jakarta,
                    labels=admin2_jakarta,
                    values=safe_colorblind_palette[1:5]) +
  xlab("District") + ylab("Incidence rate per 100,000") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     legend.position = "none") +
  ylim(0,4.15)

## weekly case lines by district
inci_matrix_district <- data.frame(i_7_district)
inci_matrix_district <- inci_matrix_district %>% dplyr::select(-weeks)
colnames(inci_matrix_district) <- c("dates","isoweeks","Central Jakarta",
                                    "East Jakarta","North Jakarta","South Jakarta",
                                    "West Jakarta")
inci_long_district <- melt(inci_matrix_district[,c(1,3,4,5,6,7)], 
                           id.vars = c("dates"))
inci_long_district$col_line <- rep(incidence_pal1(5),each=18)

p1_B <- ggplot(data=inci_long_district, aes(x=dates, y=value, group=variable)) +
  geom_line(aes(colour=variable), linewidth=1.5) + 
  scale_colour_manual(breaks=admin2_jakarta,
                      labels=admin2_jakarta,
                      values=safe_colorblind_palette[1:5]) +
  xlab("") + ylab("Weekly incidence") +
  scale_x_date(date_labels="%b %y") +
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     legend.title=element_blank()) +
  scale_y_continuous(breaks=c(0,5,10,15,20),limits=c(0,20)) +
  theme(legend.position = c(0.875, 0.655))

## by age and sex cases
age_sex <- read.csv("data/age-sex-data.csv",header = T)
age_sex <- age_sex %>% rename(Age_group=Age.group)
age_sex_inci <- read.csv("data/age-sex-cases.csv",header = T)
age_sex_inci$Age_group <- ""

age_sex$Age_group <- factor(age_sex$Age_group,levels=c("0-4","5-9","10-14",
                                                       "15-19","20-24","25-29",
                                                       "30-34","35-39","40-44", 
                                                       "45-49","50-54","55-59",">=60"))

for (i in seq_len(nrow(age_sex_inci))){
  if (is.na(age_sex_inci[i,1])==T){
    age_sex_inci[i,3] <- ""
  }
  else if (age_sex_inci[i,1]<5){
    age_sex_inci[i,3] <- "0-4"
  }
  else if (age_sex_inci[i,1]<10){
    age_sex_inci[i,3] <- "5-9"
  }
  else if (age_sex_inci[i,1]<15){
    age_sex_inci[i,3] <- "10-14"
  }
  else if (age_sex_inci[i,1]<20){
    age_sex_inci[i,3] <- "15-19"
  }
  else if (age_sex_inci[i,1]<25){
    age_sex_inci[i,3] <- "20-24"
  }
  else if (age_sex_inci[i,1]<30){
    age_sex_inci[i,3] <- "25-29"
  }
  else if (age_sex_inci[i,1]<35){
    age_sex_inci[i,3] <- "30-34"
  }
  else if (age_sex_inci[i,1]<40){
    age_sex_inci[i,3] <- "35-39"
  }
  else if (age_sex_inci[i,1]<45){
    age_sex_inci[i,3] <- "40-44"
  }
  else if (age_sex_inci[i,1]<50){
    age_sex_inci[i,3] <- "45-49"
  }
  else if (age_sex_inci[i,1]<55){
    age_sex_inci[i,3] <- "50-54"
  }
  else if (age_sex_inci[i,1]<60){
    age_sex_inci[i,3] <- "55-59"
  }
  else{
    age_sex_inci[i,3] <- ">=60"
  }
}

age_sex_inci <- age_sex_inci %>% 
  mutate(Age_group=factor(Age_group,
                          levels=c("0-4","5-9","10-14","15-19","20-24",
                                   "25-29","30-34","35-39","40-44",
                                   "45-49","50-54","55-59",">=60")))

inci_age <- as.data.frame(table(age_sex_inci$Age_group))
inci_sex <- as.data.frame(table(age_sex_inci$Sex)[-1]) # remove NA
inci_age_sex <- as.data.frame(table(age_sex_inci$Age_group,age_sex_inci$Sex)[,-1]) # remove NA

colnames(inci_age) <- c("Age_group","Freq")
colnames(inci_sex) <- c("Sex","Freq")
colnames(inci_age_sex) <- c("Age_group","Sex","Freq")

inci_age_sex <- merge(inci_age_sex,age_sex,by=c("Age_group","Sex"))
inci_age_sex$incidence_rate <- inci_age_sex$Freq/inci_age_sex$Pop*100000
inci_age_sex$Age_group <- as.character(inci_age_sex$Age_group)
inci_age_sex$Age_group <- factor(inci_age_sex$Age_group,
                                 levels=c("0-4","5-9","10-14","15-19","20-24",
                                          "25-29","30-34","35-39","40-44",
                                          "45-49","50-54","55-59",">=60"))

p1_C <- ggplot(data=inci_age_sex, aes(x=Age_group,y=incidence_rate,fill=Sex)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(breaks=c("Female","Male"),
                    labels=c("Female","Male"),
                    values=safe_colorblind_palette[10:11]) + 
  xlab("Age group") + ylab("Incidence rate per 100,000") + 
  labs(fill=NULL) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks=c(0,2,4,6,8),limits=c(0,8)) +
  theme(legend.position = c(0.9, 0.80))

p1 <- plot_grid(p1_A,p1_B,p1_C,labels=c("A","B","C"),nrow=3)

ggsave("figures/new/p1.png", p1,
       width = 16, height = 20, units = "cm")

### susceptibility estimates
susceptibility_estimates <- read.csv("data/susceptibility-estimates.csv",header = TRUE,
                                     stringsAsFactors = FALSE) %>% 
  mutate(age_group=factor(age_group,levels=c("0-4","5-9","10-14","15-19","20-24",
                                             "25-29","30-34","35-39","40-44",
                                             "45-49","50-54","55-59",">=60","Overall")))

p2 <- susceptibility_estimates %>% 
  ggplot(aes(x=age_group,y=susceptible_proportion*100,fill=scenario)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(breaks=c("Optimistic","Pessimistic"),
                    labels=c("Optimistic","Pessimistic"),
                    values=safe_colorblind_palette[c(7,9)]) + 
  labs(fill=NULL,y="Population susceptibility (%)",x="Age group") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks=seq(0,50,10),limits=c(0,50)) +
  theme(legend.position = c(0.9, 0.825))

ggsave("figures/new/p2.png", p2,
       width = 16, height = 7, units = "cm")