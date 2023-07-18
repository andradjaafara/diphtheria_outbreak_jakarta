library(tidyverse)

# long contact matrix
cm_long <- read.csv("data/contact-matrix-idn-long.csv")
age_group <- unique(cm_long$age_contactor)
cm_long <- cm_long %>% 
  mutate(age_contactor=factor(age_contactor,levels = age_group),
         age_cotactee=factor(age_cotactee,levels = age_group))

# wide contact matrix
cm_wide <- cm_long %>% 
  pivot_wider(-c(iso3c,setting,location_contact),
              names_from=age_cotactee,values_from=mean_number_of_contacts)
# write.csv(cm_wide,"data/contact-matrix-idn-wide.csv",row.names = FALSE)