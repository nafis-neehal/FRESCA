setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
nhanes <- read.csv("./Data/Processed/nhanes_hypertension.csv")

summarized_nhanes <- nhanes %>% 
  group_by(Gender, Age_Group, Race_or_Ethnicity) %>% 
  summarise(n = sum(background_n))

write.csv(summarized_nhanes, 
          "./Data/Results/NHANES_count_level3.csv")
