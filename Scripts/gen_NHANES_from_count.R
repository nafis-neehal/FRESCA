setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
source("./Scripts/common.R")

###### Create New Target Population by extrapolation using joint probabilities ######

summarized_nhanes <- nhanes %>% 
  group_by(Gender, Age_Group, Race_or_Ethnicity) %>% 
  summarise(n = sum(background_n))

summarized_nhanes$ratio <- summarized_nhanes$n / sum(summarized_nhanes$n)

sum_nhanes <- summarized_nhanes %>% select(Gender, Age_Group, Race_or_Ethnicity, 
                                           ratio) %>% ungroup()

target_size <- 20000

new_nhanes <- sample_n(sum_nhanes, size = target_size, replace = T, 
                       weight = ratio)

write.csv(new_nhanes, 
          "./Data/Results/NHANES_regen_20k.csv") 