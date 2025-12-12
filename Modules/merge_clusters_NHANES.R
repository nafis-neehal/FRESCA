library(dplyr)

setwd("/Users/nafisneehal/ESCA/")
nhanes_data <- read.csv("./Data/NHANES/NHANES_age_adj/nhanes_hypertension.csv")

gender <- nhanes_data %>% group_by(.dots = c("Gender")) %>%
  summarise(background_n = sum(background_n)) %>%
  rename(Level = Gender) %>% 
  mutate(Var = "Gender") %>% 
  select(Var, Level, background_n)

age <- nhanes_data %>% group_by(.dots = c("Age_Group")) %>%
  summarise(background_n = sum(background_n)) %>%
  rename(Level = Age_Group) %>%
  mutate(Var = "Age_Group") %>% 
  select(Var, Level, background_n)

race <- nhanes_data %>% group_by(.dots = c("Race_or_Ethnicity")) %>%
  summarise(background_n = sum(background_n)) %>%
  rename(Level = Race_or_Ethnicity) %>%
  mutate(Var = "Race_or_Ethnicity") %>% 
  select(Var, Level, background_n)

nhanes_summary <- do.call("rbind", list(gender, age, race))

#gender
#female: 30505435, male: 24594942

#age_group 
#40-59:17182561, 59+:37917816

#race_or_ethnicity
#nh_white: 38165591, nh black: 6624866, nh asian: 2155632, hispanic: 5529827, other: 2624461

# total:55100377

# bg_rate <- Rate_Calculation(9181764, 75912636)
# user_rate <- Rate_Calculation(985, 9356)
# LDI <- Log_Disparate_Impact(as.numeric(bg_rate),as.numeric(user_rate), for_plot = FALSE)
# LDI