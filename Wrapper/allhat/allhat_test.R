#### SOURCES ####
setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
source("./Wrapper/load_packages.R")
source("./Wrapper/wrapper_modules.R")
source("./Wrapper/load_data_formatted.R")
source("./Modules/equity_metrics_miao.R")

allhat_data_NA <- read.csv("./Data/ALLHAT_data/ALLHAT_processed_with_NA.csv")
allhat_data <- read.csv("./Data/ALLHAT_data/ALLHAT_processed.csv")
allhat_data$X <- NULL
allhat_data_NA$X <- NULL

#create the treated populations
ctrl <- allhat_data %>% filter(RANDASSIGN==2) #Ch
tr1 <- allhat_data %>% filter(RANDASSIGN==3) #Am
tr2 <- allhat_data %>% filter(RANDASSIGN==4) #Li

g1 <- rbind(tr1, ctrl)
g2 <- rbind(tr2, ctrl)
g1 <- g1 %>% mutate(RANDASSIGN = ifelse(RANDASSIGN==2, 0, 1))
g2 <- g2 %>% mutate(RANDASSIGN = ifelse(RANDASSIGN==2, 0, 1))

############# Equity of whole ALLHAT data ############
# var_name_list <- c('Gender', 'Age_Group', 'Race_or_Ethnicity')
# var_subgroup_list <- list(c('Male','Female'),
#                           c('40-59','59+'),
#                           c('NH White', 'NH Black', 'NH Asian', 'Hispanic', 'Other'))
# target_data_count <- read.csv("./Data/Results/NHANES_Gen/NHANES_count_level3.csv")
# 
# Ctrl_LDM_Summary <- get_LD_nth_level(ctrl, target_data_count, NULL,
#                  var_name_list, var_subgroup_list, 1)
# Ctrl_LDM_Summary <- Ctrl_LDM_Summary %>% select(Subgroup, LD)
# Tr1_LDM_Summary <- get_LD_nth_level(tr1, target_data_count, NULL,
#                                      var_name_list, var_subgroup_list, 1)
# Tr1_LDM_Summary <- Tr1_LDM_Summary %>% select(Subgroup, LD)
# Tr2_LDM_Summary <- get_LD_nth_level(tr2, target_data_count, NULL,
#                                      var_name_list, var_subgroup_list, 1)
# Tr2_LDM_Summary <- Tr2_LDM_Summary %>% select(Subgroup, LD)

############# Biased EC Distribution Test - Protected Attributes ###################

run_simulation <- function(pop, trgrp, with_target=T){

  i <- 1
  all_data <- partition_data_generic(pop, partition_seed=i)
  target_data <- data.frame(all_data[2])
  RCT <- data.frame(all_data[[1]][[1]])
  BIASED_EC <- data.frame(all_data[[1]][[2]])
  remove(all_data)

  RCT_tr11 <- RCT %>% filter(RANDASSIGN==1) %>% mutate(grp = trgrp)
  RCT_cn1 <- RCT %>% filter(RANDASSIGN==0) %>% mutate(grp = "G2:Cn(Ch)")
  biased1 <- BIASED_EC %>% mutate(grp = "G3:Biased_EC")
  tar1 <- target_data %>% mutate(grp = "G4:Target")

  if(with_target){
    fullpop <- do.call("rbind", list(RCT_tr11 %>% select(Gender, Age_Group, Race_or_Ethnicity, grp), 
                                     RCT_cn1 %>% select(Gender, Age_Group, Race_or_Ethnicity, grp), 
                                     biased1 %>% select(Gender, Age_Group, Race_or_Ethnicity, grp), 
                                     tar1 %>% select(Gender, Age_Group, Race_or_Ethnicity, grp)))
    return(fullpop)
  } 
  else{
    rctpop <- do.call("rbind", list(RCT_tr11, RCT_cn1, biased1))
    return(rctpop)
  }

}

TA_Size <- nrow(tr1)
CC_Size <- 4000
SC_Size <- TA_Size - CC_Size
IPF_maxiter <- 100
fullpop1 <- run_simulation(pop=g1, trgrp = "G1:Trt (AM)", with_target=F)
#table1(~ Gender + Age_Group + Race_or_Ethnicity | grp, data = fullpop1, overall=F)
table1(~ Age_Group + Gender + Race_or_Ethnicity + Education + 
       Prior_HypTreat + SBP + DBP + Smoker + Had_MIS + Had_CRV + Had_ASCVD + 
       Had_MajorST + Had_T2D + Had_HDLC + Had_LVH_ELCT + Had_LVH_ECHO + 
       Has_CHD + BMI | grp, data = fullpop1, overall=F)



TA_Size <- nrow(tr2)
CC_Size <- 4000
SC_Size <- TA_Size - CC_Size
IPF_maxiter <- 100
fullpop2 <- run_simulation(pop=g2, trgrp = "G1:Trt (LI)", with_target=F)
#table1(~ Gender + Age_Group + Race_or_Ethnicity | grp, data = fullpop2, overall=F)
table1(~ Age_Group + Gender + Race_or_Ethnicity + Education + 
         Prior_HypTreat + SBP + DBP + Smoker + Had_MIS + Had_CRV + Had_ASCVD + 
         Had_MajorST + Had_T2D + Had_HDLC + Had_LVH_ELCT + Had_LVH_ECHO + 
         Has_CHD + BMI | grp, data = fullpop2, overall=F)


############# Treatment Effect Test ################
# #Calculate relative risk (Coxs Proportional Hazard Model)
# cox_model1 <- coxph(Surv(EVENTDAYS, EVENT) ~ RANDASSIGN, data = g1)
# model1_coef <- coef(cox_model1)
# ci1 <- exp(confint(cox_model1))
# hazard_ratio1 <- exp(model1_coef)
# p_value1 <- summary(cox_model1)$coefficients[5]
# 
# cox_model2 <- coxph(Surv(EVENTDAYS, EVENT) ~ RANDASSIGN, data = g2)
# model2_coef <- coef(cox_model2)
# ci2 <- exp(confint(cox_model2))
# hazard_ratio2 <- exp(model2_coef)
# p_value2 <- summary(cox_model2)$coefficients[5]

