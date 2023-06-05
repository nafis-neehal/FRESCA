setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

directory <- "./Data/Results/M6/"
patt_filename <- "PATT_Summary_medbias_2k.csv"
ldm_filename <- "LDM_Summary_medbias_NC.csv"
Final_PATT_Summary <- read.csv(paste(directory, patt_filename, sep = ""))
Final_LDM_Summary <- read.csv(paste(directory, ldm_filename, sep = ""))

# ##### Life Changing Experimental Code #####

Final_PATT_Summary_copy <- Final_PATT_Summary

Final_PATT_Summary <- Final_PATT_Summary %>% mutate(TA_CC = exp(TA_CC),
                              TA_HC = exp(TA_HC),
                              TA_HC_Prop = exp(TA_HC_Prop),
                              TA_HC_IPF = exp(TA_HC_IPF),
                              TA_HC_Both = exp(TA_HC_Both))

###########################################

confidence_interval <- function(vector, interval){
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
  return(result)
}

############### DO THE EXPONENT HERE FIRST BEFORE TAKING MEAN ##############

#PATT
confint <- 0.95
mean_patt_summary <- Final_PATT_Summary %>% 
  group_by(CC_Size) %>% 
  summarise(mean(TA_CC), 
            TA_CC_Low = confidence_interval(TA_CC, confint)[1], 
            TA_CC_High = confidence_interval(TA_CC, confint)[2],
            mean(TA_HC_Prop), 
            TA_HC_Prop_Low = confidence_interval(TA_HC_Prop, confint)[1], 
            TA_HC_Prop_High = confidence_interval(TA_HC_Prop, confint)[2],
            mean(TA_HC_IPF), 
            TA_HC_IPF_Low = confidence_interval(TA_HC_IPF, confint)[1], 
            TA_HC_IPF_High = confidence_interval(TA_HC_IPF, confint)[2], 
            mean(TA_HC_Both), 
            TA_HC_Both_Low = confidence_interval(TA_HC_Both, confint)[1], 
            TA_HC_Both_High = confidence_interval(TA_HC_Both, confint)[2]) %>%
  rename(TA_CC = "mean(TA_CC)",
         TA_HC_Prop = "mean(TA_HC_Prop)",
         TA_HC_IPF = "mean(TA_HC_IPF)",
         TA_HC_Both = "mean(TA_HC_Both)") %>%
  ungroup()

#LDM
Inf_replace <- max(Final_LDM_Summary %>% mutate_all(function(x) ifelse(is.infinite(x), -1, x)) %>% 
                     select(LDM_CC:LDM_HC_Prop_IPF)) + 0.5
LDM_Summary_Inf <- Final_LDM_Summary
LDM_Summary_Inf[sapply(LDM_Summary_Inf, is.infinite)] <- Inf_replace

mean_ldm_summary <- LDM_Summary_Inf %>%
  group_by(CC_Size, Level) %>%
  summarise(median(LDM_TA),
            LDM_TA_Low = confidence_interval(unlist(LDM_TA), 0.95)[1],
            LDM_TA_High = confidence_interval(unlist(LDM_TA), 0.95)[2],
            median(LDM_CC),
            LDM_CC_Low = confidence_interval(unlist(LDM_CC), 0.95)[1],
            LDM_CC_High = confidence_interval(unlist(LDM_CC), 0.95)[2],
            median(LDM_HC_Prop),
            LDM_HC_Prop_Low = confidence_interval(unlist(LDM_HC_Prop), 0.95)[1],
            LDM_HC_Prop_High = confidence_interval(unlist(LDM_HC_Prop), 0.95)[2],
            median(LDM_HC_IPF),
            LDM_HC_IPF_Low = confidence_interval(unlist(LDM_HC_IPF), 0.95)[1],
            LDM_HC_IPF_High = confidence_interval(unlist(LDM_HC_IPF), 0.95)[2],
            median(LDM_TA_Prop_IPF),
            LDM_TA_Prop_IPF_Low = confidence_interval(unlist(LDM_TA_Prop_IPF), 0.95)[1],
            LDM_TA_Prop_IPF_High = confidence_interval(unlist(LDM_TA_Prop_IPF), 0.95)[2],
            median(LDM_HC_Prop_IPF),
            LDM_HC_Prop_IPF_Low = confidence_interval(unlist(LDM_HC_Prop_IPF), 0.95)[1],
            LDM_HC_Prop_IPF_High = confidence_interval(unlist(LDM_HC_Prop_IPF), 0.95)[2],
  ) %>%
  rename(LDM_TA = "median(LDM_TA)",
         LDM_CC = "median(LDM_CC)",
         LDM_HC_Prop = "median(LDM_HC_Prop)",
         LDM_HC_IPF = "median(LDM_HC_IPF)",
         LDM_TA_Prop_IPF = "median(LDM_TA_Prop_IPF)",
         LDM_HC_Prop_IPF = "median(LDM_HC_Prop_IPF)") %>%
  mutate(LDM_CC_Low = ifelse(LDM_CC==Inf_replace, LDM_CC, LDM_CC_Low),
         LDM_CC_High = ifelse(LDM_CC==Inf_replace, LDM_CC, LDM_CC_High)) %>%
  mutate(Plot = ifelse(LDM_CC==Inf_replace, paste(as.character(CC_Size),"*"),
                       as.character(CC_Size))) %>%
  ungroup()

mean_patt_summary <- mean_patt_summary %>% filter(CC_Size == 500)
mean_ldm_summary <- mean_ldm_summary %>% filter(CC_Size == 500)


exp(mean_patt_summary %>% select(TA_CC:TA_HC_Both_High))

mean_ldm_summary%>%select(LDM_HC_Prop:LDM_HC_Prop_High) %>%
  summarise(LDM_HC_Prop = mean(LDM_HC_Prop),
            LDM_HC_Prop_Low = mean(LDM_HC_Prop_Low),
            LDM_HC_Prop_High = mean(LDM_HC_Prop_High))



#mean_ldm_summary <- mean_ldm_summary %>% filter(CC_Size == 500, Level == 'Female')
#exp_mean_patt_summary <- mean_patt_summary %>% select(TA_CC:TA_HC_Both_High) %>% exp()
# mean(mean_ldm_summary$LDM_CC)
# confidence_interval(mean_ldm_summary$LDM_CC, 0.95)
# mean(mean_ldm_summary$LDM_HC_Prop)
# confidence_interval(mean_ldm_summary$LDM_HC_Prop, 0.95)
# mean(mean_ldm_summary$LDM_HC_IPF)
# confidence_interval(mean_ldm_summary$LDM_HC_IPF, 0.95)
# mean(mean_ldm_summary$LDM_HC_Prop_IPF)
# confidence_interval(mean_ldm_summary$LDM_HC_Prop_IPF, 0.95)






