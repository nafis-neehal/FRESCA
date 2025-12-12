#### SOURCES ####
setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
source("./Wrapper/load_packages.R")
source("./Wrapper/wrapper_modules.R")
source("./Wrapper/load_data_formatted.R")


####################################### Make Equity Summary ######################################

eq1 <- read.csv("./Data/Results/M6/ALLHAT_Results/Equity/Matchit3/EQUITY_New_full_4000_0.csv")
eq2 <- read.csv("./Data/Results/M6/ALLHAT_Results/Equity/Matchit3/EQUITY_New_full_4000_500.csv")
eq3 <- read.csv("./Data/Results/M6/ALLHAT_Results/Equity/Matchit3/EQUITY_New_full_4000_1000.csv")
eq4 <- read.csv("./Data/Results/M6/ALLHAT_Results/Equity/Matchit3/EQUITY_New_full_4000_2000.csv")
eq5 <- read.csv("./Data/Results/M6/ALLHAT_Results/Equity/Matchit3/EQUITY_New_full_4000_4000.csv")

eq1$CC_Size <- 0
eq2$CC_Size <- 500
eq3$CC_Size <- 1000
eq4$CC_Size <- 2000
eq5$CC_Size <- 4000

eq_df <- do.call("rbind", list(eq1,eq2,eq3,eq4,eq5))

########## If we want to replace INF with high value #########
target_column <- 'LD'
if (any(is.infinite(eq_df[, target_column]))) {
  eq_df[eq_df[, target_column] == Inf, target_column] <- NA
}
##############################################################

population_order <- c("Adj_TA","Only_CC", "Only_Prop", "Only_Equity", "Both")
summary_eq <- eq_df %>% group_by(CC_Size, combination, Population) %>% 
  summarise(Median   = mean(LD, na.rm=T),
            CI_Lower = CI(LD[complete.cases(LD)])[3],
            CI_Upper = CI(LD[complete.cases(LD)])[1],
            .groups='drop') %>%
  mutate(Population = factor(Population, levels = population_order)) %>%
  arrange(combination, CC_Size, Population) %>%
  select(combination, CC_Size, Population, Median, CI_Lower, CI_Upper) %>% 
  as.data.frame()

write.csv(summary_eq,
          "./Data/Results/M6/ALLHAT_Results/Equity/Matchit3/EQUITY_New_summary_4k.csv",
          row.names=F)   
  
####################################### Make ASMD Summary ######################################

smd1 <- read.csv("./Data/Results/M6/ALLHAT_Results/ASMD_Final/Matchit3/ASMD_new_weighted_4000_0.csv")
smd2 <- read.csv("./Data/Results/M6/ALLHAT_Results/ASMD_Final/Matchit3/ASMD_new_weighted_4000_500.csv")
smd3 <- read.csv("./Data/Results/M6/ALLHAT_Results/ASMD_Final/Matchit3/ASMD_new_weighted_4000_1000.csv")
smd4 <- read.csv("./Data/Results/M6/ALLHAT_Results/ASMD_Final/Matchit3/ASMD_new_weighted_4000_2000.csv")
smd5 <- read.csv("./Data/Results/M6/ALLHAT_Results/ASMD_Final/Matchit3/ASMD_new_weighted_4000_4000.csv")

smd1$CC_Size <- 0
smd2$CC_Size <- 500
smd3$CC_Size <- 1000
smd4$CC_Size <- 2000
smd5$CC_Size <- 4000

smd_df <- do.call("rbind", list(smd1,smd2,smd3,smd4,smd5))

########## If we want to replace INF with high value #########
# target_column <- 'LD'
# if (any(is.infinite(eq_df[, target_column]))) {
#   eq_df[eq_df[, target_column] == Inf, target_column] <- 2
# }
##############################################################

control_order <- c("only_cc","only_biased_ec", "only_sc", "only_prop", "only_equity", "both")
summary_smd <- smd_df %>% group_by(CC_Size, combination, Controls) %>% 
  summarise(Mean     = CI(SMD)[2],
            CI_Lower = CI(SMD)[3],
            CI_Upper = CI(SMD)[1],
            .groups='drop') %>%
  mutate(Controls = factor(Controls, levels =  control_order)) %>%
  arrange(combination, CC_Size, Controls) %>%
  select(combination, CC_Size, Controls, Mean, CI_Lower, CI_Upper) %>% 
  as.data.frame()

#summary2 <- summary_smd %>% select(CC_Size, Controls, Mean) %>% filter(Controls!="only_sc")

write.csv(summary_smd,
          "./Data/Results/M6/ALLHAT_Results/ASMD_new_weighted_summary.csv",
          row.names=F)   

  
  