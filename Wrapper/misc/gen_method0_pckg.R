#### SOURCES ####
setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
source("./Wrapper/load_packages.R")
source("./Wrapper/wrapper_modules.R")
source("./Wrapper/load_data_formatted.R")


#params
num_exams <- 50
IPF_maxiter <- 100
weight_method <- "gen"

#Whole Sprint data
baseline_SPRINT_data <- read.csv("./Data/Processed/SPRINT_example.csv")
baseline_SPRINT_data <- na.omit(baseline_SPRINT_data)
primary_outcome_data <- read.csv("./Data/Processed/primarysurv.csv")
data <- baseline_SPRINT_data %>%
  inner_join(primary_outcome_data, by="MASKID") %>%
  mutate(CVDHISTORY  = factor(CVDHISTORY),
         CVDPOINTS   = if_else(CVDPOINTS<=19,'Moderate', 'High'),
         GFRESTIMATE = if_else(GFRESTIMATE>=60, 'Normal', 'Disease')) %>%
  select(MASKID:GFRESTIMATE, EVENTDAYS, EVENTDAYS_POSTI, EVENT)

target <- read.csv("./Data/Results/NHANES_Gen/NHANES_regen_20k.csv")

#data <- create_dummies(data, column = "Race_or_Ethnicity")
#target <- create_dummies(target, column = "Race_or_Ethnicity")

phr_estimate_list <- c()
for (i in 1:num_exams){
  set.seed(i)
  cat("Running Seed",i, "\n")
  new_rct <- sample_n(data, size=nrow(data), replace = T)
  adjusted_weights <- get_adjusted_weights(new_rct, target,
                                           proba_estimate_method="lr", 
                                           weight_method=weight_method)
  phr_estimate <- get_adjusted_effect(new_rct, adjusted_weights)
  phr_estimate_list <- c(phr_estimate_list, exp(phr_estimate))
}

meanCI <- CI(phr_estimate_list)
paste("Mean:", meanCI["mean"])
paste("Lower:", meanCI["lower"])
paste("Upper:", meanCI["upper"])

####### Result (100 runs) - IPTW Re-Weighting ######
#"Mean: 0.762557780181857"
#"Lower: 0.748128695242101"
#"Upper: 0.776986865121613"

####### Result (100 runs) - IPF Re-Weighting######
#"Mean: 0.799793761624687"
#"Lower: 0.781592237249894"
#"Upper: 0.81799528599948"
