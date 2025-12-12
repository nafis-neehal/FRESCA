#### SOURCES ####
setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
source("./Wrapper/load_packages.R")
source("./Wrapper/wrapper_modules.R")
source("./Wrapper/load_data_formatted.R")


#params
num_exams <- 50
IPF_maxiter <- 100
weight_method <- "ipf"
outcome_days_column <- "EVENTDAYS" 
outcome_status_column <- "EVENT"
scale_to <- 100

############### Step 1: Data Import #############
baseline_SPRINT_data <- read.csv("./Data/Processed/SPRINT_example.csv")
baseline_SPRINT_data <- na.omit(baseline_SPRINT_data)
primary_outcome_data <- read.csv("./Data/Processed/primarysurv.csv")
sprint_data <- baseline_SPRINT_data %>%
  inner_join(primary_outcome_data, by="MASKID") %>%
  mutate(CVDHISTORY  = factor(CVDHISTORY),
         CVDPOINTS   = if_else(CVDPOINTS<=19,'Moderate', 'High'),
         GFRESTIMATE = if_else(GFRESTIMATE>=60, 'Normal', 'Disease')) %>%
  select(MASKID:GFRESTIMATE, EVENTDAYS, EVENT)

target <- read.csv("./Data/Results/NHANES_Gen/NHANES_regen_20k.csv")

run_main_simulation <- function(num_seed, pop){
  phr_estimate_list <- c()
  for (i in 1:num_seed){
    set.seed(i)
    cat("Running Seed",i, "\n")
    new_RCT <- sample_n(pop, size=nrow(pop), replace = T)
    
    
    adjusted_weights <- get_adjusted_weights(new_RCT, target,
                                             proba_estimate_method="lr",
                                             weight_method=weight_method)
    if(sum(unlist(adjusted_weights))!=scale_to){
      adjusted_weights <- scale_vector(unlist(adjusted_weights), scale_to)
    }
    
    phr_estimate <- get_adjusted_effect(new_RCT, outcome_days_column, outcome_status_column, 
                                        adjusted_weights)
    phr_estimate_list <- c(phr_estimate_list, exp(phr_estimate))
  }
  return(phr_estimate_list)
}

phr_estimates <- run_main_simulation(num_seed = num_exams, pop = sprint_data)
meanCI <- CI(phr_estimates)
paste("Mean:", meanCI["mean"])
paste("Lower:", meanCI["lower"])
paste("Upper:", meanCI["upper"])

#ipfw = 0.797 [0.772, 0.823]
#amia = 0.798 [0.781, 0.817]


