#### SOURCES ####
setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
source("./Wrapper/load_packages.R")
source("./Wrapper/wrapper_modules.R")
source("./Wrapper/load_data_formatted.R")

allhat_data <- read.csv("./Data/ALLHAT_data/ALLHAT_processed_secondary.csv")
allhat_data$X <- NULL

#params
num_exams <- 50
IPF_maxiter <- 100
weight_method <- "ipf"
scale_to <- 100

#create the treated populations
ctrl <- allhat_data %>% filter(RANDASSIGN==2) #Ch
tr1 <- allhat_data %>% filter(RANDASSIGN==3) #Am
tr2 <- allhat_data %>% filter(RANDASSIGN==4) #Li

## Split the Trial population into TA, CC and EC
## Change RANDASSIGN to 0 for controls and 1 for treated
g1 <- rbind(tr1, ctrl)
g1 <- g1 %>% mutate(RANDASSIGN = ifelse(RANDASSIGN==2, 0, 1))
g2 <- rbind(tr2, ctrl)
g2 <- g2 %>% mutate(RANDASSIGN = ifelse(RANDASSIGN==2, 0, 1))

target <- read.csv("./Data/Results/NHANES_Gen/NHANES_regen_20k.csv")

run_main_simulation <- function(num_seed, pop){
  phr_estimate_list <- c()
  for (i in 1:num_seed){
    set.seed(i)
    cat("Running Seed",i, "\n")
    new_rct <- sample_n(pop, size=nrow(pop), replace = T)
    adjusted_weights <- get_adjusted_weights(new_rct, target,
                                             proba_estimate_method="lr", 
                                             weight_method=weight_method)
    # if(sum(unlist(adjusted_weights))!=scale_to){
    #   adjusted_weights <- scale_vector(unlist(adjusted_weights), scale_to)
    # }
    phr_estimate <- get_adjusted_effect(new_rct, outcome_days_column, outcome_status_column, adjusted_weights)
    phr_estimate_list <- c(phr_estimate_list, exp(phr_estimate))
  }
  return(phr_estimate_list)
}

outcome_days_column <- "HF_EVENTDAYS"
outcome_status_column <- "HF_EVENT"

phr_estimates <- run_main_simulation(num_seed = num_exams, pop = g1)
meanCI <- CI(phr_estimates)
paste("Mean:", meanCI["mean"])
paste("Lower:", meanCI["lower"])
paste("Upper:", meanCI["upper"])

#####
#Mean: 1.3868293793936
#Lower: 1.36074981678313"
#Upper: 1.41290894200407
#


# phr_estimates <- run_main_simulation(num_seed = num_exams, pop = g2)
# meanCI <- CI(phr_estimates)
# paste("Mean:", meanCI["mean"])
# paste("Lower:", meanCI["lower"])
# paste("Upper:", meanCI["upper"])




