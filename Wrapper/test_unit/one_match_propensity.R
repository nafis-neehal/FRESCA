setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
source("./Wrapper/load_packages.R")
source("./Wrapper/wrapper_modules.R")
source("./Wrapper/load_data_formatted.R")

TA_Size <- 2000
CC_Size <- 1000
SC_Size <- 1000
IPF_maxiter <- 100

all_data <- read_data(partition_seed=1)
target_data <- data.frame(all_data[2])
RCT <- data.frame(all_data[[1]][[1]])
BIASED_EC <- data.frame(all_data[[1]][[2]])

RCT <- create_dummies(RCT, column="Race_or_Ethnicity")
BIASED_EC <- create_dummies(BIASED_EC, column="Race_or_Ethnicity")

get_matched_data <- function(population, pscore){
  p <- population %>% mutate(Age_Group = factor(Age_Group),
                             Gender = factor(Gender),
                             Education = factor(Education),
                             Smoker = factor(Smoker),
                             SBP = factor(SBP),
                             FPG = factor(FPG),
                             TC = factor(TC),
                             CVDPOINTS = factor(CVDPOINTS),
                             CVDHISTORY = factor(CVDHISTORY),
                             GFRESTIMATE = factor(GFRESTIMATE))
  
  m.out <- matchit(RANDASSIGN ~ Age_Group + Gender + NH_Asian + NH_Black + NH_White + Hispanic + 
                     Other + Education +
                     Smoker + SBP + FPG + TC + CVDHISTORY + CVDPOINTS + SERUMCREAT +
                     GFRESTIMATE, data=p, distance = pscore, replace = T, caliper = 0.1)
  
  #return(m.out)
  
  matched_data <- match.data(m.out)
  return(matched_data %>% filter(RANDASSIGN==0)) #contains weights
}


#slice data to only keep propensity adjustment variables
trial_data <- RCT %>% select(MASKID, RANDASSIGN, Age_Group, Gender, NH_Asian, NH_Black, NH_White, 
                             Hispanic, Other, Education,
                             Smoker, SBP, FPG, TC, CVDHISTORY, CVDPOINTS, 
                             SERUMCREAT, GFRESTIMATE)

ec_data <- BIASED_EC %>% select(MASKID, RANDASSIGN, Age_Group, Gender, NH_Asian, NH_Black, NH_White, 
                                Hispanic, Other, Education,
                                Smoker, SBP, FPG, TC, CVDHISTORY, CVDPOINTS, 
                                SERUMCREAT, GFRESTIMATE)

#step 1: get the propensity model
propensity_model <- get_propensity_model(model_name = "logistic_regression",
                                         population=(rbind(ec_data %>% mutate(inTrial=0),
                                                           trial_data %>% mutate(inTrial=1))
                                         ))

#step 2: get propensity score
# Treated samples will have propensity score = 1
ps <- get_propensity_score(propensity_model, rbind(trial_data %>% filter(RANDASSIGN==1), ec_data)) #contains MASKID, RANDASSIGN and PRS

#step 3: get matched SC
synthetic_controls <- get_matched_data(rbind(trial_data %>% filter(RANDASSIGN==1), ec_data), unlist(ps$pr_score))


