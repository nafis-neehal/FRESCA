setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
source("./Scripts/common.R")
source("./Scripts/run_patt_scenario.R")
source("./Modules/p_val_add.R")
source("./Modules/equity_metrics_miao.R")

get_ECWorld_Bias_weights_extra <- function(dat, maxIter){
  #fix the target distributions from TP
  targets <- list()
  targets$Race_or_Ethnicity <- tibble(
    '1' = 20.2, #black original = 28, change = 20.2
    '0' = 54, #white original = 57, change = 54
    '3' = 17.2,  #hispanic original = 12, change = 17.2
    '4' = 1.7, #other orig = 0.3, change = remaining
    '2' = 6.9   #asian orig = 1, change = 6.9
  )
  targets$Gender <- tibble(
    'Male' = 72.6, #change=72.6
    'Female' = 27.4   #change=27.4
  )
  targets$CVDPOINTS <- tibble(
    'High' = 45.0,
    'Moderate' = 55.0
  )
  targets$CVDHISTORY <- tibble(
    '1' = 45.000,
    '0' = 55.000
  )
  
  var_list <- c("Race_or_Ethnicity", 'Gender', 'CVDPOINTS', 'CVDHISTORY')
  
  dat <- dat %>% select(unlist(var_list))
  
  #Age_Group = recode(Age_Group, '40-59'=0, '59+'=1),
  dat2 <- dat %>% mutate(Race_or_Ethnicity = recode(Race_or_Ethnicity,
                                                    'NH White' = 0,
                                                    'NH Black' = 1,
                                                    'NH Asian' = 2,
                                                    'Hispanic' = 3,
                                                    'Other' = 4))
  
  dat2 <- as_tibble(dat2)
  result <- ipu(dat2, targets, max_iterations = maxIter)
  weights <- result$weight_tbl$weight
  
  return(weights)
}

get_ECWorld_Bias_weights_med <- function(dat, maxIter){
  #fix the target distributions from TP
  targets <- list()
  targets$Age_Group <- tibble(
    '1' = 40, #59+ =   40   #NHANES=69
    '0' = 60  #40-59 = 60   #NHANES=31
  )
  targets$Race_or_Ethnicity <- tibble(
    '1' = 20.2, #black original = 28, change = 20.2
    '0' = 54, #white original = 57, change = 54
    '3' = 17.2,  #hispanic original = 12, change = 17.2
    '4' = 1.7, #other orig = 0.3, change = remaining
    '2' = 6.9   #asian orig = 1, change = 6.9
  )
  targets$Gender <- tibble(
    'Male' = 72.6, #change=72.6
    'Female' = 27.4   #change=27.4
  )
  targets$GFRESTIMATE <- tibble(
    'Normal' = 60.0,
    'Disease' = 40.0
  )
  
  var_list <- c("Race_or_Ethnicity", 'Gender', 'Age_Group', 'GFRESTIMATE')
  
  dat <- dat %>% select(unlist(var_list))
  
  dat2 <- dat %>% mutate(Age_Group = recode(Age_Group, '40-59'=0, '59+'=1),
                         Race_or_Ethnicity = recode(Race_or_Ethnicity,
                                                    'NH White' = 0,
                                                    'NH Black' = 1,
                                                    'NH Asian' = 2,
                                                    'Hispanic' = 3,
                                                    'Other' = 4))
  
  dat2 <- as_tibble(dat2)
  result <- ipu(dat2, targets, max_iterations = maxIter)
  weights <- result$weight_tbl$weight
  
  return(weights)
}

TA_World <- data %>% filter(RANDASSIGN==1)
TA_size <- 1000
total_seeds <- 15 #number of total scenarios
total_success <- 30 #number of inner bootstraps
IPF_maxiter <- 100
step_size <- 250
randomization_ratio <- 1
num_col <- 7
CC_size_list <- seq(TA_size, 0, -step_size)
EC_size_list <- randomization_ratio*TA_size - CC_size_list

seed_value <- 1
set.seed(seed_value)
TA <- sample_n(TA_World %>% filter(RANDASSIGN==1), TA_size)

set.seed(seed_value)
CC_World <- sample_n(data %>% filter(RANDASSIGN==0), TA_size)

external_population <- setdiff(data, rbind(TA_World, CC_World))
EC_world <- external_population %>% filter(RANDASSIGN==0) # no bias

W_EC_bias_ipf_extra <- get_ECWorld_Bias_weights_extra(dat = EC_world, maxIter = IPF_maxiter)
W_EC_bias_ipf_med <- get_ECWorld_Bias_weights_med(dat = EC_world, maxIter = IPF_maxiter)

Biased_EC_world_extra <- EC_world[sample(seq_len(nrow(EC_world)),
                                   size = nrow(EC_world),
                                   replace = TRUE, prob = W_EC_bias_ipf_extra),]
Biased_EC_world_med <- EC_world[sample(seq_len(nrow(EC_world)),
                                         size = nrow(EC_world),
                                         replace = TRUE, prob = W_EC_bias_ipf_med),]


TA_World$group <- "TA"
CC_World$group <- "CC"
EC_world$group <- "EC (no bias)"
Biased_EC_world_med$group <- "Biased EC (Sample 1)"
Biased_EC_world_extra$group <- "Biased EC (Sample 2)"

population <- do.call("rbind", list(TA_World, CC_World, EC_world, Biased_EC_world_med, 
                                    Biased_EC_world_extra))

population$Education[population$Education=='<HSG'] <- "Below HSG"
population$CVDHISTORY <- factor(population$CVDHISTORY, levels=c(0,1),
                                labels=c("No", "Yes"))

population$group <- factor(population$group, levels = c("TA", "CC", "EC (no bias)", 
                                                        "Biased EC (Sample 1)", "Biased EC (Sample 2)"))

label(population$Age_Group) <- "Age Group"
label(population$Gender) <- "Gender"
label(population$Race_or_Ethnicity) <- "Race or Ethnicity"
label(population$Education) <- "Educational Status"
label(population$Smoker) <- "Smoker?"
label(population$FPG) <- "Fasting Glucose level"
label(population$TC) <- "Total Cholesterol"
label(population$SBP) <- "Average of 3 sitting Systolic BP"
label(population$CVDHISTORY) <- "Has Clinical or Subclinical CVD"
label(population$CVDPOINTS) <- "Framingham Risk Score"
label(population$SERUMCREAT) <- "Serum Creatinine mg/dL"
label(population$GFRESTIMATE) <- "Estimated GFR within past 6 months"

table1(~ Age_Group + Gender + Race_or_Ethnicity + Education | group, 
       data=population, overall=F)

table1(~ Smoker + FPG + TC + SBP + CVDHISTORY + CVDPOINTS + SERUMCREAT +
                  GFRESTIMATE | group, data=population, overall=F)









