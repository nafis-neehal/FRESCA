######## Generic Helper Functions #########

scale_vector <- function(vector, n) {
  # Calculate the sum of the vector
  sum_vector <- sum(vector)
  
  # Scale each element of the vector
  scaled_vector <- (vector / sum_vector) * n
  
  return(scaled_vector)
}

###########################################

#### READ DATA FILES #####
read_data <- function(partition_seed, import_target=T){
  nhanes <- read.csv("./Data/Processed/nhanes_hypertension.csv")
  baseline_SPRINT_data <- read.csv("./Data/Processed/SPRINT_example.csv")
  baseline_SPRINT_data <- na.omit(baseline_SPRINT_data)
  primary_outcome_data <- read.csv("./Data/Processed/primarysurv.csv")
  data <- baseline_SPRINT_data %>%
    inner_join(primary_outcome_data, by="MASKID") %>%
    mutate(CVDHISTORY  = factor(CVDHISTORY),
           CVDPOINTS   = if_else(CVDPOINTS<=19,'Moderate', 'High'),
           GFRESTIMATE = if_else(GFRESTIMATE>=60, 'Normal', 'Disease')) %>%
    select(MASKID:GFRESTIMATE, EVENTDAYS, EVENTDAYS_POSTI, EVENT)
  
  dat <- get_rct_biased_ec_data(seed_value=partition_seed, data=data, TA_Size=TA_Size, CC_Size=CC_Size, 
                                IPF_maxiter=IPF_maxiter)
  
  #get the target dataframe ready
  if(import_target==F){
    return(dat)
  }
  else{
    target <- read.csv("./Data/Results/NHANES_Gen/NHANES_regen_20k.csv")
    return (list(dat, target))
  }
}

partition_data_generic <- function(data, partition_seed, import_target=T){
  dat <- get_rct_biased_ec_data(seed_value=partition_seed, data=data, TA_Size=TA_Size, CC_Size=CC_Size, 
                                IPF_maxiter=IPF_maxiter)
  #get the target dataframe ready
  if(import_target==F){
    return(dat)
  }
  else{
    target <- read.csv("./Data/Results/NHANES_Gen/NHANES_regen_20k.csv")
    return (list(dat, target))
  }
}

############# Functions to bias External Controls for ALLHAT trial #############
get_BiasAllhat <- function(dat, maxIter){
  #note: using military population ratio
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
  targets$Age_Group <- tibble(
    '1' = 40, #59+ change =   40   #NHANES=69
    '0' = 60  #40-59 change = 60   #NHANES=31
  )
  targets$Smoker <- tibble(
    '0' = 40,
    '1' = 60
  )
  targets$Had_MajorST <- tibble(
    "NO" = 40,
    "Yes" = 60
  )
  targets$Had_HDLC <- tibble(
    "NO" = 40,
    "Yes" = 60
  )
  
  var_list <- c("Race_or_Ethnicity", 'Gender', 'Age_Group', 'Smoker', 'Had_MajorST', 'Had_HDLC')
  
  dat <- dat %>% select(unlist(var_list))
  
  #recoding these two variables because subgroup names are not suitable to access using dots
  dat2 <- dat %>% mutate(Age_Group = recode(Age_Group, '40-59'=0, '59+'=1), 
                         Race_or_Ethnicity = recode(Race_or_Ethnicity,
                                                    'NH White' = 0,
                                                    'NH Black' = 1,
                                                    'NH Asian' = 2,
                                                    'Hispanic' = 3,
                                                    'Other' = 4),
                         Smoker = recode(Smoker, "No smoke"=0, "Smoke"=1))
  
  dat2 <- as_tibble(dat2)
  result <- ipu(dat2, targets, max_iterations = maxIter)
  weights <- result$weight_tbl$weight
  
  return(weights)
}

############# Functions to bias External Controls for SPRINT trial #############
get_extraBiasSprint <- function(dat, maxIter){
  #note: using military population ratio
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

get_mediumBiasSprint <- function(dat, maxIter){
  #fix the target distributions from TP
  targets <- list()
  targets$Age_Group <- tibble(
    '1' = 40, #59+ change =   40   #NHANES=69
    '0' = 60  #40-59 change = 60   #NHANES=31
  )
  targets$GFRESTIMATE <- tibble(
    'Normal' = 60.0,
    'Disease' = 40.0
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

################################################################################
###### Function to borrow synthetic controls from EC #############
get_matched_controls <- function(method_name, RCT, BIASED_EC, vars_to_match, N){
  
  # #slice data to only keep propensity adjustment variables
  # trial_data <- RCT %>% select(MASKID, RANDASSIGN, Age_Group, Gender, Race_or_Ethnicity, Education,
  #                              Smoker, SBP, FPG, TC, CVDHISTORY, CVDPOINTS, 
  #                              SERUMCREAT, GFRESTIMATE)
  # 
  # ec_data <- BIASED_EC %>% select(MASKID, RANDASSIGN, Age_Group, Gender, Race_or_Ethnicity, Education,
  #                                 Smoker, SBP, FPG, TC, CVDHISTORY, CVDPOINTS, 
  #                                 SERUMCREAT, GFRESTIMATE)
  
  trial_data <- RCT
  ec_data <- BIASED_EC
  
  if(method_name=="top_n_propensity"){
    matched_controls_df <- get_top_n_propensity_from_EC(trial_data, ec_data, N)
  }
  else if(method_name=="one_to_one_with_replace_caliper"){
    matched_controls_df <- get_one_to_one_match_from_EC(trial_data, ec_data,
                                                        replace = T, caliper = 1)
  }
  else if(method_name=="simple_matchit"){
    matched_controls_df <- get_simple_matchit4(trial_data, ec_data, vars_to_match, N)
  }
  
  return(matched_controls_df)
  #return (BIASED_EC[BIASED_EC$MASKID %in% matched_controls_df$MASKID,])
}

# get_matched_controls_with_dummies <- function(method_name, RCT, BIASED_EC, N){
#   
#   #slice data to only keep propensity adjustment variables
#   trial_data <- RCT %>% select(MASKID, RANDASSIGN, Age_Group, Gender, NH_Asian, NH_Black, NH_White, 
#                                Hispanic, Other, Education,
#                                Smoker, SBP, FPG, TC, CVDHISTORY, CVDPOINTS, 
#                                SERUMCREAT, GFRESTIMATE)
#   
#   ec_data <- BIASED_EC %>% select(MASKID, RANDASSIGN, Age_Group, Gender, NH_Asian, NH_Black, NH_White, 
#                                   Hispanic, Other, Education,
#                                   Smoker, SBP, FPG, TC, CVDHISTORY, CVDPOINTS, 
#                                   SERUMCREAT, GFRESTIMATE)
#   
#   if(method_name=="top_n_propensity"){
#     matched_controls_df <- get_top_n_propensity_from_EC(trial_data, ec_data, N)
#   }
#   else if(method_name=="one_to_one_with_replace_caliper"){
#     matched_controls_df <- get_one_to_one_match_from_EC(trial_data, ec_data,
#                                                         replace = T, caliper = 1)
#   }
#   
#   return (BIASED_EC[BIASED_EC$MASKID %in% matched_controls_df$MASKID,])
# }

############# Functions to generate matched controls #############

#this function will declare a propensity model based on population
get_propensity_model <- function(model_name, population){
  if (model_name=="logistic_regression"){
    new_population <- population[,!names(population) %in% c("MASKID", "RANDASSIGN")]
    model <- glm(inTrial ~ ., data = new_population, family = binomial())
  }
  return (model)
}

#this function will get the propensity scores of the EC population based on the 
#propensity model
get_propensity_score <- function(Model, population){
  predict_population <- population[,!names(population) %in% c("MASKID", "RANDASSIGN")]
  prs_df <- data.frame(MASKID = population %>% select(MASKID), 
                       RANDASSIGN = population %>% select(RANDASSIGN),
                       pr_score = predict(Model, predict_population, type = "response"))
  return (prs_df)
}

#this function will
## 1. declare a propensity model based on trial data
## 2. use that model to get the propensity score of all external data
## 3. return controls with top N highest propensity score
#### Allowing duplicate patients existing in Biased EC due to biasing <<<
get_top_n_propensity_from_EC <- function(trial_data, ec_data, N){
  #step 1: get the propensity model
  propensity_model <- get_propensity_model(model_name = "logistic_regression",
                                           population=(rbind(ec_data %>% mutate(inTrial=0),
                                                             trial_data %>% mutate(inTrial=1))
                                                       ))
  
  #step 2: use the model to predict ps of external data
  propensity_score_of_all_EC_df <- get_propensity_score(Model=propensity_model,
                                                        population = (ec_data %>% mutate(inTrial=0)))
  
  
  #step 3: select top N patients from EC Control group
  top_N_ec_patients <- top_n(propensity_score_of_all_EC_df, N, pr_score) %>% 
    arrange(desc(pr_score))
  
  return (top_N_ec_patients)
}

get_one_to_one_match_from_EC <- function(trial_data, ec_data, replace, caliper){
  
  #step 1: get the propensity model using the whole trial data and the biased_ec_data
  propensity_model <- get_propensity_model(model_name = "logistic_regression",
                                           population=(rbind(ec_data %>% mutate(inTrial=0),
                                                             trial_data %>% mutate(inTrial=1))
                                           ))
  
  #step 2: get propensity score
  # Treated samples will have propensity score = 1
  # contains MASKID, RANDASSIGN and PRS
  # only with Treated and EC Data
  ps <- get_propensity_score(propensity_model, rbind(trial_data %>% filter(RANDASSIGN==1), 
                                                     ec_data)) 
  
  #step 3: get matched controls
  #matched controls between Treated and EC data
  population <- rbind(trial_data %>% filter(RANDASSIGN==1), ec_data)
  
  p <- population %>% mutate(Age_Group = factor(Age_Group),
                             Gender = factor(Gender),
                             Race_or_Ethnicity = factor(Race_or_Ethnicity),
                             Education = factor(Education),
                             Smoker = factor(Smoker),
                             SBP = factor(SBP),
                             FPG = factor(FPG),
                             TC = factor(TC),
                             CVDPOINTS = factor(CVDPOINTS),
                             CVDHISTORY = factor(CVDHISTORY),
                             GFRESTIMATE = factor(GFRESTIMATE))
  
  m.out <- matchit(RANDASSIGN ~ Age_Group + Gender + Race_or_Ethnicity + Education +
                     Smoker + SBP + FPG + TC + CVDHISTORY + CVDPOINTS + SERUMCREAT +
                     GFRESTIMATE, data=p, distance = unlist(ps$pr_score), 
                   replace = replace, caliper = caliper)
  
  # m.out <- matchit(RANDASSIGN ~ Age_Group + Gender + NH_Asian + NH_Black + NH_White + 
  #                    Hispanic + Other + Education +
  #                    Smoker + SBP + FPG + TC + CVDHISTORY + CVDPOINTS + SERUMCREAT +
  #                    GFRESTIMATE, data=p, distance = unlist(ps$pr_score), 
  #                    replace = replace, caliper = caliper)
  
  matched_data <- match.data(m.out)
  return(matched_data %>% filter(RANDASSIGN==0)) #contains weights
}

get_simple_matchit <- function(trial_data, ec_data, vars_to_match){
  
  #generate the matching vars list
  
  matching_vars_list <- vars_to_match
  
  All_strings <- lapply(matching_vars_list, as.character)
  var_string <- paste(All_strings, collapse = " + ")
  response_var <- "RANDASSIGN"
  formula <- reformulate(termlabels = var_string, response = response_var)
  
  population <- rbind(trial_data %>% filter(RANDASSIGN==1),
                      ec_data)
  
  ####### Additional code to maintain Fixed SC size #########
  # tr <- trial_data %>% filter(RANDASSIGN==1)
  # cn <- trial_data %>% filter(RANDASSIGN==0)
  # tr_sample <- sample_n(tr, nrow(tr) - nrow(cn), replace = F)
  # population <- rbind(tr_sample, ec_data)
  ########################################################
  
  m1 <- matchit(formula, data = population, method = "nearest", distance = "glm", replace = T)
  matched_data <- match.data(m1)
  
  return(matched_data %>% filter(RANDASSIGN==0)) #return the matched controls data
  
}

get_simple_matchit2 <- function(trial_data, ec_data, vars_to_match){
  
  #generate the matching vars list
  
  matching_vars_list <- vars_to_match
  
  All_strings <- lapply(matching_vars_list, as.character)
  var_string <- paste(All_strings, collapse = " + ")
  trial_data$inTrial <- 1
  ec_data$inTrial <- 0
  response_var <- "inTrial"
  formula <- reformulate(termlabels = var_string, response = response_var)
  
  population <- rbind(trial_data, ec_data)
  
  m1 <- matchit(formula, data = population, method = "nearest", distance = "glm", replace = T)
  matched_data <- match.data(m1)
  
  return(matched_data %>% filter(inTrial==0)) #return the matched controls data
  
}

get_simple_matchit3 <- function(trial_data, ec_data, vars_to_match, N){ #similar to NC
  
  #generate the matching vars list
  
  if(N!=0){ #we need to borrow some synthetic controls
    matching_vars_list <- vars_to_match
    
    All_strings <- lapply(matching_vars_list, as.character)
    var_string <- paste(All_strings, collapse = " + ")
    trial_data$inTrial <- 1
    ec_data$inTrial <- 0
    response_var <- "inTrial"
    new_formula <- reformulate(termlabels = var_string, response = response_var)
    
    training_population <- do.call("rbind", list(trial_data, ec_data)) 
    
    ta_sample <- sample_n(trial_data %>% filter(RANDASSIGN==1), size=N, replace = F)
    predict_population <- rbind(ta_sample, ec_data)
    
    #step 1: get the propensity model
    # propensity_model <- get_propensity_model(model_name = "logistic_regression",
    #                                          population=training_population)
    propensity_model <- glm(new_formula, data = training_population, family = binomial())
    
    #step 2: use the model to predict ps of external data
    # propensity_score_of_predict_population <- get_propensity_score(Model=propensity_model,
    #                                                       population = predict_population)
    pr_score = predict(propensity_model, predict_population, type = "response")
    
    #step 3: get matched external con†rols
    m1 <- matchit(new_formula, data = predict_population, method = "nearest", distance = pr_score, replace = F)
    matched_data <- match.data(m1)
    return(matched_data %>% filter(inTrial==0))
  }
}
  # m1 <- matchit(formula, data = population, method = "nearest", distance = "glm", replace = T)
  # matched_data <- match.data(m1)
  # 
  # #got the propensity score from the computed model, saved as distance column in matched_data dataframe
  # ta_sample <- sample_n(matched_data$RANDASSIGN==1, size=N, replace = F)
  # ec_data <
  # population2 <- 
  # m2 <- matchit(formula, data = , method = "nearest", distance = matched)
  
  #return(matched_data %>% filter(inTrial==0)) #return the matched controls data
  

get_simple_matchit4 <- function(trial_data, ec_data, vars_to_match, N){ #similar to NC
  
  #generate the matching vars list
  
  if(N!=0){ #we need to borrow some synthetic controls
    matching_vars_list <- vars_to_match
    
    All_strings <- lapply(matching_vars_list, as.character)
    var_string <- paste(All_strings, collapse = " + ")
    trial_data$inTrial <- 1
    ec_data$inTrial <- 0
    response_var <- "inTrial"
    new_formula <- reformulate(termlabels = var_string, response = response_var)
    
    training_population <- do.call("rbind", list(trial_data, ec_data)) 
    
    ta_sample <- trial_data %>% filter(RANDASSIGN==1) #using all TA
    predict_population <- rbind(ta_sample, ec_data)
    
    #step 1: get the propensity model
    # propensity_model <- get_propensity_model(model_name = "logistic_regression",
    #                                          population=training_population)
    propensity_model <- glm(new_formula, data = training_population, family = binomial())
    
    #step 2: use the model to predict ps of external data
    # propensity_score_of_predict_population <- get_propensity_score(Model=propensity_model,
    #                                                       population = predict_population)
    pr_score = predict(propensity_model, predict_population, type = "response")
    
    #step 3: get matched external con†rols
    m1 <- matchit(new_formula, data = predict_population, method = "nearest", distance = pr_score, replace = T)
    matched_data <- match.data(m1)
    
    synthetic_controls <- matched_data %>% filter(inTrial==0)
    synthetic_controls$inTrial <- NULL
    sc_weights <- synthetic_controls$weights
    synthetic_controls <- sample_n(synthetic_controls, size = N, replace = T, weight = sc_weights)
    synthetic_controls$distance <- NULL
    synthetic_controls$weighs <- NULL
    synthetic_controls$RANDASSIGN <- 0
    
    return(synthetic_controls)
  }
}

################################################################################


############# Functions to calculate generalized weights for the RCT population #############

get_gen_trial_weights_v1 <- function(trial, target, proba_estimate_method){
  target$EVENTDAYS <- NA
  target$EVENT <- NA
  
  trial$inTrial <- 1
  target$inTrial <- 0
  
  merged <- rbind(trial %>% select(Gender, Age_Group, Race_or_Ethnicity, EVENTDAYS, EVENT, inTrial),
                  target %>% select(Gender, Age_Group, Race_or_Ethnicity, EVENTDAYS, EVENT, inTrial))
  
  merged <- merged %>% mutate(Gender = recode(Gender, "Male"=0,"Female"=1),
                              Age_Group = recode(Age_Group, "40-59"=0, "59+"=1),
                              Race_or_Ethnicity = recode(Race_or_Ethnicity,
                                                         'NH White' = 0,
                                                         'NH Black' = 1,
                                                         'NH Asian' = 2,
                                                         'Hispanic' = 3,
                                                         'Other' = 4))
  
  covariates <- c("Gender", "Age_Group", "Race_or_Ethnicity")
  
  assess_object <- assess(trial="inTrial", selection_covariates = covariates,
                          data = merged, selection_method = proba_estimate_method,
                          trimpop = F)
  
  all_weights <- assess_object$weights
  weights <- all_weights[1:nrow(trial)]
  
  return (weights)
  
}

get_gen_trial_weights_v2 <- function(trial, target, proba_estimate_method){
  target$EVENTDAYS <- NA
  target$EVENT <- NA
  
  trial$inTrial <- 1
  target$inTrial <- 0
  
  merged <- rbind(trial %>% select(Gender, Age_Group, EVENTDAYS, EVENT, inTrial),
                  target %>% select(Gender, Age_Group, EVENTDAYS, EVENT, inTrial))
  
  merged <- merged %>% mutate(Gender = recode(Gender, "Male"=0,"Female"=1),
                              Age_Group = recode(Age_Group, "40-59"=0, "59+"=1))
  
  covariates <- c("Age_Group", "Gender")
  
  assess_object <- assess(trial="inTrial", selection_covariates = covariates,
                          data = merged, selection_method = proba_estimate_method,
                          trimpop = F)
  
  all_weights <- assess_object$weights
  weights <- all_weights[1:nrow(trial)]
  
  return (weights)
  
}

get_gen_trial_weights_v3 <- function(trial, target, proba_estimate_method){
  target$EVENTDAYS <- NA
  target$EVENT <- NA
  
  trial$inTrial <- 1
  target$inTrial <- 0
  
  merged <- rbind(trial %>% select(Gender, Race_or_Ethnicity, EVENTDAYS, EVENT, inTrial),
                  target %>% select(Gender, Race_or_Ethnicity, EVENTDAYS, EVENT, inTrial))
  
  merged <- merged %>% mutate(Gender = recode(Gender, "Male"=0,"Female"=1),
                              Race_or_Ethnicity = recode(Race_or_Ethnicity,
                                                         'NH White' = 0,
                                                         'NH Black' = 1,
                                                         'NH Asian' = 2,
                                                         'Hispanic' = 3,
                                                         'Other' = 4))
  
  covariates <- c("Gender", "Race_or_Ethnicity")
  
  assess_object <- assess(trial="inTrial", selection_covariates = covariates,
                          data = merged, selection_method = proba_estimate_method,
                          trimpop = F)
  
  all_weights <- assess_object$weights
  weights <- all_weights[1:nrow(trial)]
  
  return (weights)
  
}

get_gen_trial_weights_v4 <- function(trial, target, proba_estimate_method){
  target$EVENTDAYS <- NA
  target$EVENT <- NA
  
  trial$inTrial <- 1
  target$inTrial <- 0
  
  selected_cols <- c("Gender", "Age_Group", "NH_White", "NH_Black", "NH_Asian", "Hispanic", "Other",
                     "EVENTDAYS", "EVENT", "inTrial")
  
  merged <- rbind(trial %>% select(selected_cols),
                  target %>% select(selected_cols))
  
  merged <- merged %>% mutate(Gender = recode(Gender, "Male"=0,"Female"=1),
                              Age_Group = recode(Age_Group, "40-59"=0, "59+"=1))
  
  covariates <- c("Gender", "Age_Group", "NH_White", "NH_Black", "NH_Asian", "Hispanic", "Other")
  
  assess_object <- assess(trial="inTrial", selection_covariates = covariates,
                          data = merged, selection_method = proba_estimate_method,
                          trimpop = F)
  
  all_weights <- assess_object$weights
  weights <- all_weights[1:nrow(trial)]
  
  return (weights)
  
}

get_IPF_weights <- function(dat, maxIter=100){
  #fix the target distributions from TP
  targets <- list()
  targets$Age_Group <- tibble(
    '0' = 31.184,
    '1' = 68.816
  )
  targets$Gender <- tibble(
    'Female' = 55.363,
    'Male' = 44.637
  )
  targets$Race_or_Ethnicity <- tibble(
    '0' = 69.266,
    '1' = 12.023,
    '2' = 3.912,
    '3' = 10.035,
    '4' = 4.763
  )
  
  var_list <- c("Age_Group", "Gender", "Race_or_Ethnicity")
  
  dat <- dat %>% select(unlist(var_list))
  
  dat2 <- dat %>% mutate(Age_Group = recode(Age_Group,
                                            '40-59' = 0,
                                            '59+' = 1),
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

################################################################################
# get_cox_effect <- function(population, weights){
#   if (!is.null(weights)){
#     fit <- coxph(Surv(CVD_EVENTDAYS, CVD_EVENT)~RANDASSIGN, data=population, weights = weights)
#   }
#   else{
#     fit <- coxph(Surv(CVD_EVENTDAYS, CVD_EVENT)~RANDASSIGN, data=population)
#   }
#   estimate <- coef(fit)['RANDASSIGN']
#   return (estimate)
# }

get_cox_effect <- function(population, outcome_days, outcome_event, weights){
  
  # Extract survival time and status vectors from the dataframe
  survival_time <- population[[outcome_days]]
  survival_status <- population[[outcome_event]]
  
  # Create a survival object
  survival_data <- Surv(time = survival_time, event = survival_status)
  
  if (!is.null(weights)){
    fit <- coxph(survival_data~RANDASSIGN, data=population, weights = weights)
  }
  else{
    fit <- coxph(survival_data~RANDASSIGN, data=population)
  }
  
  estimate <- coef(fit)['RANDASSIGN']
  return (estimate)
}

################################################################################

#### ADJUST HYBRID RCT DATA TO MATCH TARGET #####

get_adjusted_weights <- function(hybrid_RCT, target, proba_estimate_method, weight_method="gen"){
  #get the adjustment weights for the hybrid_RCT
  #methods exist = lr, rf, lass0
  #https://rdrr.io/github/katiecoburn/generalizeR/man/assess.html
  if (weight_method=="gen"){
    RCT_weights <- get_gen_trial_weights_v1(hybrid_RCT, target, 
                                       proba_estimate_method = proba_estimate_method)
  } else if (weight_method=="ipf"){
    RCT_weights <- get_IPF_weights(hybrid_RCT, maxIter = IPF_maxiter)
  }
  return (RCT_weights)
}

get_adjusted_effect <- function(hybrid_RCT, time, status, RCT_weights){
  
  #treatment effect using the hybrid RCT dataset
  treatment_effect <- get_cox_effect(hybrid_RCT, time, status, RCT_weights)
  return(treatment_effect)
  
}

################################################################################

#### Adjusted Equity Generation Helper functions #####

generate_subgroups <- function(V, n) {
  if (n == 1) {
    # Flatten all vectors into a single list
    result <- unlist(V)
  } else {
    # Generate combinations of elements from different vectors up to level n
    result <- combn(V, n, function(x) {
      apply(expand.grid(x), 1, paste, collapse = "/")
    }, simplify = FALSE)
  }
  
  # Remove duplicates
  result <- unique(unlist(result))
  
  return(result)
}

subgroup_in_which_var <- function(var_name_list, var_subgroup_list, name){
  for (i in 1:length(var_subgroup_list)){
    list_members <- var_subgroup_list[i]
    if (unlist(name) %in% unlist(list_members)){
      return(var_name_list[i])
    }
  }
}

get_LD_nth_level <- function(trial_pop, target_pop, trial_weights, var_name_list, var_subgroup_list, n){
  
  #create LDM summary matrix
  df <- data.frame(matrix(ncol = 4, nrow = 0))
  x <- c("Subgroup", "Background_Rate", "Observed_Rate", "LD")
  colnames(df) <- x
  
  #combine trial data and weights into one df for slicing easiness
  trial_pop$weights <- trial_weights
  
  #generate var combinations on n-th level subgroup 
  output <- generate_subgroups(var_subgroup_list, n)
  
  for (i in 1:length(output)){ #do action for each n-th level subgroup
    
    #trial_copy_grp_member <- trial_pop
    trial_copy_grp_member <- sample_n(trial_pop, size = nrow(trial_pop), replace = T, weight = trial_weights)
    target_copy_grp_member <- target_pop
    
    #take each of the subgroup combination
    #Female/NH Black -> ("Female", "NH Black")
    sub_group_name <- as.list(strsplit(output[i], "/")[[1]])
    for (j in 1:length(sub_group_name)){
      var_name <- subgroup_in_which_var(var_name_list, var_subgroup_list, sub_group_name[j])
      trial_copy_grp_member <- trial_copy_grp_member %>% filter(!!sym(var_name)==sub_group_name[j])
      target_copy_grp_member <- target_copy_grp_member %>% filter(!!sym(var_name)==sub_group_name[j])
    }
    
    if((n==1)|(n==length(var_name_list))){ #if we use first level or lowest-level subgroup
      if (is.null(trial_weights)){
        trial_grp_member_count <- nrow(trial_copy_grp_member)
        trial_total_count <- nrow(trial_pop)
      }
      else{
        # trial_grp_member_count <- ceiling(sum(trial_copy_grp_member$weights))
        # trial_total_count <- ceiling(sum(trial_weights))
        trial_grp_member_count <- nrow(trial_copy_grp_member)
        trial_total_count <- nrow(trial_pop)
      }

      target_grp_member_count <- sum(target_copy_grp_member$n)
      target_total_count <- sum(target_pop$n)
      
      #cat(target_grp_member_count, target_total_count, trial_grp_member_count, trial_total_count, "\n", sep=" ")
      
      bg_rate <- Rate_Calculation(target_grp_member_count, target_total_count)
      user_rate <- Rate_Calculation(trial_grp_member_count, trial_total_count)
      LDM <- Log_Disparate_Impact(as.numeric(round(bg_rate, 4)), as.numeric(round(user_rate,4)), for_plot = FALSE)
      df[nrow(df) + 1,] <- c(output[i], as.numeric(round(bg_rate, 4)), as.numeric(round(user_rate,4)), 
                              as.numeric(round(abs(LDM),4)))
    }
  }
  return(df)
}

################################################################################

create_dummies <- function(df, column){
  df_with_dummy <- dummy_cols(df, select_columns = column, 
                              remove_selected_columns = T)
  for ( col in 1:ncol(df_with_dummy)){
    #removes main var name from dummy columns
    colnames(df_with_dummy)[col] <-  sub(paste(column,"_", sep = ""), "", 
                                         colnames(df_with_dummy)[col]) 
    #replaces space with _ in each dummy column
    colnames(df_with_dummy)[col] <-  sub(" ", "_",                     
                                         colnames(df_with_dummy)[col])
  }
  
  return(df_with_dummy)
}

#used only for equity measurement
recode_race <- function(df){
  df$NH_Asian <- ifelse(df$NH_Asian==1,"NH_Asian_1", "NH_Asian_0")
  df$NH_White <- ifelse(df$NH_White==1,"NH_White_1", "NH_White_0")
  df$NH_Black <- ifelse(df$NH_Black==1,"NH_Black_1", "NH_Black_0")
  df$Hispanic <- ifelse(df$Hispanic==1,"Hispanic_1", "Hispanic_0")
  df$Other    <- ifelse(df$Other==1,"Other_1", "Other_0")
  return (df)
}



################################################################################
## Can Handle CC Size = 0
get_ASMD <- function(pop1, pop2, seed, controls, var_list){
  pop1$group <- "c1"
  if (nrow(pop2)!=0){
    pop2$group <- "c2"
    v <- CreateTableOne(var_list, strata = "group", data = rbind(pop1, pop2), test=F)
    df <- data.frame(Seed        = seed,
                     Controls    = controls,
                     Variable    = rownames(ExtractSmd(v)),
                     SMD         = as.numeric(ExtractSmd(v)))
  }
  else{
    df <- data.frame(Seed        = seed,
                     Controls    = controls,
                     Variable    = var_list,
                     SMD         = NA)
  }
  return(df)
}


# get_weighted_mean <- function(values, weights) {
#   print(values)
#   print(weights)
#   print(sum(values * weights))
#   wm <- sum(values * weights) / sum(weights)
#   return(wm)
# }


##### Don;t need this <- because we won't be calculated weighted var for the controls, 
##### as we are just using the treatment arm
# get_weighted_var <- function(values, weights, mean) {
#   return(sum(weights * (values - mean)^2) / sum(weights))
# }




get_ASMD_new <- function(pop1, pop2, seed, controls, var_list, control_weights){
  
  
  #smd_data <- data.frame(Seed = numeric(0), Controls = character(0), Variable = character(0), SMD = numeric(0), stringsAsFactors = FALSE)
  
  if(nrow(pop2)==0){
    
    df <- data.frame(Seed = seed, 
                     Controls = controls, 
                     Variable = var_list, 
                     SMD = NA)
    return(df)
  }
  else{
    
    pop1$tr <- 1
    pop2$tr <- 0
    combined_pop <- rbind(pop1, pop2)
    cmd <- combined_pop %>% select(var_list) #selected vars
    
    
    res <- bal.tab(cmd, treat = combined_pop$tr, s.d.denom = "treat", binary="std", 
                   abs=T, weights = c(rep(1, nrow(pop1)),control_weights))
    
    smd_data <- data.frame(Seed = seed,
                           Controls = controls,
                           Variable = rownames(res$Balance),
                           SMD = res$Balance$Diff.Adj)
    return(smd_data)
  }
}


check_subgroup_counts <- function(data, columns, unique_counts) {
  for (i in seq_along(columns)) {
    col <- columns[i]
    unique_vals <- unique(data[[col]])
    
    if (length(unique_vals) != unique_counts[i]) {
      return(FALSE)
    }
  }
  
  return(TRUE)
}

