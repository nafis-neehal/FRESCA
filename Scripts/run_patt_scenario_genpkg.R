setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
source("./Scripts/common.R")
source("./Modules/p_val_add.R")
source("./Modules/equity_metrics_miao.R")

#add new nhanes data
target_data <- read.csv("./Data/Results/NHANES_regen_20k.csv")

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

get_cox_effect <- function(population, weights){
  if (!is.null(weights)){
    fit <- coxph(Surv(EVENTDAYS, EVENT)~RANDASSIGN, data=population, weights = weights)
  }
  else{
    fit <- coxph(Surv(EVENTDAYS, EVENT)~RANDASSIGN, data=population)
  }
  estimate <- coef(fit)['RANDASSIGN']
  return (estimate)
}

get_gen_trial_weights <- function(trial, target){
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
  
  covariates <- c("Age_Group", "Gender", "Race_or_Ethnicity")
  
  assess_object <- assess(trial="inTrial", selection_covariates = covariates,
                          data = merged, selection_method = "lr",
                          trimpop = F)
  
  all_weights <- assess_object$weights
  weights <- all_weights[1:nrow(trial)]
  
  return (weights)
  
}

get_propensity_model <- function(population){
  m_ps <- glm(inTrial ~ Age_Group + Gender + Race_or_Ethnicity + Education +
                Smoker + SBP + FPG + TC + CVDHISTORY + CVDPOINTS + SERUMCREAT +
                GFRESTIMATE, data = population, family = binomial())
  return (m_ps)
}

get_propensity_score <- function(Model, population){
  prs_df <- data.frame(MASKID = population %>% select(MASKID), 
                       RANDASSIGN = population %>% select(RANDASSIGN), 
                       pr_score = predict(Model, population, type = "response"))
  return (prs_df)
}

get_top_n_patients <- function(Model, population, top_n){
  prs_df <- data.frame(MASKID = population %>% select(MASKID),
                       pr_score = predict(Model, population, type="response"))
  slice_df <- head(population[order(population$pr_score, decreasing=T),], n=top_n)
  return (slice_df)
}

get_matched_data <- function(population, pscore){
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
                     GFRESTIMATE, data=p, distance = pscore, replace = T)
  
  matched_data <- match.data(m.out) 
  
  return(matched_data %>% filter(RANDASSIGN==0)) #contains weights
}

get_asmd_dataframe <- function(TA, CC, EC_world, Biased_EC_world){
  
  TrialPop <- rbind(TA, CC)
  TrialPop$inTrial <- "Yes" 
  
  ExtPop <- EC_world
  ExtPop$inTrial <- "No"
  
  BExtPop <- Biased_EC_world
  BExtPop$inTrial <- "No"
  
  Pop <- rbind(TrialPop, ExtPop)
  BPop <- rbind(TrialPop, BExtPop)
  
  vars <- c("Age_Group", "Gender", "Race_or_Ethnicity", "Education",
            "Smoker", "SBP", "FPG", "TC", "CVDHISTORY", "CVDPOINTS", "SERUMCREAT",
            "GFRESTIMATE")
  
  unbiasedTab <- CreateTableOne(vars = vars, strata = "inTrial", data = Pop, test=FALSE)
  biasedTab <- CreateTableOne(vars = vars, strata = "inTrial", data = BPop, test=FALSE)
  
  dataPlot <- data.frame(variable   = rownames(ExtractSmd(unbiasedTab)),
                         Unbiased   = as.numeric(ExtractSmd(unbiasedTab)),
                         Biased     = as.numeric(ExtractSmd(biasedTab)))
  
  return (dataPlot)
  
}

run_one_patt_scenario <- function(seed_value, CC_Size){
  
  PATT_Summary <- data.frame(matrix(ncol = num_col, nrow = 0))
  colsPATT <- c("Seed", "CC_Size", "TA_CC", "TA_HC", "TA_HC_Prop", 
                "TA_HC_IPF", "TA_HC_Both")
  colnames(PATT_Summary) <- colsPATT
  
  set.seed(seed_value)
  TA <- sample_n(TA_World %>% filter(RANDASSIGN==1), TA_size)
  #W_TA_ipf <- get_IPF_weights(dat = TA, maxIter = IPF_maxiter)
  
  set.seed(seed_value)
  CC_World <- sample_n(data %>% filter(RANDASSIGN==0), TA_size)
  
  external_population <- setdiff(data, rbind(TA_World, CC_World))
  EC_world <- external_population %>% filter(RANDASSIGN==0)
  
  W_EC_bias_ipf <- get_ECWorld_Bias_weights(dat = EC_world, maxIter = IPF_maxiter)
  Biased_EC_world <- EC_world[sample(seq_len(nrow(EC_world)),
                                     size = nrow(EC_world),
                                     replace = TRUE, prob = W_EC_bias_ipf),]
  #Biased_EC_world <- EC_world
  
  ATT_TA_CC_list <- 0
  ATT_TA_HC_list <- 0
  ATT_TA_HC_prop_list <- 0
  ATT_TA_HC_IPF_list <- 0
  ATT_TA_HC_both_list <- 0
  
  flag <- 1
  n_plus_sum <- 1
  success_counter <-1
  inner_seed <- 1
  
  set.seed(seed_value)
  CC <- sample_n(CC_World, CC_Size) #fixed CC
  
  #calculate ASMD <- quantifying bias <- collect data for love plot
  #TA+CC VS EC_World <- unadjusted
  #TA + CC VS Biased_EC_World <- biased
  if(CC_Size==TA_size & seed_value==1){
    asmd_df <- get_asmd_dataframe(TA, CC, EC_world, Biased_EC_world)
    
    #save the data
    setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
    directory <- "./Data/Results/M6/"
    filename <- "ASMD_extrabias_genpkg_treat2k.csv"
    write.csv(asmd_df, paste(directory, filename, sep = ""), row.names = FALSE)
  }
  
  #train propensity here TA, CC, Biased EC World = EHR Records available
  #create the propensity score model - propensity to get in the trial
  M <- get_propensity_model(do.call("rbind", list(TA %>% mutate(inTrial = 1),
                                                  CC %>% mutate(inTrial = 1),
                                                  Biased_EC_world %>% mutate(inTrial = 0)))) #push within internal loop in specific scenario
  
  
  #start bootstrap after this point
  while (success_counter <= total_success){ ## Inner Bootstraps
    
    print(paste("n=", success_counter))
    
    #<- Collect TA_CC (no prop)
    ATT_TA_CC_list <- ATT_TA_CC_list + get_cox_effect(rbind(TA, CC), weights = NULL)
    
    for (n_plus_sum in 0:4){ #for each inner_seed value I'll try 5 times to get a pop without missing subcategory
      set.seed(inner_seed + n_plus_sum)
      #unmatched_EC <- sample_n(Biased_EC_world, randomization_ratio * TA_size - CC_Size) #pseudo-random
      recruit_needed <- randomization_ratio * TA_size - CC_Size
      unmatched_EC <- Biased_EC_world
      
      if (nrow(unmatched_EC)!=0){
        #get propensity score of this unmatched EC sample
        ps <- get_propensity_score(M, rbind(TA, unmatched_EC)) #contains MASKID, RANDASSIGN and PRS
        
        #get matched EC with weights
        matched_EC_df <- get_matched_data(rbind(TA, unmatched_EC), unlist(ps$pr_score))
        matched_EC <- matched_EC_df %>% select(MASKID:EVENT)
        matched_EC_weights <- matched_EC_df$weights #contains only EC weights
      }
      
      #create matched and unmatched HC
      unmatched_HC <- rbind(CC, unmatched_EC) 
      
      if (nrow(unmatched_EC)!=0){
        matched_HC <- rbind(CC, matched_EC) #just contains maskid:EVENT
        matched_HC_weights <- c(rep(1, nrow(CC)), matched_EC_weights)
        matched_HC_df <- matched_HC
        matched_HC_df$weights <- matched_HC_weights
      }
      else{
        matched_HC <- unmatched_HC
        matched_HC_df <- matched_HC
      }
      
      flag <- check_subgroup_counts(gender_sbgrps = 2, race_sbgrps = 5,
                                    cvd_sbgrps = 2, frs_sbgrps = 2, TA = TA, HC = matched_HC)

      # flag <- check_subgroup_counts(age_sbgrps = 2, gender_sbgrps = 2, race_sbgrps = 5, 
      #                               gfr_sbgrps = 2, TA = TA, HC = matched_HC)
      
      if (flag==1){
        break
      }
      else{
        flag<- 0
        print("Changing seed inside 1..n loop due to missing subcategory!")
      }
    }
    
    if(flag==0){
      inner_seed <- inner_seed + 5
      next
    }
    
    success_counter <- success_counter + 1
    
    #<- Collect TA_HC (no prop)
    ATT_TA_HC_list <- ATT_TA_HC_list + get_cox_effect(rbind(TA, unmatched_HC), weights = NULL)
    
    #<- Collect TA_HC_prop
    if (nrow(unmatched_EC)!=0){ #unmatched EC has rows
      ATT_TA_HC_prop_list <- ATT_TA_HC_prop_list + get_cox_effect(rbind(TA, matched_HC), 
                                                                  weights = c(rep(1, nrow(TA)), 
                                                                              matched_HC_weights))
    }
    else{ #unmatched EC has no rows
      ATT_TA_HC_prop_list <- ATT_TA_HC_prop_list + get_cox_effect(rbind(TA, matched_HC), 
                                                                  weights = NULL)
    }
    
    #make prediction for TA + Unmatched HC
    unmatched_trial <- rbind(TA, unmatched_HC)
    W_unmatched_adjust_weights <- get_gen_trial_weights(unmatched_trial, target_data) #size Treat + Control
    
    #make prediction for TA + Matched HC
    matched_trial <- rbind(TA, matched_HC)
    W_matched_adjust_weights <- get_gen_trial_weights(matched_trial, target_data) # size Treat + Control
    
    #TA_HC Equity adjustment
    ATT_TA_HC_IPF_list <- ATT_TA_HC_IPF_list + get_cox_effect(rbind(TA, unmatched_HC), 
                                                              weights = W_unmatched_adjust_weights) #Unmatched data
    
    #TA_HC Both Adjustment
    if(nrow(unmatched_EC)!=0){  # <<<<<<<<<<< ------ ERROR is here in coxph 
      propensity_weights <- c(rep(1, nrow(TA)), matched_HC_weights)
      new_weights <- (0.1 * propensity_weights + 0.9 * W_matched_adjust_weights)
      ATT_TA_HC_both_list <- ATT_TA_HC_both_list + get_cox_effect(rbind(TA, matched_HC), 
                                                                  weights = new_weights) #matched HC
    }
    else{ #when CC = full size, so no need to recruit from EC
      ATT_TA_HC_both_list <- ATT_TA_HC_both_list + get_cox_effect(rbind(TA,matched_HC), 
                                                                  weights = W_matched_adjust_weights) #unmatched HC
    }
    inner_seed <- inner_seed + 1
    
  }
  #add dataframe
  row <- as.data.frame(list(seed_value, CC_Size, 
                            ATT_TA_CC_list/total_success, 
                            ATT_TA_HC_list/total_success, 
                            ATT_TA_HC_prop_list/total_success, 
                            ATT_TA_HC_IPF_list/total_success, 
                            ATT_TA_HC_both_list/total_success))
  colnames(row) <- colsPATT
  PATT_Summary <- rbind(PATT_Summary, row) 
  return(PATT_Summary)
}


######### Mini Experiment for generating Table1 for TA_World, CC_World, EC_World
# tw_copy <- TA_World
# tw_copy$group <- "TA_World"
# ccw_copy <- CC_World
# ccw_copy$group <- "CC_World"
# ecw_copy <- Biased_EC_world
# ecw_copy$group <- "Biased_EC_World"
# 
# population <- do.call("rbind", list(tw_copy, ccw_copy, ecw_copy))
# population$Education[population$Education=='<HSG'] <- "Below HSG"
# population$CVDHISTORY <- factor(population$CVDHISTORY, levels=c(0,1),
#                                 labels=c("No", "Yes"))
# 
# table1(~ Age_Group + Gender + Race_or_Ethnicity + Education +
#          Smoker + FPG + TC + CVDHISTORY + CVDPOINTS + SERUMCREAT +
#          GFRESTIMATE | group, data=population)



