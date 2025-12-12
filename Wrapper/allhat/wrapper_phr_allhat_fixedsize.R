#### SOURCES ####
setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
source("./Wrapper/load_packages.R")
source("./Wrapper/wrapper_modules.R")
source("./Wrapper/load_data_formatted.R")

allhat_data <- read.csv("./Data/ALLHAT_data/ALLHAT_processed_secondary.csv")
allhat_data$X <- NULL

### Factorize variables
vars_to_convert <- c("Age_Group", "Gender", "Race_or_Ethnicity", "Education", 
                     "Prior_HypTreat","Smoker","Had_MIS", "Had_CRV","Had_ASCVD",
                     "Had_MajorST","Had_T2D","Had_HDLC", "Had_LVH_ELCT","Had_LVH_ECHO","Has_CHD")

### Matching variables
vars_to_match <- c("Age_Group", "Gender", "Race_or_Ethnicity", "Education", 
                   "Prior_HypTreat", "SBP", "DBP", "Smoker","Had_MIS", "Had_CRV","Had_ASCVD",
                   "Had_MajorST","Had_T2D","Had_HDLC", "Had_LVH_ELCT","Had_LVH_ECHO","Has_CHD", "BMI")

#convert variables to factors
for (var in vars_to_convert) {
  allhat_data[[var]] <- as.factor(allhat_data[[var]])
}

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

#### Run Main Simulation #####
run_main_simulation <- function(num_seed, pop){
  
  #phr_estimate_list <- c()
  only_cc <- c()
  only_syn <- c()
  only_prop <- c()
  only_equity <- c()
  both <- c()
  seeds_used <- c()
  
  for (i in 1:num_seed){
    
    zero_cc <- F
    zero_sc <- F 
    
    cat("\nIteration:", i, "\n")
    
    #Step 1: create data
    all_data <- partition_data_generic(pop, partition_seed=i)
    target_data <- data.frame(all_data[2])
    RCT <- data.frame(all_data[[1]][[1]])
    BIASED_EC <- data.frame(all_data[[1]][[2]])
    
    remove(all_data)
    cat("Data partitioned...\n")
    
    if(CC_Size == 0){
      zero_cc <- T 
    }
    if(SC_Size == 0){
      zero_sc <- T
    }
    
    ###### only cc effect -> append
    only_cc <- c(only_cc, exp(get_adjusted_effect(RCT, outcome_days_column, outcome_status_column, NULL)))
    ###############################
    
    #Step 2: Borrow Synthetic Controls
    if(zero_sc==F){ #we are borrowing some SC
      synthetic_controls <- get_matched_controls("simple_matchit",
                                                 RCT, BIASED_EC,
                                                 vars_to_match, SC_Size)
      # sc_weights <- synthetic_controls$weights
      # synthetic_controls <- sample_n(synthetic_controls, size = SC_Size, replace = T, weight = sc_weights)
      # 
      # synthetic_controls$distance <- NULL
      # synthetic_controls$weighs <- NULL
      # synthetic_controls$RANDASSIGN <- 0
      synthetic_controls <- synthetic_controls %>% select(colnames(RCT))
      
      only_syn <- c(only_syn, exp(get_adjusted_effect(rbind(RCT %>% filter(RANDASSIGN==1), synthetic_controls), 
                                                      outcome_days_column, outcome_status_column, NULL)))
      
      cat(nrow(synthetic_controls), "Synthetic Controls borrowed...\n", sep = " ")
      
      hybrid_RCT <- rbind(RCT, synthetic_controls)
    }
    
    else{ #when zero_sc == True, that means no synthetic controls borrowed, all CC
      only_syn <- c(only_syn, NA)
      
      hybrid_RCT <- RCT
    }
    
    ###### only propensity effect -> append
    only_prop <- c(only_prop, exp(get_adjusted_effect(hybrid_RCT, outcome_days_column, outcome_status_column, NULL)))
    
    #### Check if subgroups are all there before IPF adjustment
    columns_to_check <- c("Age_Group", "Gender", "Race_or_Ethnicity")
    unique_counts_to_check <- c(2, 2, 5)
    flag1 <- check_subgroup_counts(RCT %>% filter(RANDASSIGN == 0), columns_to_check, unique_counts_to_check) # only CC
    flag2 <- check_subgroup_counts(hybrid_RCT %>% filter(RANDASSIGN == 0), columns_to_check, unique_counts_to_check) # SC + SC
    
    if ((zero_cc==F) && (flag1==F || flag2==F)){
      only_cc <- head(only_cc, -1)
      only_syn <- head(only_syn, -1)
      only_prop <- head(only_prop, -1)
      cat("Missing subgroup found in Age/Race/Gender with Non-Zero CC \n", sep=" ")
      next
    }
    else if((zero_cc==T) && (flag2==F)){
      only_cc <- head(only_cc, -1)
      only_syn <- head(only_syn, -1)
      only_prop <- head(only_prop, -1)
      cat("Missing subgroup found in Age/Race/Gender with Zero CC \n", sep=" ")
      next
    }
    
    seeds_used <- c(seeds_used, i)
    
    
    #Step 3: Measure adjusted Treatment Effect
    #Use lr/rf. Lasso takes a lot of time - when using IPTW
    
    #------> 3.1 Get adjusted weights for <Treatment arm>
    
    tr_adjusted_weights <- get_adjusted_weights(RCT %>% filter(RANDASSIGN == 1), target_data,
                                                proba_estimate_method="lr", weight_method=weight_method)
    
    # if(sum(unlist(tr_adjusted_weights))!=scale_to){
    #   tr_adjusted_weights <- scale_vector(unlist(tr_adjusted_weights), scale_to)
    # }
    
    #------> 3.2 Get adjusted weights for <Just CC>
    if(zero_cc==F){
      cc_adjusted_weights <- get_adjusted_weights(RCT %>% filter(RANDASSIGN == 0), target_data,
                                                  proba_estimate_method="lr", weight_method=weight_method)
      # if(sum(unlist(cc_adjusted_weights))!=scale_to){
      #   cc_adjusted_weights <- scale_vector(unlist(cc_adjusted_weights), scale_to)
      # }
    }
    
    #------> 3.3 Get adjusted weights for <Hybrid CC, SC>
    hybrid_adjusted_weights <- get_adjusted_weights(hybrid_RCT %>% filter(RANDASSIGN == 0), target_data,
                                                    proba_estimate_method="lr", weight_method=weight_method)
    # if(sum(unlist(hybrid_adjusted_weights))!=scale_to){
    #   hybrid_adjusted_weights <- scale_vector(unlist(hybrid_adjusted_weights), scale_to)
    # }
    
    cat("Adjusted weights generated...\n")
    
    #Step 4: Measure Adjusted Effects
    
    ###### only equity effect -> append
    #------> 4.1 Calculate only equity adjusted effect TA + CC + Tr weights + CC Weights
    if (zero_cc==F){
      only_equity <- c(only_equity, exp(get_adjusted_effect(RCT, outcome_days_column, outcome_status_column, c(tr_adjusted_weights, cc_adjusted_weights))))
    }
    else{
      only_equity <- c(only_equity, NA) #when there is no CC, CC_Size = 0, only equity would be NA
    }
    
    #------> 4.2 Calculate both adjusted effect TA + Hybrid CC, SC + Tr weights + Hybrid CC, SC weights
    both <- c(both, exp(get_adjusted_effect(hybrid_RCT, outcome_days_column, outcome_status_column, c(tr_adjusted_weights, hybrid_adjusted_weights))))
    
    cat("Adjusted effect estimated...\n", sep = " ")
    
  }
  
  final_df <- data.frame(seed = seeds_used,
                         CC_only = only_cc,
                         SYN_only = only_syn,
                         PROP_only = only_prop,
                         EQUITY_only = only_equity,
                         BOTH = both)
  
  return (final_df)
}

######### mini-helpers #######

calc_mean_ci <- function(x) {
  if(all(is.na(x))) {
    mean_val <- NA
    lower_ci <- NA
    upper_ci <- NA
  }
  else{
    mean_val <- mean(x, na.rm = TRUE)
    ci <- t.test(x, conf.level = 0.95)$conf.int
    lower_ci <- ci[1]
    upper_ci <- ci[2]
  }
  return(c(mean = mean_val, lower_ci = lower_ci, upper_ci = upper_ci))
}

get_summary <- function(df){
  summary_df <- df %>%
    select(CC_only:BOTH) %>%
    summarize_all(.funs = calc_mean_ci)
  return(summary_df)
}

get_results_df <- function(TA_Size, CC_Size, population, combi_str){
  
  ############ 50% reduction ############
  start <- Sys.time()
  
  result <- run_main_simulation(num_seed=num_seed, pop=population)
  result$combination <- combi_str
  end <- Sys.time()
  
  paste("Total time elapsed:", end-start)
  
  return(result)
}

get_summary_df <- function(result, combi_str){
  summary_df <- round(get_summary(result), digits = 2)
  summary_df$combination <- combi_str
  summary_df$values <- c("mean", "lower_CI", "upper_CI")
  
  return(summary_df)
}

############################ Global: Big Umbrella in Loop Varying CC Size #############
TA_Size_list <- c(3500, 3000, 2500, 2000)
CC_Size_list <- c(500, 1000, 1500, 2000)
SC_Size_list <- c(3000, 2000, 1000, 0)

for (i in seq(1, length(TA_Size_list), 1)){
  
  #### STUDY PARAMETERS ####
  IPF_maxiter <- 100
  num_seed <- 50
  scale_to <- 100
  TA_Size <- TA_Size_list[i]
  CC_Size <- CC_Size_list[i]
  SC_Size <- SC_Size_list[i]
  outcome_days_column <- "HF_EVENTDAYS"
  outcome_status_column <- "HF_EVENT"
  outcome_name <- "HF"
  weight_method <- "ipf" #"gen" or "ipf"
  
  ################################## Trial 1: Li vs CH ####################################
  
  combination <- "TR1"
  cat("\nTrial:", combination,"vs Controls \n")
  cat("Starting simulation for TA =", TA_Size, "and CC =", CC_Size,"\n", sep = " ")
  population <- g1
  result1 <- get_results_df(TA_Size, CC_Size, population, combination) 
  summary_df1 <- get_summary_df(result1, combination)
  
  ################################## Trial 2: Li vs CH ####################################
  
  # combination <- "TR2"
  # cat("\nTrial:", combination,"vs Controls \n")
  # cat("Starting simulation for TA =", TA_Size, "and CC =", CC_Size,"\n", sep = " ")
  # population <- g2
  # result2 <- get_results_df(TA_Size, CC_Size, population, combination) 
  # summary_df2 <- get_summary_df(result2, combination)
  
  
  ####################### SAVE Results ####################
  
  file_string_variable <- paste(TA_Size,"_",CC_Size,".csv", sep = "")
  cat("Saving File", file_string_variable, "...\n")
  
  write.csv(do.call("rbind", list(result1)),
            paste('./Data/Results/M6/ALLHAT_Results/Samecost/trial_size_4k/Matchit3/',outcome_name,'_samecost_full_', file_string_variable, sep = ""),
            row.names=F) 
  write.csv(do.call("rbind", 
                    list(summary_df1)) %>% select(combination, values, CC_only:BOTH),
            paste('./Data/Results/M6/ALLHAT_Results/Samecost/trial_size_4k/Matchit3/',outcome_name,'_samecost_summary_', file_string_variable, sep= ""),
            row.names=F) 
}


