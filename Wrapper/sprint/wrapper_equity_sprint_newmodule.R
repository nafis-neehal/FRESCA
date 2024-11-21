#### SOURCES ####
setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
source("./Wrapper/load_packages.R")
source("./Wrapper/wrapper_modules.R")
source("./Wrapper/load_data_formatted.R")
source("./Wrapper/load_data_formatted.R")
source("./Modules/equity_metrics_miao.R")

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

### Factorize variables
vars_to_convert <- c("Age_Group", "Gender", "Race_or_Ethnicity", "Education", 
                     "Smoker", "SBP", "FPG", "TC", "CVDPOINTS", "GFRESTIMATE")

### Matching variables
vars_to_match <- c("Age_Group", "Gender", "Race_or_Ethnicity", "Education", 
                   "Smoker", "SBP", "FPG", "TC", "CVDHISTORY", "CVDPOINTS", "SERUMCREAT", "GFRESTIMATE")

#convert variables to factors
for (var in vars_to_convert) {
  sprint_data[[var]] <- as.factor(sprint_data[[var]])
}

####Target Data####
target_data <- read.csv("./Data/Results/NHANES_Gen/NHANES_regen_20k.csv")
target_data$X <- NULL
target_data$ratio <- NULL
###################

add_equity_df <- function(pop, target, weights, var_names, var_subgroups, level, seed, arm_string, container_df){
  ldm_summary <- get_LD_nth_level(pop, target, weights,
                                  var_names, var_subgroups, level) #try 1 and 3 for now
  ldm_summary$Seed <- seed
  ldm_summary$Population <- arm_string
  
  return(ldm_summary)
}


#### Run Main Simulation #####
run_main_simulation <- function(num_seed, pop){
  
  TR_summary <- data.frame(matrix(ncol = 6, nrow = 0))
  Only_CC_summary <- data.frame(matrix(ncol = 6, nrow = 0))
  #Only_SYN_summary <- data.frame(matrix(ncol = 6, nrow = 0))
  Only_Prop_summary <- data.frame(matrix(ncol = 6, nrow = 0))
  Only_Equity_summary <- data.frame(matrix(ncol = 6, nrow = 0))
  Both_summary <- data.frame(matrix(ncol = 6, nrow = 0))
  
  x <- c("Subgroup", "Background_Rate", "Observed_Rate", "LD", "Seed", "Population")
  
  colnames(TR_summary) <- x
  colnames(Only_CC_summary) <- x
  #colnames(Only_SYN_summary) <- x
  colnames(Only_Prop_summary) <- x
  colnames(Only_Equity_summary) <- x
  colnames(Both_summary) <- x
  
  for (i in 1:num_seed){
    
    zero_cc <- F
    zero_sc <- F 
    
    cat("\nIteration:", i, "\n")
    
    #Step 1: create data
    all_data <- partition_data_generic(pop, partition_seed=i, import_target = F)
    target_data_count <- read.csv("./Data/Results/NHANES_Gen/NHANES_count_level3.csv")
    
    if (read_target==T){
      RCT <- data.frame(all_data[[1]][[1]])
      BIASED_EC <- data.frame(all_data[[1]][[2]])
    } else{
      RCT <- data.frame(all_data[[1]])
      BIASED_EC <- data.frame(all_data[[2]])
    }
    remove(all_data)
    
    if(CC_Size == 0){
      zero_cc <- T 
    }
    if(SC_Size == 0){
      zero_sc <- T
    }
    
    cat("Data Partitioned...\n")
    
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
      
      cat(nrow(synthetic_controls), "Synthetic Controls borrowed...\n", sep = " ")
      
      hybrid_RCT <- rbind(RCT, synthetic_controls)
    }
    
    else{ #when zero_sc == True, that means no synthetic controls borrowed, all CC
      hybrid_RCT <- RCT
    }
    
    #### Check if subgroups are all there before IPF adjustment
    columns_to_check <- c("Age_Group", "Gender", "Race_or_Ethnicity")
    unique_counts_to_check <- c(2, 2, 5)
    flag1 <- check_subgroup_counts(RCT %>% filter(RANDASSIGN == 0), columns_to_check, unique_counts_to_check) # only CC
    flag2 <- check_subgroup_counts(hybrid_RCT %>% filter(RANDASSIGN == 0), columns_to_check, unique_counts_to_check) # SC + SC
    
    if ((zero_cc==F) && (flag1==F || flag2==F)){
      cat("Missing subgroup found in Age/Race/Gender with Non-Zero CC \n", sep=" ")
      next
    }
    else if((zero_cc==T) && (flag2==F)){
      cat("Missing subgroup found in Age/Race/Gender with Zero CC \n", sep=" ")
      next
    }
    
    #Step 3: Calculate Adjusted Weights
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
    
    #Step 4: Measure Adjusted Equity
    
    #pop, target, weights, var_names, var_subgroups, level, seed, arm_string, container_df
    #------> 4.1 Adjusted Treatment Summary
    TR_ldm_summary <- add_equity_df(pop= RCT %>% filter(RANDASSIGN == 1), 
                                    target = target_data_count, weights = tr_adjusted_weights,
                                    var_names = var_name_list, var_subgroups = var_subgroup_list,
                                    level = 1, seed = i, arm_string = "Adj_TA", container_df = TR_summary)
    TR_summary <- rbind(TR_summary, TR_ldm_summary)
    
    #------> 4.2 Calculate equity of Only CC
    Only_CC_ldm_summary <- add_equity_df(pop= RCT %>% filter(RANDASSIGN == 0), 
                                         target = target_data_count, weights = NULL,
                                         var_names = var_name_list, var_subgroups = var_subgroup_list,
                                         level = 1, seed = i, arm_string = "Only_CC", container_df = Only_CC_summary)
    Only_CC_summary <- rbind(Only_CC_summary, Only_CC_ldm_summary)
    
    #------> 4.3 Calculate equity of Only Equity
    if(zero_cc==F){
      Only_Equity_ldm_summary <- add_equity_df(pop= RCT %>% filter(RANDASSIGN == 0), 
                                               target = target_data_count, weights = cc_adjusted_weights,
                                               var_names = var_name_list, var_subgroups = var_subgroup_list,
                                               level = 1, seed = i, arm_string = "Only_Equity", container_df = Only_Equity_summary)
    }
    else{
      Only_Equity_ldm_summary <- add_equity_df(pop= RCT %>% filter(RANDASSIGN == 0), 
                                               target = target_data_count, weights = NULL,
                                               var_names = var_name_list, var_subgroups = var_subgroup_list,
                                               level = 1, seed = i, arm_string = "Only_Equity", container_df = Only_Equity_summary)
    }
    
    Only_Equity_summary <- rbind(Only_Equity_summary, Only_Equity_ldm_summary)
    
    #------> 4.4 Calculate equity of Only Propensity
    Only_Prop_ldm_summary <- add_equity_df(pop= hybrid_RCT %>% filter(RANDASSIGN == 0), 
                                           target = target_data_count, weights = NULL,
                                           var_names = var_name_list, var_subgroups = var_subgroup_list,
                                           level = 1, seed = i, arm_string = "Only_Prop", container_df = Only_Prop_summary)
    Only_Prop_summary <- rbind(Only_Prop_summary, Only_Prop_ldm_summary)
    
    #------> 4.5 Calculate equity of Both
    Both_ldm_summary <- add_equity_df(pop= hybrid_RCT %>% filter(RANDASSIGN == 0), 
                                      target = target_data_count, weights = hybrid_adjusted_weights,
                                      var_names = var_name_list, var_subgroups = var_subgroup_list,
                                      level = 1, seed = i, arm_string = "Both", container_df = Both_summary)
    Both_summary <- rbind(Both_summary, Both_ldm_summary)
    cat("Adjusted equity estimated...\n", sep = " ")
    
  }
  
  final_df <- do.call("rbind", list(TR_summary, Only_CC_summary, Only_Prop_summary, Only_Equity_summary, Both_summary))
  return (final_df)
}


CC_Size_list <- c(1000, 0, 500, 1500, 2000)
for (CC_Size in CC_Size_list){
  
  #### STUDY PARAMETERS ####
  IPF_maxiter <- 100
  num_seed <- 50
  scale_to <- 100
  TA_Size <- 2000
  read_target <- F
  SC_Size <- TA_Size - CC_Size
  weight_method <- "ipf" #"gen" or "ipf"
  
  var_name_list <- c('Gender', 'Age_Group', 'Race_or_Ethnicity')
  var_subgroup_list <- list(c('Male','Female'),
                            c('40-59','59+'),
                            c('NH White', 'NH Black', 'NH Asian', 'Hispanic', 'Other'))
  
  eq_str <- "EQUITY_new"
  
  ################################## Trial: Tr vs CH ####################################
  
  combination <- "TR1"
  cat("\nTrial:", combination,"vs Controls \n")
  cat("Starting simulation for TA =", TA_Size, "and CC =", CC_Size,"\n", sep = " ")
  population <- sprint_data
  
  start <- Sys.time()
  summary1 <- run_main_simulation(num_seed = num_seed, pop = population)
  end <- Sys.time()
  paste("Total time elapsed:", end-start)
  summary1$combination <- combination
  
  ####################### SAVE Results ####################
  
  file_string_variable <- paste(TA_Size,"_",CC_Size,".csv", sep = "")
  cat("Saving File", file_string_variable, "...\n")
  
  write.csv(summary1,
            paste('./Data/Results/M6/New_SPRINT_Results/equity/Matchit3/',eq_str,'_full_', file_string_variable, sep = ""),
            row.names=F)
  
}




