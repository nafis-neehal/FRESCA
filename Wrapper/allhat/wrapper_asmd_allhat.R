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
  
  seeds_used <- c()
  
  #### ADDED ASMD CALCULATION ###
  columns <- c("Seed","Combination","Variable", "SMD") 
  asmd_df <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
  colnames(asmd_df) <- columns
  
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
    
    ## DD Snapshot TA vs Only CC ## <---------
    new_asmd_df <- get_ASMD_new(pop1 = RCT %>% filter(RANDASSIGN == 1),
                            pop2 = RCT %>% filter(RANDASSIGN == 0),
                            seed = i,
                            controls = "only_cc",
                            #var_list = colnames(RCT %>% select(Age_Group:BMI)),
                            var_list = colnames(RCT %>% select(vars_to_match)),
                            control_weights = rep(1, nrow(RCT %>% filter(RANDASSIGN == 0))))
    
    asmd_df <- rbind(asmd_df, new_asmd_df)
    ###############################
    
    ## DD Snapshot TA vs Only Biased EC ## <---------
    new_asmd_df <- get_ASMD_new(pop1 = RCT %>% filter(RANDASSIGN == 1),
                            pop2 = BIASED_EC,
                            seed = i,
                            controls = "only_biased_ec",
                            #var_list = colnames(RCT %>% select(Age_Group:BMI)),
                            var_list = colnames(RCT %>% select(vars_to_match)),
                            control_weights = rep(1, nrow(BIASED_EC)))
    
    asmd_df <- rbind(asmd_df, new_asmd_df)
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
      
      ## DD Snapshot TA vs Only SC ## <---------
      new_asmd_df <- get_ASMD_new(pop1 = RCT %>% filter(RANDASSIGN == 1),
                              pop2 = synthetic_controls,
                              seed = i,
                              controls = "only_sc",
                              #var_list = colnames(RCT %>% select(Age_Group:BMI)),
                              var_list = colnames(RCT %>% select(vars_to_match)),
                              control_weights = rep(1, nrow(synthetic_controls)))
      
      asmd_df <- rbind(asmd_df, new_asmd_df)
      ###############################
      
      cat(nrow(synthetic_controls), "Synthetic Controls borrowed...\n", sep = " ")
      
      hybrid_RCT <- rbind(RCT, synthetic_controls)
    }
    
    else{ #when zero_sc == True, that means no synthetic controls borrowed, all CC
      
      new_asmd_df <- data.frame(Seed        = i,
                                Controls    = "only_sc",
                                #var_list = colnames(RCT %>% select(Age_Group:BMI)),
                                Variable = colnames(RCT %>% select(vars_to_match)),
                                SMD         = NA)
      asmd_df <- rbind(asmd_df, new_asmd_df)
      
      hybrid_RCT <- RCT
    }
    
    ## DD Snapshot TA vs Hybrid CC, SC ## <---------
    new_asmd_df <- get_ASMD_new(pop1 = RCT %>% filter(RANDASSIGN == 1),
                            pop2 = hybrid_RCT %>% filter(RANDASSIGN == 0),
                            seed = i,
                            controls = "only_prop",
                            #var_list = colnames(RCT %>% select(Age_Group:BMI)),
                            var_list = colnames(RCT %>% select(vars_to_match)),
                            control_weights = rep(1, nrow(hybrid_RCT %>% filter(RANDASSIGN == 0))))
    
    asmd_df <- rbind(asmd_df, new_asmd_df)
    ###############################
    
    #### Check if subgroups are all there before IPF adjustment
    columns_to_check <- c("Age_Group", "Gender", "Race_or_Ethnicity")
    unique_counts_to_check <- c(2, 2, 5)
    flag1 <- check_subgroup_counts(RCT %>% filter(RANDASSIGN == 0), columns_to_check, unique_counts_to_check) # only CC
    flag2 <- check_subgroup_counts(hybrid_RCT %>% filter(RANDASSIGN == 0), columns_to_check, unique_counts_to_check) # SC + SC
    
    if ((zero_cc==F) && (flag1==F || flag2==F)){
      asmd_df <- asmd_df %>% filter(Seed!=i) #rollback without adding asmd seeds of this iteration
      cat("Missing subgroup found in Age/Race/Gender with Non-Zero CC \n", sep=" ")
      next
    }
    else if((zero_cc==T) && (flag2==F)){
      asmd_df <- asmd_df %>% filter(Seed!=i) #rollback without adding asmd seeds of this iteration
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
    
    #Step 4: Measure Adjusted SMD
    ## DD Snapshot TA vs Weight Adjusted CC ## <---------
    if(zero_cc==F){
      # adjusted_CC_pop <- sample_n(RCT %>% filter(RANDASSIGN==0), size = nrow(RCT %>% filter(RANDASSIGN==0)), 
      #                             weight = cc_adjusted_weights, replace = T)
      
      new_asmd_df <- get_ASMD_new(pop1 = RCT %>% filter(RANDASSIGN == 1),
                              pop2 = RCT %>% filter(RANDASSIGN==0),
                              seed = i,
                              controls = "only_equity",
                              #var_list = colnames(RCT %>% select(Age_Group:BMI)),
                              var_list = colnames(RCT %>% select(vars_to_match)),
                              control_weights = cc_adjusted_weights)
      
      asmd_df <- rbind(asmd_df, new_asmd_df)
    }
    else{ #when no cc, doing this will auto check cc size and create NULL rows
      new_asmd_df <- get_ASMD_new(pop1 = RCT %>% filter(RANDASSIGN == 1),
                              pop2 = RCT %>% filter(RANDASSIGN == 0), ### passing only CC here as a bypass, automatically will be managed inside
                              seed = i,
                              controls = "only_equity",
                              #var_list = colnames(RCT %>% select(Age_Group:BMI)),
                              var_list = colnames(RCT %>% select(vars_to_match)),
                              control_weights = rep(1, nrow(RCT %>% filter(RANDASSIGN == 0))))
      
      asmd_df <- rbind(asmd_df, new_asmd_df)
    }
    ###############################
    
    ## DD Snapshot TA vs Weight Adjusted Hybrid ## <---------
    # adjusted_HC_pop <- sample_n(hybrid_RCT %>% filter(RANDASSIGN==0), size = nrow(hybrid_RCT %>% filter(RANDASSIGN==0)), 
    #                             weight = hybrid_adjusted_weights, replace = T)
    
    new_asmd_df <- get_ASMD_new(pop1 = RCT %>% filter(RANDASSIGN == 1),
                            pop2 = hybrid_RCT %>% filter(RANDASSIGN==0),
                            seed = i,
                            controls = "both",
                            #var_list = colnames(RCT %>% select(Age_Group:BMI)),
                            var_list = colnames(RCT %>% select(vars_to_match)),
                            control_weights = hybrid_adjusted_weights)
    
    asmd_df <- rbind(asmd_df, new_asmd_df)
    ###############################
    
    
    cat("All SMD Calculated...\n", sep = " ")
    
  }
  
  return (asmd_df)
}

######### mini-helpers #######

get_results_df <- function(TA_Size, CC_Size, population, combi_str){
  
  start <- Sys.time()
  
  result <- run_main_simulation(num_seed=num_seed, pop=population)
  result_asmd <- result
  result_asmd$combination <- combi_str
  
  end <- Sys.time()
  
  paste("Total time elapsed:", end-start)
  
  return(result_asmd)
}

############################ Global: Big Umbrella in Loop Varying CC Size #############

#CC_Size_list <- c(2000, 0, 500, 1000, 4000)
#CC_Size_list <- c(0, 1000, 2000, 3000, 4000, 5000)

TA_Size_list <- c(7000, 8000, 7500, 6500, 6000, 5500, 5000)
CC_Size_list <- c(3000, 2000, 2500, 3500, 4000, 4500, 5000)
SC_Size_list <- c(4000, 6000, 5000, 3000, 2000, 1000, 0)


# for (CC_Size in CC_Size_list){
for (i in seq(1, length(TA_Size_list), 1)){
  
  #### STUDY PARAMETERS ####
  IPF_maxiter <- 100
  num_seed <- 50
  scale_to <- 100
  #TA_Size <- 5000
  #SC_Size <- TA_Size - CC_Size
  TA_Size <- TA_Size_list[i]
  CC_Size <- CC_Size_list[i]
  SC_Size <- SC_Size_list[i]
  
  weight_method <- "ipf" #"gen" or "ipf"
  
  ################################## Trial 1: Li vs CH ####################################
  
  combination <- "TR1"
  cat("\nTrial:", combination,"vs Controls \n")
  cat("Starting simulation for TA =", TA_Size, "and CC =", CC_Size,"\n", sep = " ")
  population <- g1
  result1 <- get_results_df(TA_Size, CC_Size, population, combination) 
  result1_asmd <- result1
  
  ################################## Trial 2: Li vs CH ####################################
  
  # combination <- "TR2"
  # cat("\nTrial:", combination,"vs Controls \n")
  # cat("Starting simulation for TA =", TA_Size, "and CC =", CC_Size,"\n", sep = " ")
  # population <- g2
  # result2 <- get_results_df(TA_Size, CC_Size, population, combination) 
  # result2_asmd <- result2
  
  ####################### SAVE Results ####################
  
  file_string_variable <- paste(TA_Size,"_",CC_Size,".csv", sep = "")
  cat("Saving File", file_string_variable, "...\n")
  
  # write.csv(rbind(result1_asmd),
  #           paste('./Data/Results/M6/ALLHAT_Results/ASMD_Final/Matchit3/ASMD_new_weighted_', file_string_variable, sep= ""),
  #           row.names=F)
  write.csv(rbind(result1_asmd),
            paste('./Data/Results/M6/ALLHAT_Results/Samecost/trial_size_10k/ASMD/ASMD_new_weighted_', file_string_variable, sep= ""),
            row.names=F) 
  
}

##############################
# controls_order <- c("only_cc", "only_biased_ec", "only_sc", "hybrid")
# new_df <- ASMD_CVD_4k2k %>% group_by(combination, Controls) %>% 
#   summarise(Mean     = mean(SMD), 
#             CI_Lower = CI(SMD)[3],
#             CI_Upper = CI(SMD)[1],
#             .groups='drop') %>% 
#   mutate(Controls = factor(Controls, levels = controls_order)) %>%
#   arrange(combination, Controls) %>%
#   as.data.frame()
# 
# ASMD_CVD_4k2k
##############################

