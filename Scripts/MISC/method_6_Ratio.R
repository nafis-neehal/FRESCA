setwd("/Users/nafisneehal/ESCA")
source("./Modules/p_val_add.R")
source("./Modules/equity_metrics_miao.R")
source("./Scripts/common.R")

directory <- "./Data/Results/PATT_m6/CC_EC_Vary_Size/matchCCEC/"
filename <- "1500TA_vary_CC_1_5EC_matchCCEC_seed10.csv"

# Hyperparameters
total_seed <- 10
seed_counter <- 1
seed_value <- 1
TA_size <- 1500
IPF_maxiter <- 2000
randomization_ratio <- 1.5

#Global Variables
df0 <- data.frame(matrix(ncol = 4, nrow = 0))
cols <- c("seed", "EC_Sample_Size", "Mean_Patt", "Var_Patt")
colnames(df0) <- cols

#TA Buffer define
TA_World <- data %>% filter(RANDASSIGN==1)

while(seed_counter <= total_seed) {
  
  print("------------------------------------------------")
  cat(paste0("Seed has been set to ", seed_value, " at Level 0\n"))
  cat(paste0("Scenario No: ", seed_counter, "\n"))
  
  #Save ATT in inner loop1
  ATT_list_loop1 <- list()
  LDM_summary_list_loop1 <- list()
  VAR_list_loop1 <- list()
  
  LDM_CC_loop1 <- list()
  LDM_HC_loop1 <- list()
  
  LDM_IPF_TA_loop1 <- list()
  LDM_IPF_HC_loop1 <- list()
  
  success_counter <- 1
  i <- 1
  s <- 1
  
  #step: 6a
  set.seed(seed_value)
  CC_World <- sample_n(data %>% filter(RANDASSIGN==0), TA_size)
  
  set.seed(seed_value)
  TA <- sample_n(TA_World %>% filter(RANDASSIGN==1), TA_size) #TA is set for loop 2 and 3
  
  external_population <- setdiff(data, rbind(TA_World, CC_World))
  EC_world <- external_population %>% filter(RANDASSIGN==0)
  
  #check for various CC size------------------>>>>>>>>>
  CC_size_list <- seq(TA_size, 0, -50)
  EC_size_list <- randomization_ratio*TA_size - CC_size_list
  total_success <- length(EC_size_list)
  
  #step: 6a(i)
  W_EC_ipf <- get_ECWorld_Bias_weights(dat = EC_world, maxIter = IPF_maxiter)
  EC_sample_size_list_loop1 <- list()
  CC_sample_size_list_loop1 <- list()
  
  flag <- 1
  
  while(success_counter <= total_success){
    
    print("---------------------------------------->")
    cat(paste0("Seed has been set to ", i, " at Level 1\n"))
    
    #Vary CC, EC Sample size <- Loop 1
    set.seed(i)
    CC <- sample_n(CC_World %>% filter(RANDASSIGN==0), CC_size_list[s]) #CC is set for loop 3
    
    set.seed(i)
    EC_ipf_sample <- EC_world[sample(seq_len(nrow(EC_world)), size = EC_size_list[s], 
                                     replace = TRUE, prob = W_EC_ipf),] #EC is set for loop 3
    
    #step: 6a(i)(1, 2)
    HC <- get_Matched_HC_from_CC_EC(TA, CC, EC_ipf_sample)
    
    #Need a checker to check for missing categories
    flag <- check_subgroup_counts(age_sbgrps = 2, gender_sbgrps = 2, race_sbgrps = 5, 
                                  TA = TA, HC = HC)
    if(flag==0 & i<=10){
      #end_iteration <- end_iteration + 1
      cat(paste0("Iteration ", i, " aborted due to missing categories. Continuting to next iteration. \n"))
      i <- i+1
      next
    }
    else if(flag==0 & i>10){
      seed_value <- seed_value + 1
      break
    }
    
    #step: 6a(i)(2)
    #HC <- rbind(CC, EC)
    
    #step: 6a(i)(3)
    W_TA_ipf <- get_IPF_weights(dat = TA, maxIter = IPF_maxiter)
    W_HC_ipf <- get_IPF_weights(dat = HC, maxIter = IPF_maxiter)
    
    #Save ATT in inner loop2
    ATT_list_loop2 <- list()
    LDM_summary_list_loop2 <- list()
    
    #step: 6a(i)(4)
    print("--------------------------------->")
    for (n in 1:20) {
      
      set.seed(n)
      TA_IPF_sample <- TA[sample(seq_len(nrow(TA)), size = nrow(TA), replace = TRUE, prob = W_TA_ipf),]
      
      set.seed(n)
      HC_IPF_sample <- HC[sample(seq_len(nrow(HC)), size = nrow(TA), replace = TRUE, prob = W_HC_ipf),]
      
      #step: 6a(i)(4)(a)
      merged_sample <- rbind(TA_IPF_sample, HC_IPF_sample)
      
      m.out <- matchit(RANDASSIGN ~ Age_Group + Gender + Race_or_Ethnicity + Education +
                         Smoker + SBP + FPG + TC + CVDHISTORY + CVDPOINTS + SERUMCREAT +
                         GFRESTIMATE, data=merged_sample, method = "nearest", distance = "glm")
      matched_data <- match.data(m.out)
      
      #run linear model to get the treatment effect
      fit <- lm(SEATSYS ~ RANDASSIGN + Age_Group + Gender + Race_or_Ethnicity + Education + 
                  Smoker + SBP + FPG + TC + CVDHISTORY + CVDPOINTS + SERUMCREAT +
                  GFRESTIMATE, data = matched_data, weights = weights)
      
      ##### result calculation in each iteration
      #calculate ATT
      ATT <- coeftest(fit, vcov. = vcovCL, cluster = ~subclass)["RANDASSIGN",,drop = FALSE][1]
      
      #store treatment effect from each run
      ATT_list_loop2 <- c(ATT_list_loop2, ATT)
      
      #store LDM from each run for LDM - TA, CC, EC_Candidate, HC, IPF_TA, IPF_HC
      LDM_summary <- get_LDM_summary(TA, CC, EC_ipf_sample, HC, TA_IPF_sample, HC_IPF_sample)
      
      #store LDM summary for each round as a list item
      LDM_summary_list_loop2 <- c(LDM_summary_list_loop2, list(LDM_summary))
      
      #print iteration end announcement
      #cat(paste0("Iteration ", n, " done in Level 2!! \n"))
      
    }
    
    ################### Result Presentation ######################
    print("---------------------------------------->")
    cat(paste0("Results from Level 1: \n"))
    cat(paste0("Estimated Mean is ",mean(unlist(ATT_list_loop2)), "\n", 
               "Variance is ", var(unlist(ATT_list_loop2)), "\n", 
               "Standard Deviation is ", sd(unlist(ATT_list_loop2)), "\n"))
    
    ldm_sum <- as.data.frame(LDM_summary_list_loop2[1]) %>% select(LDM_TA:LDM_IPF_HC) %>% abs()
    if (length(LDM_summary_list_loop2)>=2){
      for (ldm in LDM_summary_list_loop2[2:length(LDM_summary_list_loop2)]){
        ldm <- as.data.frame(ldm)
        ldm_dat <- ldm %>% select(LDM_TA:LDM_IPF_HC) %>% abs()
        ldm_sum <- ldm_sum + ldm_dat
      }
    }
    
    ldm_mean <- ldm_sum/length(LDM_summary_list_loop2)
    ldm_mean <- cbind(LDM_summary%>%select(VAR, Level), ldm_mean)
    
    
    ldm_mean2_cc <- ldm_mean %>% select(LDM_CC) %>% abs() %>% sum() / length(ldm_mean)
    ldm_mean2_hc <- ldm_mean %>% select(LDM_HC) %>% abs() %>% sum() / length(ldm_mean)
    ldm_mean2_ipf_ta <- ldm_mean %>% select(LDM_IPF_TA) %>% abs() %>% sum() / length(ldm_mean)
    ldm_mean2_ipf_hc <- ldm_mean %>% select(LDM_IPF_HC) %>% abs() %>% sum() / length(ldm_mean)
    
    #this is what gets saved for one EC Sample Size and one seed
    ATT_list_loop1 <- c(ATT_list_loop1, mean(unlist(ATT_list_loop2)))
    LDM_summary_list_loop1 <- c(LDM_summary_list_loop1, list(ldm_mean))
    VAR_list_loop1 <- c(VAR_list_loop1, var(unlist(ATT_list_loop2)))
    EC_sample_size_list_loop1 <- c(EC_sample_size_list_loop1, EC_size_list[s])
    CC_sample_size_list_loop1 <- c(CC_sample_size_list_loop1, CC_size_list[s])
    
    LDM_CC_loop1 <- c(LDM_CC_loop1, ldm_mean2_cc)
    LDM_HC_loop1 <- c(LDM_HC_loop1, ldm_mean2_hc)
    LDM_IPF_TA_loop1 <- c(LDM_IPF_TA_loop1, ldm_mean2_ipf_ta)
    LDM_IPF_HC_loop1 <- c(LDM_IPF_HC_loop1, ldm_mean2_ipf_hc)
    
    i <- i + 1
    success_counter <- success_counter + 1
    s <- s+1
    
  }
  
  if(flag == 0){
    seed_value <- seed_value + 1
    next
  }
  
  df_ATT_summary_loop1 <- data.frame(seed = seed_value,
                                     EC_Sample_Size = unlist(EC_sample_size_list_loop1),
                                     CC_Sample_Size = unlist(CC_sample_size_list_loop1),
                                     Mean_Patt = unlist(ATT_list_loop1),
                                     Var_Patt = unlist(VAR_list_loop1),
                                     LDM_CC = unlist(LDM_CC_loop1),
                                     LDM_HC = unlist(LDM_HC_loop1),
                                     LDM_IPF_TA = unlist(LDM_IPF_TA_loop1),
                                     LDM_IPF_HC = unlist(LDM_IPF_HC_loop1))
  
  df0 <- rbind(df0, df_ATT_summary_loop1)
  
  print(df_ATT_summary_loop1)
  
  seed_counter <- seed_counter + 1
  seed_value <- seed_value + 1
  
}

#print result in console
df0 %>% group_by(EC_Sample_Size) %>% summarise(mean(Mean_Patt))

#Save results
write.csv(df0, paste(directory, filename, sep = ""), row.names = FALSE)


