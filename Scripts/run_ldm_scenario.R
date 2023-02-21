setwd("/home/neehan/data/Nafis/ESCA_Git")
source("./Scripts/common.R")
source("./Modules/p_val_add.R")
source("./Modules/equity_metrics_miao.R")

TA_World <- data %>% filter(RANDASSIGN==1)
TA_size <- 1000
total_seeds <- 20
IPF_maxiter <- 100
step_size <- 100
randomization_ratio <- 1
num_col <- 10
CC_size_list <- seq(TA_size-step_size, 100, -step_size)
EC_size_list <- randomization_ratio*TA_size - CC_size_list

confidence_interval <- function(vector, interval) {
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

get_propensity_model <- function(population){
  m_ps <- glm(inTrial ~ Age_Group + Gender + Race_or_Ethnicity + Education +
                Smoker + FPG + TC + CVDHISTORY + CVDPOINTS + SERUMCREAT +
                GFRESTIMATE, data = population, family = binomial())
  return (m_ps)
}

get_propensity_score <- function(Model, population){
  prs_df <- data.frame(MASKID = population %>% select(MASKID), 
                       RANDASSIGN = population %>% select(RANDASSIGN), 
                       pr_score = predict(Model, population, type = "response"))
  return (prs_df)
}

get_matched_data <- function(population, pscore){
  p <- population %>% mutate(Age_Group = factor(Age_Group),
                             Gender = factor(Gender),
                             Race_or_Ethnicity = factor(Race_or_Ethnicity),
                             Education = factor(Education),
                             Smoker = factor(Smoker),
                             FPG = factor(FPG),
                             TC = factor(TC),
                             GFRESTIMATE = factor(GFRESTIMATE))
  
  m.out <- matchit(RANDASSIGN ~ Age_Group + Gender + Race_or_Ethnicity + Education +
                     Smoker + FPG + TC + CVDHISTORY + CVDPOINTS + SERUMCREAT +
                     GFRESTIMATE, data=p, distance = pscore, replace = T)
  
  matched_data <- match.data(m.out) 
  
  return(matched_data %>% filter(RANDASSIGN==0)) #contains weights
  
}


run_one_ldm_scenario <- function(seed_value, CC_Size){
  
  LDM_Summary <- data.frame(matrix(ncol = 10, nrow = 0))
  colsLDM <- c("Seed", "CC_Size", "Var", "Level", "LDM_TA", "LDM_CC", "LDM_HC_Prop", "LDM_HC_IPF", 
               "LDM_TA_Prop_IPF", "LDM_HC_Prop_IPF" )
  colnames(LDM_Summary) <- colsLDM
  
  set.seed(seed_value)
  TA <- sample_n(TA_World %>% filter(RANDASSIGN==1), TA_size)
  W_TA_ipf <- get_IPF_weights(dat = TA, maxIter = IPF_maxiter)
  LDM_TA <- get_LDM_from_count(TA) %>% select(LDM) %>% rename(LDM_TA = LDM)
  
  set.seed(seed_value)
  CC_World <- sample_n(data %>% filter(RANDASSIGN==0), TA_size)
  
  external_population <- setdiff(data, rbind(TA_World, CC_World))
  EC_world <- external_population %>% filter(RANDASSIGN==0)
  
  W_EC_bias_ipf <- get_ECWorld_Bias_weights(dat = EC_world, maxIter = IPF_maxiter)
  Biased_EC_world <- EC_world[sample(seq_len(nrow(EC_world)), 
                                     size = nrow(EC_world), 
                                     replace = TRUE, prob = W_EC_bias_ipf),]
  
  #create the propensity score model - propensity to get in the trial
  M <- get_propensity_model(do.call("rbind", list(TA %>% mutate(inTrial = 1),
                                                  CC_World %>% mutate(inTrial = 1),
                                                  Biased_EC_world %>% mutate(inTrial = 0))))
  
  LDM_CC_list <- list()
  LDM_HC_Prop_list <- list()
  LDM_HC_IPF_list <- list()
  LDM_TA_Prop_IPF_list <- list()
  LDM_HC_Prop_IPF_list <- list()
  
  flag <- 1
  n_plus_sum <- 1
  success_counter <-1
  total_success <- 10
  inner_seed <- 1
  
  while (success_counter <= total_success){ #loops 10 times
    
    print(paste("n=", success_counter))
    
    set.seed(inner_seed)
    CC <- sample_n(CC_World %>% filter(RANDASSIGN==0), CC_Size)
    
    #<- Collect LDM_CC 
    LDM_CC_three <- get_LDM_from_count(CC) %>% select(VAR, Level, LDM)
    VAR <- LDM_CC_three$VAR
    Level <- LDM_CC_three$Level 
    
    LDM_CC <- LDM_CC_three %>% select(LDM) %>% rename(LDM_CC = LDM)
    LDM_CC_list <- c(LDM_CC_list, LDM_CC)
    
    for (n_plus_sum in 0:4){ #for each inner_seed value I'll try 5 times to get a pop without missing subcategory
      
      set.seed(inner_seed + n_plus_sum)
      unmatched_EC <- sample_n(Biased_EC_world, randomization_ratio * TA_size - CC_Size)
      
      #get propensity score of this unmatched EC sample
      ps <- get_propensity_score(M, rbind(TA, unmatched_EC)) #contains MASKID, RANDASSIGN and PRS
      
      #get matched EC with weights
      matched_EC_df <- get_matched_data(rbind(TA, unmatched_EC), unlist(ps$pr_score))
      matched_EC <- matched_EC_df %>% select(MASKID:SEATSYS)
      matched_EC_weights <- matched_EC_df$weights #contains only EC weights
      
      
      #create matched and unmatched HC
      unmatched_HC <- rbind(CC, unmatched_EC)
      matched_HC <- rbind(CC, matched_EC) #just contains maskid:seatsys
      matched_HC_weights <- c(rep(1, nrow(CC)), matched_EC_weights)
      matched_HC_df <- matched_HC
      matched_HC_df$weights <- matched_HC_weights
      
      flag <- check_subgroup_counts(age_sbgrps = 2, gender_sbgrps = 2, race_sbgrps = 5, 
                                    gfr_sbgrps = 2, TA = TA, HC = matched_HC)
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
    
    #<- Collect LDM_HC_Prop
    LDM_HC_Prop <- get_LDM_from_count(matched_HC) %>% select(LDM) %>% rename(LDM_HC_Prop = LDM)
    LDM_HC_Prop_list <- c(LDM_HC_Prop_list, LDM_HC_Prop)
    
    #Apply IPF on TA and matched HC
    W_unmatched_HC_ipf <- get_IPF_weights(dat = unmatched_HC, maxIter = IPF_maxiter) #no propensity
    W_matched_HC_ipf <- get_IPF_weights(dat = matched_HC, maxIter = IPF_maxiter) #propensity
    
    #<- Collect TA_HC_IPF c[0]
    #<- Collect TA_HC_both c[1]
    set.seed(inner_seed)
    TA_IPF_sample <- TA[sample(seq_len(nrow(TA)), size = nrow(TA), replace = TRUE, prob = W_TA_ipf),]
    
    set.seed(inner_seed)
    unmatched_HC_IPF_sample <- unmatched_HC[sample(seq_len(nrow(unmatched_HC)), 
                                                   size = nrow(TA), replace = TRUE, prob = W_unmatched_HC_ipf),]
    
    set.seed(inner_seed)
    matched_HC_IPF_sample_df <- matched_HC_df[sample(seq_len(nrow(matched_HC_df)), 
                                                     size = nrow(TA), replace = TRUE, prob = W_matched_HC_ipf),]
    
    matched_HC_IPF_sample <- matched_HC_IPF_sample_df %>% select(MASKID:SEATSYS)
    matched_HC_IPF_sample_weights <- matched_HC_IPF_sample_df$weights
    
    #<- Collect LDM_TA_IPF, LDM_HC_IPF
    LDM_TA_Prop_IPF <- get_LDM_from_count(TA_IPF_sample) %>% select(LDM) %>% rename(LDM_TA_Prop_IPF = LDM)
    LDM_TA_Prop_IPF_list <- c(LDM_TA_Prop_IPF_list, LDM_TA_Prop_IPF)
    
    LDM_HC_IPF <- get_LDM_from_count(unmatched_HC_IPF_sample) %>% select(LDM) %>% rename(LDM_HC_IPF = LDM)
    LDM_HC_IPF_list <- c(LDM_HC_IPF_list, LDM_HC_IPF)
    
    LDM_HC_Prop_IPF <- get_LDM_from_count(matched_HC_IPF_sample) %>% select(LDM) %>% rename(LDM_HC_Prop_IPF = LDM)
    LDM_HC_Prop_IPF_list <- c(LDM_HC_Prop_IPF_list, LDM_HC_Prop_IPF)
    
    inner_seed <- inner_seed + 1
    
  }
  
  #add dataframe
  row <- as.data.frame(list(seed_value, CC_Size, VAR, Level,
                            abs(LDM_TA),
                            abs(rowMeans(as.data.frame(LDM_CC_list))),
                            abs(rowMeans(as.data.frame(LDM_HC_Prop_list))),
                            abs(rowMeans(as.data.frame(LDM_HC_IPF_list))),
                            abs(rowMeans(as.data.frame(LDM_TA_Prop_IPF_list))),
                            abs(rowMeans(as.data.frame(LDM_HC_Prop_IPF_list)))))
  colnames(row) <- colsLDM
  LDM_Summary <- rbind(LDM_Summary, row)
  return(LDM_Summary)
}
