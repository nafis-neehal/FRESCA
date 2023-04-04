setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
source("./Modules/p_val_add.R")
source("./Modules/equity_metrics_miao.R")
source("./Scripts/common.R")

#Save ATT for next bootstraps
LDM_CC_list <- list()
LDM_summary_list <- list()

total_success <- 50
success_counter <- 1
i <- 1

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

LDM_Summary <- data.frame(matrix(ncol = 4, nrow = 0))
colsLDM <- c("Seed", "Var", "Level","LDM_CC_IPF")
colnames(LDM_Summary) <- colsLDM


while(success_counter <= total_success){
  
  #step: a
  TA <- data %>% filter(RANDASSIGN==1)
  CC <- data %>% filter(RANDASSIGN==0)
  
  # set.seed(i)
  # TA_sample <- sample_frac(TA, 0.90)
  # 
  # set.seed(i)
  # CC_sample <- sample_frac(CC, 0.90)
  
  HC <- CC

  #step: e
  W_TA_ipf <- get_IPF_weights(dat = TA, maxIter = 1000)
  W_HC_ipf <- get_IPF_weights(dat = HC, maxIter = 1000)

  #step: f
  set.seed(i)
  TA_IPF_sample <- TA[sample(seq_len(nrow(TA)), size = nrow(TA), replace = TRUE, prob = W_TA_ipf),]

  set.seed(i)
  HC_IPF_sample <- HC[sample(seq_len(nrow(HC)), size = nrow(HC), replace = TRUE, prob = W_HC_ipf),]

  LDM_CC_three <- get_LDM_from_count(HC_IPF_sample) %>% select(VAR, Level, LDM)
  VAR <- LDM_CC_three$VAR
  Level <- LDM_CC_three$Level
  LDM_CC <- LDM_CC_three %>% select(LDM) %>% rename(LDM_CC = LDM)
  
  row <- as.data.frame(list(success_counter, VAR, Level, LDM_CC))
  colnames(row) <- colsLDM
  LDM_Summary <- rbind(LDM_Summary, row)
  
  
  #step: g
  #merged_sample <- rbind(TA_IPF_sample, HC_IPF_sample)
  # merged_sample <- rbind(TA_sample, CC_sample)
  # 
  # #run linear model to get the treatment effect
  # ATT <- get_cox_effect(merged_sample, weights = NULL)
  # 
  # #store treatment effect from each run
  # ATT_list <- c(ATT_list, ATT)
  
  #print iteration end announcement
  cat(paste0("Iteration ", i, " done!! \n"))
  
  i <- i + 1
  success_counter <- success_counter + 1
  
}

mean_LDM_Summary <- LDM_Summary %>% group_by(Level) %>% summarise(median(abs(LDM_CC_IPF)))


####### Note #######

#-0.225 [-0.247, -0.202] PATT Estimate Unadj
#-0.295 [-0.302, -0.287] SATT Estimate Unadj

################### Result Presentation ######################
# print("------------------------------------------------")
# cat(paste0("Estimated Mean is ",mean(unlist(ATT_list)), "\n", 
#            "Variance is ", var(unlist(ATT_list)), "\n", 
#            "Standard Deviation is ", sd(unlist(ATT_list)), "\n",
#            "Confidence Interval is", confidence_interval(unlist(ATT_list), 0.95)))

write.csv(as.data.frame(LDM_Summary),
          "./Data/Results/PLDM_ground50.csv")
