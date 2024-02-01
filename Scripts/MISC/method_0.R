setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
source("./Modules/p_val_add.R")
source("./Modules/equity_metrics_miao.R")
source("./Scripts/common.R")

#Save ATT for next bootstraps
ATT_list <- list()
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

while(success_counter <= total_success){
  
  #step: a
  TA <- data %>% filter(RANDASSIGN==1)
  CC <- data %>% filter(RANDASSIGN==0)
  
  set.seed(i)
  #TA_sample <- sample_frac(TA, 0.95)
  TA_sample <- sample_n(TA, size = nrow(TA), replace = T)
  
  set.seed(i)
  #CC_sample <- sample_frac(CC, 0.95)
  CC_sample <- sample_n(CC, size=nrow(CC), replace=T)
  
  HC_sample <- CC_sample

  #step: e
  W_TA_ipf <- get_IPF_weights(dat = TA_sample, maxIter = 2000)
  W_HC_ipf <- get_IPF_weights(dat = HC_sample, maxIter = 2000)

  #step: f
  set.seed(i)
  TA_IPF_sample <- TA_sample[sample(seq_len(nrow(TA_sample)), size = nrow(TA_sample), 
                                    replace = TRUE, prob = W_TA_ipf),]

  set.seed(i)
  HC_IPF_sample <- HC_sample[sample(seq_len(nrow(HC_sample)), size = nrow(HC_sample), 
                                    replace = TRUE, prob = W_HC_ipf),]
  
  #step: g
  merged_sample <- rbind(TA_IPF_sample, HC_IPF_sample)
  #merged_sample <- rbind(TA_sample, CC_sample)
  
  #run linear model to get the treatment effect
  ATT <- get_cox_effect(merged_sample, weights = NULL)
  
  #store treatment effect from each run
  ATT_list <- c(ATT_list, ATT)
  
  #print iteration end announcement
  cat(paste0("Iteration ", i, " done!! \n"))
  
  i <- i + 1
  success_counter <- success_counter + 1
  
}

####### Note #######

#-0.225 [-0.247, -0.202] PATT Estimate Unadj
#-0.295 [-0.302, -0.287] SATT Estimate Unadj

hr_list <- lapply(ATT_list, exp)

################### Result Presentation ######################
print("------------------------------------------------")
cat(paste0("Estimated Mean is ",mean(unlist(hr_list)), "\n", 
           "Variance is ", var(unlist(hr_list)), "\n", 
           "Standard Deviation is ", sd(unlist(hr_list)), "\n",
           "Confidence Interval is", confidence_interval(unlist(hr_list), 0.95)))

write.csv(as.data.frame(unlist(ATT_list)) %>% rename(PATT_ground = 'unlist(ATT_list)'), 
          "./Data/Results/PATT_ground500.csv") 
