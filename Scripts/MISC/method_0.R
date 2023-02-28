setwd("/home/neehan/data/Nafis/ESCA_Primary")
source("./Modules/p_val_add.R")
source("./Modules/equity_metrics_miao.R")
source("./Scripts/common.R")

#Save ATT for next bootstraps
ATT_list <- list()
LDM_summary_list <- list()

total_success <- 50
success_counter <- 1
i <- 1

get_cox_effect <- function(population, weights){
  if (!is.null(weights)){
    fit <- coxph(Surv(EVENTDAYS, EVENT)~RANDASSIGN + Age_Group + Gender + Race_or_Ethnicity + Education +
                   Smoker + SBP + FPG + TC + CVDHISTORY + CVDPOINTS + SERUMCREAT +
                   GFRESTIMATE, data=population, weights = weights)
  }
  else{
    fit <- coxph(Surv(EVENTDAYS, EVENT)~RANDASSIGN + Age_Group + Gender + Race_or_Ethnicity + Education +
                   Smoker + SBP + FPG + TC + CVDHISTORY + CVDPOINTS + SERUMCREAT +
                   GFRESTIMATE, data=population)
  }
  estimate <- coef(fit)['RANDASSIGN']
  return (estimate)
}

while(success_counter <= total_success){
  
  #step: a
  TA <- datap %>% filter(RANDASSIGN==1)
  CC <- datap %>% filter(RANDASSIGN==0)
  HC <- CC
  
  #step: e
  W_TA_ipf <- get_IPF_weights(dat = TA, maxIter = 1000)
  W_HC_ipf <- get_IPF_weights(dat = HC, maxIter = 1000)
  
  #step: f
  set.seed(i)
  TA_IPF_sample <- TA[sample(seq_len(nrow(TA)), size = nrow(TA), replace = TRUE, prob = W_TA_ipf),]
  
  set.seed(i)
  HC_IPF_sample <- HC[sample(seq_len(nrow(HC)), size = nrow(HC), replace = TRUE, prob = W_HC_ipf),]
  
  #step: g
  merged_sample <- rbind(TA_IPF_sample, HC_IPF_sample)
  
  #run linear model to get the treatment effect
  ATT <- get_cox_effect(merged_sample, weights = NULL)
  
  #store treatment effect from each run
  ATT_list <- c(ATT_list, ATT)
  
  #print iteration end announcement
  cat(paste0("Iteration ", i, " done!! \n"))
  
  i <- i + 1
  success_counter <- success_counter + 1
  
}

################### Result Presentation ######################
print("------------------------------------------------")
cat(paste0("Estimated Mean is ",mean(unlist(ATT_list)), "\n", 
           "Variance is ", var(unlist(ATT_list)), "\n", 
           "Standard Deviation is ", sd(unlist(ATT_list)), "\n"))

write.csv(as.data.frame(unlist(ATT_list)) %>% rename(PATT_ground = 'unlist(ATT_list)'), 
          "./Data/Results/PATT_ground50_adj.csv")
