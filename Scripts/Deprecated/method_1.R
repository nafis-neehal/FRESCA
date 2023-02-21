setwd("/Users/nafisneehal/ESCA")
source("./Modules/p_val_add.R")
source("./Modules/equity_metrics_miao.R")
source("./Scripts/common.R")


#Save ATT for next bootstraps
ATT_list <- list()
LDM_summary_list <- list()

total_success <- 50
success_counter <- 1
i <- 1

while(success_counter <= total_success){
  
  set.seed(i)
  
  #step: a
  TA <- sample_n(data %>% filter(RANDASSIGN==1), 1000)
  CC <- sample_n(data %>% filter(RANDASSIGN==0), 500)
  external_population <- setdiff(data, rbind(TA, CC))
  EC_world <- external_population %>% filter(RANDASSIGN==0)
  
  #step: b
  #EC_candidate <- sample_n(EC_world, 1000) #choosing randomly, can instead choose to be different than CC
  EC_candidate <- EC_world %>%              #biased EC candidate
    group_by(Age_Group, Gender) %>% 
    sample_n(size=250, replace = TRUE) %>%
    ungroup()
  
  #step: c
  EC <- get_Matched_EC (TA, nrow(CC), EC_candidate)
  
  #Need a checker to check for missing categories
  flag <- check_subgroup_counts(age_sbgrps = 2, gender_sbgrps = 2, race_sbgrps = 5, 
                                TA = TA, EC=EC, CC = CC)
  if(flag==0){
    #end_iteration <- end_iteration + 1
    cat(paste0("Iteration ", i, " aborted due to missing categories. Continuting to next iteration. \n"))
    i <- i + 1
    next
  }
  
  #step: d
  HC <- rbind(CC, EC)
  
  #step: g
  merged_sample <- rbind(TA, HC)
  
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
  ATT_list <- c(ATT_list, ATT)
  
  #store LDM from each run for LDM - TA, CC, EC_Candidate, EC, HC, IPF_TA, IPF_HC
  LDM_summary <- get_LDM_summary(TA, CC, EC_candidate, EC, HC, EC_world, EC_world)
  
  LDM_summary <- LDM_summary %>% select(VAR:LDM_HC)
  
  #store LDM summary for each round as a list item
  LDM_summary_list <- c(LDM_summary_list, list(LDM_summary))
  
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

write.csv(as.data.frame(unlist(ATT_list)) %>% rename(PATT_m1 = 'unlist(ATT_list)'), 
          "./Data/Results/PATT_m1.csv")

ldm_sum <- as.data.frame(LDM_summary_list[1]) %>% select(LDM_TA:LDM_HC)
if (length(LDM_summary_list)>=2){
  for (ldm in LDM_summary_list[2:length(LDM_summary_list)]){
    ldm <- as.data.frame(ldm)
    ldm_dat <- ldm %>% select(LDM_TA:LDM_HC)
    ldm_sum <- ldm_sum + ldm_dat
  }
}

ldm_mean <- ldm_sum/length(LDM_summary_list)
ldm_mean <- cbind(LDM_summary%>%select(VAR, Level), ldm_mean)

#### Plot Colored Table #####
customGreen = "#d6e6e8" # Representative
customRedNeg = "#cf846e" #Highly Underrepresented
customRedPos = '#697ba3' #Highly Overrepresented
customYellowNeg = '#e6bcac' #Underrepresented
customYellowPos = '#a6b0cc' #Overrepresented

final_dat <- datatable(ldm_mean, options = list(scrollX = T)) %>% 
  formatStyle(c('LDM_TA', "LDM_CC", "LDM_EC_candidate", "LDM_EC", "LDM_HC"),
              backgroundColor = styleInterval(c(-0.51, -0.22, 0.22, 0.51), 
                                              c(customRedNeg, customYellowNeg, customGreen,
                                                customYellowPos, customRedPos))) %>%
  formatRound(c('LDM_TA', "LDM_CC", "LDM_EC_candidate", "LDM_EC", "LDM_HC"), digits = 3)

final_dat