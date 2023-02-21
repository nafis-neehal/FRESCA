setwd("/home/neehan/data/Nafis/ESCA")
source("./Modules/p_val_add.R")
source("./Modules/equity_metrics_miao.R")
source("./Scripts/common.R")


#Save ATT for next bootstraps
ATT_list <- list()
LDM_summary_list <- list()

total_success <- 20
success_counter <- 1
i <- 1

while(success_counter <= total_success){
  
  #step: a
  TA <- data %>% filter(RANDASSIGN==1)
  CC <- data %>% filter(RANDASSIGN==0)
  
  HC <- CC
  
  #Need a checker to check for missing categories
  flag <- check_subgroup_counts(age_sbgrps = 2, gender_sbgrps = 2, race_sbgrps = 5,
                                gfr_sbgrps = 2,
                                TA = TA, HC=HC)
  if(flag==0){
    #end_iteration <- end_iteration + 1
    cat(paste0("Iteration ", i, " aborted due to missing categories. Continuting to next iteration. \n"))
    i <- i + 1
    next
  }
  
  #step: e
  W_TA_ipf <- get_IPF_weights(dat = TA, maxIter = 50)
  W_HC_ipf <- get_IPF_weights(dat = HC, maxIter = 50)
  
  #step: f
  set.seed(i)
  TA_IPF_sample <- TA[sample(seq_len(nrow(TA)), size = nrow(TA), replace = TRUE, prob = W_TA_ipf),]
  
  set.seed(i)
  HC_IPF_sample <- HC[sample(seq_len(nrow(HC)), size = nrow(HC), replace = TRUE, prob = W_HC_ipf),]
  
  #step: g
  merged_sample <- rbind(TA_IPF_sample, HC_IPF_sample)
  
  merged_sample <- merged_sample %>% mutate(Age_Group = factor(Age_Group),
                             Gender = factor(Gender),
                             Race_or_Ethnicity = factor(Race_or_Ethnicity),
                             Education = factor(Education),
                             Smoker = factor(Smoker),
                             FPG = factor(FPG),
                             TC = factor(TC),
                             GFRESTIMATE = factor(GFRESTIMATE))
  
  m.out <- matchit(RANDASSIGN ~ Age_Group + Gender + Race_or_Ethnicity + Education +
                     Smoker + FPG + TC + CVDHISTORY + CVDPOINTS + SERUMCREAT +
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
  
  # #store LDM from each run for LDM - TA, CC, EC_Candidate, EC, HC, IPF_TA, IPF_HC
  # LDM_summary <- get_LDM_summary(TA, CC, CC, CC, HC, TA_IPF_sample, HC_IPF_sample)
  # LDM_summary <- LDM_summary %>% select(VAR:LDM_TA, LDM_HC, LDM_IPF_TA, LDM_IPF_HC)
  # 
  # #store LDM summary for each round as a list item
  # LDM_summary_list <- c(LDM_summary_list, list(LDM_summary))
  
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
          "./Data/Results/PATT_ground20.csv")

gtdat <- read.csv("./Data/Results/PATT_ground20.csv")

# ldm_sum <- as.data.frame(LDM_summary_list[1]) %>% select(LDM_TA:LDM_IPF_HC)
# if (length(LDM_summary_list)>=2){
#   for (ldm in LDM_summary_list[2:length(LDM_summary_list)]){
#     ldm <- as.data.frame(ldm)
#     ldm_dat <- ldm %>% select(LDM_TA:LDM_IPF_HC)
#     ldm_sum <- ldm_sum + ldm_dat
#   }
# }

# ldm_mean <- ldm_sum/length(LDM_summary_list)
# ldm_mean <- cbind(LDM_summary%>%select(VAR, Level), ldm_mean)

#### Plot Colored Table #####
# customGreen = "#d6e6e8" # Representative
# customRedNeg = "#cf846e" #Highly Underrepresented
# customRedPos = '#697ba3' #Highly Overrepresented
# customYellowNeg = '#e6bcac' #Underrepresented
# customYellowPos = '#a6b0cc' #Overrepresented
# 
# final_dat <- datatable(ldm_mean, options = list(scrollX = T)) %>% 
#   formatStyle(c('LDM_TA', "LDM_HC", "LDM_IPF_TA", "LDM_IPF_HC"),
#               backgroundColor = styleInterval(c(-0.51, -0.22, 0.22, 0.51), 
#                                               c(customRedNeg, customYellowNeg, customGreen,
#                                                 customYellowPos, customRedPos))) %>%
#   formatRound(c('LDM_TA', "LDM_HC", "LDM_IPF_TA", "LDM_IPF_HC"), digits = 3)
# 
# final_dat
