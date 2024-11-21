setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
source("./Wrapper/load_packages.R")
source("./Wrapper/wrapper_modules.R")
source("./Wrapper/load_data_formatted.R")
source("./Modules/equity_metrics_miao.R")

#### STUDY PARAMETERS ####
TA_Size <- 2000
CC_Size <- 1000
SC_Size <- 1000
IPF_maxiter <- 100
scale_to <- 100
num_seed <- 50
read_target <- F
weight_method <- "ipf" #"gen" or "ipf"

var_name_list <- c('Gender', 'Age_Group', 'Race_or_Ethnicity')
var_subgroup_list <- list(c('Male','Female'),
                          c('40-59','59+'),
                          c('NH White', 'NH Black', 'NH Asian', 'Hispanic', 'Other'))

# var_name_list <- c('Gender','Race_or_Ethnicity')
# var_subgroup_list <- list(c('Male','Female'),
#                           c('NH White', 'NH Black', 'NH Asian', 'Hispanic', 'Other'))

# var_name_list <- c('Gender', 'Age_Group')
# var_subgroup_list <- list(c('Male','Female'),
#                           c('40-59','59+'))

#Whole Sprint data
baseline_SPRINT_data <- read.csv("./Data/Processed/SPRINT_example.csv")
baseline_SPRINT_data <- na.omit(baseline_SPRINT_data)
primary_outcome_data <- read.csv("./Data/Processed/primarysurv.csv")
data <- baseline_SPRINT_data %>%
  inner_join(primary_outcome_data, by="MASKID") %>%
  mutate(CVDHISTORY  = factor(CVDHISTORY),
         CVDPOINTS   = if_else(CVDPOINTS<=19,'Moderate', 'High'),
         GFRESTIMATE = if_else(GFRESTIMATE>=60, 'Normal', 'Disease')) %>%
  select(MASKID:GFRESTIMATE, EVENTDAYS, EVENTDAYS_POSTI, EVENT)

####Target Data####
target_data <- read.csv("./Data/Results/NHANES_Gen/NHANES_regen_20k.csv")
target_data$X <- NULL
target_data$ratio <- NULL
###################

run_main_simulation <- function(num_seed){
  
  TR_summary <- data.frame(matrix(ncol = 5, nrow = 0))
  CN_summary <- data.frame(matrix(ncol = 5, nrow = 0))
  x <- c("Subgroup", "Background_Rate", "Observed_Rate", "LD", "Seed")
  colnames(TR_summary) <- x
  colnames(CN_summary) <- x
  
  for (i in 1:num_seed){
    
    cat("Iteration", i, "\n")
    set.seed(i)
    #hybrid_RCT <- sample_n(data, size=nrow(data), replace = T)
    
    #Step 1: create data
    all_data <- read_data(partition_seed=i, import_target=read_target)
    target_data_count <- read.csv("./Data/Results/NHANES_Gen/NHANES_count_level3.csv")

    if (read_target==T){
      RCT <- data.frame(all_data[[1]][[1]])
      BIASED_EC <- data.frame(all_data[[1]][[2]])
    } else{
      RCT <- data.frame(all_data[[1]])
      BIASED_EC <- data.frame(all_data[[2]])
    }

    remove(all_data)

    #Step 2: Borrow Synthetic Controls
    synthetic_controls <- get_matched_controls(method_name = "one_to_one_with_replace_caliper",
                                               RCT = RCT, BIASED_EC = BIASED_EC,
                                               N = SC_Size)

    hybrid_RCT <- rbind(RCT, synthetic_controls)

    cat("Synthetic Control borrowing done...\n")
    
    #Step 3: Measure Equity in Hybrid Treatment Arm
    trial_TR_weights_scaled <- get_adjusted_weights(hybrid_RCT %>% filter(RANDASSIGN==1), target_data,
                                             proba_estimate_method="lr", weight_method=weight_method)

    if(sum(unlist(trial_TR_weights_scaled))!=scale_to){
      trial_TR_weights_scaled <- scale_vector(unlist(trial_TR_weights_scaled), scale_to)
    }

    TR_ldm_summary <- get_LD_nth_level(hybrid_RCT %>% filter(RANDASSIGN==1), target_data_count, 
                                       trial_TR_weights_scaled,
                                       var_name_list, var_subgroup_list, 1) #try 1 and 3 for now

    TR_ldm_summary$Seed <- i

    TR_summary <- rbind(TR_summary, TR_ldm_summary)
    
    cat("Treatment group equity measured...\n")
    
    #Step 4: Measure Equity in Hybrid Control Arm
    trial_CN_weights_scaled <- get_adjusted_weights(hybrid_RCT %>% filter(RANDASSIGN==0), target_data,
                                             proba_estimate_method="lr", weight_method=weight_method)

    if(sum(unlist(trial_CN_weights_scaled))!=scale_to){
      trial_CN_weights_scaled <- scale_vector(unlist(trial_CN_weights_scaled), scale_to)
    }

    CN_ldm_summary <- get_LD_nth_level(hybrid_RCT %>% filter(RANDASSIGN==0), target_data_count, 
                                       trial_CN_weights_scaled,
                                       var_name_list, var_subgroup_list, 1) #try 1 and 3 for now

    CN_ldm_summary$Seed <- i

    CN_summary <- rbind(CN_summary, CN_ldm_summary)
    
    cat("Control group equity measured...\n")
    
  }
  return(list(TR_summary, CN_summary))
}

start <- Sys.time()
summary <- run_main_simulation(num_seed = num_seed)
end <- Sys.time()
paste("Total time elapsed:", end-start)

tr_summary <- as.data.frame(summary[1])
cn_summary <- as.data.frame(summary[2])

tr_summary_df <- tr_summary %>% group_by(Subgroup) %>% summarise(mean_LD = mean(as.numeric(LD)),
                                                lower_ci_LD = CI(as.numeric(LD))["lower"],
                                                upper_ci_LD = CI(as.numeric(LD))["upper"])

cn_summary_df <- cn_summary %>% group_by(Subgroup) %>% summarise(mean_LD = mean(as.numeric(LD)),
                                                                 lower_ci_LD = CI(as.numeric(LD))["lower"],
                                                                 upper_ci_LD = CI(as.numeric(LD))["upper"])


