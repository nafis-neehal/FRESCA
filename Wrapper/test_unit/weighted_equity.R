setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
source("./Wrapper/load_packages.R")
source("./Wrapper/wrapper_modules.R")
source("./Wrapper/load_data_formatted.R")
source("./Modules/equity_metrics_miao.R")

#### STUDY PARAMETERS ####
TA_Size <- 2000
CC_Size <- 2000
SC_Size <- 1000
IPF_maxiter <- 100
num_seed <- 100
read_target <- F
var_name_list <- c('Gender', 'Age_Group', 'Race_or_Ethnicity')
var_subgroup_list <- list(c('Male','Female'),
                          c('40-59','59+'),
                          c('NH White', 'NH Black', 'NH Asian', 'Hispanic', 'Other'))

# var_name_list <- c('Gender', 'Age_Group')
# var_subgroup_list <- list(c('Female','Male'),  
#                           c('40-59','59+'))

# var_name_list <- c('Gender', 'Race_or_Ethnicity')
# var_subgroup_list <- list(c('Female','Male'),
#                           c('NH White', 'NH Black', 'NH Asian', 'Hispanic', 'Other'))

####Target Data####
target_data <- read.csv("./Data/Results/NHANES_Gen/NHANES_regen_200k.csv")
target_data$X <- NULL
target_data$ratio <- NULL
###################

all_data <- read_data(partition_seed=1, import_target=read_target)
target_data_count <- read.csv("./Data/Results/NHANES_Gen/NHANES_count_level3.csv")

if (read_target==T){
  RCT <- data.frame(all_data[[1]][[1]])
  BIASED_EC <- data.frame(all_data[[1]][[2]])
} else{
  RCT <- data.frame(all_data[[1]])
  BIASED_EC <- data.frame(all_data[[2]])
}

trial_TR_weights <- get_adjusted_weights(RCT %>% filter(RANDASSIGN==1), target_data,
                                      "lr", weight_method="gen")

trial_TR_weights_scaled <- scale_vector(unlist(trial_TR_weights),  100)

remove(all_data)

generate_subgroups <- function(V, n) {
  if (n == 1) {
    # Flatten all vectors into a single list
    result <- unlist(V)
  } else {
    # Generate combinations of elements from different vectors up to level n
    result <- combn(V, n, function(x) {
      apply(expand.grid(x), 1, paste, collapse = "/")
    }, simplify = FALSE)
  }

  # Remove duplicates
  result <- unique(unlist(result))

  return(result)
}

######### Test #########
# V <- list(c('a', 'b'), c('c', 'd'), c('e', 'f', 'g', 'h', 'i'))
# n <- 3
# output <- generate_subgroups(V, n)
# p <- as.list(strsplit(output[i], "")[[1]])
# p <- p[p != " "]
# print(output)
########################

subgroup_in_which_var <- function(var_name_list, var_subgroup_list, name){
  for (i in 1:length(var_subgroup_list)){
    list_members <- var_subgroup_list[i]
    if (unlist(name) %in% unlist(list_members)){
      return(var_name_list[i])
    }
  }
}


get_LD_nth_level <- function(trial_pop, target_pop, trial_weights, var_name_list, var_subgroup_list, n){
  
  #create LDM summary matrix
  df <- data.frame(matrix(ncol = 4, nrow = 0))
  x <- c("Subgroup", "Background_Rate", "Observed_Rate", "LD")
  colnames(df) <- x
  
  #combine trial data and weights into one df for slicing easiness
  trial_pop$weights <- trial_weights
  
  #generate var combinations on n-th level subgroup
  output <- generate_subgroups(var_subgroup_list, n)
  
  for (i in 1:length(output)){ #do action for each n-th level subgroup
    print(output[i])
    
    trial_copy_grp_member <- trial_pop
    target_copy_grp_member <- target_pop
    
    #take each of the subgroup combination
    #Female/NH Black -> ("Female", "NH Black")
    sub_group_name <- as.list(strsplit(output[i], "/")[[1]])
    for (j in 1:length(sub_group_name)){
      var_name <- subgroup_in_which_var(var_name_list, var_subgroup_list, sub_group_name[j])
      trial_copy_grp_member <- trial_copy_grp_member %>% filter(!!sym(var_name)==sub_group_name[j])
      target_copy_grp_member <- target_copy_grp_member %>% filter(!!sym(var_name)==sub_group_name[j])
    }
    
    if((n==1)|(n==length(var_name_list))){ #if we use first level or lowest-level subgroup
      if (is.null(trial_weights)){
        trial_grp_member_count <- nrow(trial_copy_grp_member)
        trial_total_count <- nrow(trial_pop)
      }
      else{
        trial_grp_member_count <- ceiling(sum(trial_copy_grp_member$weights))
        trial_total_count <- ceiling(sum(trial_weights))
      }
      #trial_total_count <- nrow(trial_pop) 
      target_grp_member_count <- sum(target_copy_grp_member$n)
      target_total_count <- sum(target_pop$n)
      
      bg_rate <- Rate_Calculation(target_grp_member_count, target_total_count)
      user_rate <- Rate_Calculation(trial_grp_member_count, trial_total_count)
      LDM <- Log_Disparate_Impact(as.numeric(round(bg_rate, 4)), as.numeric(round(user_rate,4)), for_plot = FALSE)
      df[nrow(df) + 1,] <- c(output[i], as.numeric(round(bg_rate, 4)), as.numeric(round(user_rate,4)), 
                              as.numeric(round(abs(LDM),4)))
    }
  }
  return (df)
}

##############
#  first level subgroup counts sums up to 2500
#  second/intermediate level subgroup counts -> how to deal, later question!!!
#  hint: try thinking about splitting the whole population into the subgroups and variables.
#  example: if Female:Black, split into F:B, F:W, M:B, M:L that's the total population
#  lowest level subgroup counts sums up to 2500
###############

#ISSUE: generalize package gives weights for each lowest level subgroup
#Solution: instead of counting nrow(RCT), sum the available weights

#new_RCT <- sample_n(RCT, nrow(RCT), replace = T, weight = trial_weights)
ldm_summary <- get_LD_nth_level(RCT %>% filter(RANDASSIGN==1), target_data_count, trial_TR_weights, 
                 var_name_list, var_subgroup_list, 1) #try 1 and 3 for now

# table(target_data$Race_or_Ethnicity)/nrow(target_data) #created target data
# val <- target_data_count %>% filter(Race_or_Ethnicity=="NH Asian") %>% summarise(sum(n))
# val / sum(target_data_count$n)

