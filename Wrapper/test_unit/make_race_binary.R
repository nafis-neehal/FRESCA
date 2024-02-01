setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
source("./Wrapper/load_packages.R")
source("./Wrapper/wrapper_modules.R")
source("./Wrapper/load_data_formatted.R")
source("./Modules/equity_metrics_miao.R")

TA_Size <- 2000
CC_Size <- 1000
SC_Size <- 1000
IPF_maxiter <- 100
scale_to <- 100
num_seed <- 2
read_target <- F
weight_method <- "gen" #"gen" or "ipf"

create_dummies <- function(df, column){
  df_with_dummy <- dummy_cols(df, select_columns = column, 
                              remove_selected_columns = T)
  for ( col in 1:ncol(df_with_dummy)){
    colnames(df_with_dummy)[col] <-  sub(paste(column,"_", sep = ""), "", 
                                         colnames(df_with_dummy)[col]) #removes main var name from column
    colnames(df_with_dummy)[col] <-  sub(" ", "_",                     #removes space with _ in each column
                                         colnames(df_with_dummy)[col])
  }
  
  return(df_with_dummy)
}

recode_race <- function(df){
  df$NH_Asian <- ifelse(df$NH_Asian==1,"NH_Asian_1", "NH_Asian_0")
  df$NH_White <- ifelse(df$NH_White==1,"NH_White_1", "NH_White_0")
  df$NH_Black <- ifelse(df$NH_Black==1,"NH_Black_1", "NH_Black_0")
  df$Hispanic <- ifelse(df$Hispanic==1,"Hispanic_1", "Hispanic_0")
  df$Other    <- ifelse(df$Other==1,"Other_1", "Other_0")
  return (df)
}


all_data <- read_data(partition_seed=1, import_target=read_target)
target_data_count <- read.csv("./Data/Results/NHANES_Gen/NHANES_count_level3.csv")

####Target Data####
target_data <- read.csv("./Data/Results/NHANES_Gen/NHANES_regen_20k.csv")
target_data$X <- NULL
target_data$ratio <- NULL
###################
RCT <- data.frame(all_data[[1]])
BIASED_EC <- data.frame(all_data[[2]])

### Create dummies ###
target_data_dummies <- create_dummies(target_data, column = "Race_or_Ethnicity")
target_data_count_dummy <- create_dummies(target_data_count, column = "Race_or_Ethnicity")
target_data_count_dummy <- target_data_count_dummy[, c(2,3,5,6,7,8,9,4)]
RCT_dummy <- create_dummies(RCT, column = "Race_or_Ethnicity")
BIASED_EC_dummy <- create_dummies(BIASED_EC, column = "Race_or_Ethnicity")
######################

### Remove unnecessary dataframes 
remove(RCT, BIASED_EC, target_data, target_data_count, all_data)

### Create adjusted weights + scale
RCT_dummy_weights <- get_adjusted_weights(RCT_dummy, target_data_dummies,
                                               proba_estimate_method="lr", weight_method=weight_method)
if(sum(unlist(RCT_dummy_weights))!=scale_to){
  RCT_dummy_weights <- scale_vector(unlist(RCT_dummy_weights), scale_to)
}

####### Make sure that each subgroup of a feature is uniquely identifiable #######
# This is needed for LDM Calculation only, not needed for weight generation
# Find a workaround later if needed
# RCT_dummy$NH_Asian <- ifelse(RCT_dummy$NH_Asian==1,"NH_Asian_1", "NH_Asian_0")
# RCT_dummy$NH_White <- ifelse(RCT_dummy$NH_White==1,"NH_White_1", "NH_White_0")
# RCT_dummy$NH_Black <- ifelse(RCT_dummy$NH_Black==1,"NH_Black_1", "NH_Black_0")
# RCT_dummy$Hispanic <- ifelse(RCT_dummy$Hispanic==1,"Hispanic_1", "Hispanic_0")
# RCT_dummy$Other    <- ifelse(RCT_dummy$Other==1,"Other_1", "Other_0")

RCT_dummy <- recode_race(RCT_dummy)
target_data_count_dummy <- recode_race(target_data_count_dummy)
                         
# target_data_count_dummy$NH_Asian <- ifelse(target_data_count_dummy$NH_Asian==1,"NH_Asian_1", "NH_Asian_0")
# target_data_count_dummy$NH_White <- ifelse(target_data_count_dummy$NH_White==1,"NH_White_1", "NH_White_0")
# target_data_count_dummy$NH_Black <- ifelse(target_data_count_dummy$NH_Black==1,"NH_Black_1", "NH_Black_0")
# target_data_count_dummy$Hispanic <- ifelse(target_data_count_dummy$Hispanic==1,"Hispanic_1", "Hispanic_0")
# target_data_count_dummy$Other    <- ifelse(target_data_count_dummy$Other==1,"Other_1", "Other_0")


var_name_list <- c('Gender', 'Age_Group', 'NH_White', 'NH_Black', 'NH_Asian', 'Hispanic', 'Other')
var_subgroup_list <- list(c('Male','Female'),
                          c('40-59','59+'),
                          c('NH_White_0','NH_White_1'),
                          c('NH_Black_0','NH_Black_1'),
                          c('NH_Asian_0','NH_Asian_1'),
                          c('Hispanic_0','Hispanic_1'),
                          c('Other_0','Other_1'))

ldm_summary <- get_LD_nth_level( RCT_dummy, target_data_count_dummy, RCT_dummy_weights,
                                   var_name_list, var_subgroup_list, 1) #try 1 for now, 3 has some issues


