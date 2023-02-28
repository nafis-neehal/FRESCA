suppressPackageStartupMessages({
  library(dplyr) #
  library(table1) #
  library(MatchIt) #
  library(lmtest) 
  library(sandwich)
  library(ggplot2) #
  library(cowplot)
  library(data.table)#
  library(tidyr) #
  library(tibble) #
  library(ipfr) # 
  library(rlang) #
  library(gridExtra) #
  library(cobalt) #
  library(patchwork) #
  library(marginaleffects) #
  library(doSNOW) #
  library(doParallel) #
  library(foreach) #
  library(ggpubr)
  library(gtools) #gives character based sort in mixedsort() #
  library(tableone) #for calculating ASMD
  library(reshape2) #for making long form data with datamelt
  library(survival)
  library(survminer)
  library(gtsummary)
})


get_IPF_weights <- function(dat, maxIter){
  #fix the target distributions from TP
  targets <- list()
  targets$Age_Group <- tibble(
    '0' = 31.184,
    '1' = 68.816
  )
  targets$Gender <- tibble(
    'Female' = 55.363,
    'Male' = 44.637
  )
  targets$Race_or_Ethnicity <- tibble(
    '0' = 69.266,
    '1' = 12.023,
    '2' = 3.912,
    '3' = 10.035,
    '4' = 4.763
  )

  var_list <- c("Age_Group", "Gender", "Race_or_Ethnicity")

  dat <- dat %>% select(unlist(var_list))
  
  dat2 <- dat %>% mutate(Age_Group = recode(Age_Group,
                                             '40-59' = 0,
                                             '59+' = 1),
                          Race_or_Ethnicity = recode(Race_or_Ethnicity,
                                                     'NH White' = 0,
                                                     'NH Black' = 1,
                                                     'NH Asian' = 2,
                                                     'Hispanic' = 3,
                                                     'Other' = 4))
  
  dat2 <- as_tibble(dat2)
  result <- ipu(dat2, targets, max_iterations = maxIter)
  weights <- result$weight_tbl$weight
  
  return(weights)
  
}

#note: using military population ratio
get_ECWorld_Bias_weights <- function(dat, maxIter){
  #fix the target distributions from TP
  targets <- list()
  # targets$Age_Group <- tibble(
  #   '1' = 40, #59+ =   40   #NHANES=69
  #   '0' = 60  #40-59 = 60   #NHANES=31
  # )
  targets$Race_or_Ethnicity <- tibble(
    '1' = 20.2, #black original = 28, change = 20.2
    '0' = 54, #white original = 57, change = 54
    '3' = 17.2,  #hispanic original = 12, change = 17.2
    '4' = 1.7, #other orig = 0.3, change = remaining
    '2' = 6.9   #asian orig = 1, change = 6.9
  )
  targets$Gender <- tibble(
    'Male' = 72.6, #change=72.6
    'Female' = 27.4   #change=27.4
  )
  # targets$GFRESTIMATE <- tibble(
  #   'Normal' = 60.0,
  #   'Disease' = 40.0
  # )
  targets$CVDPOINTS <- tibble(
    'High' = 45.0,
    'Moderate' = 55.0
  )
  targets$CVDHISTORY <- tibble(
    '1' = 45.000,
    '0' = 55.000
  )
  
  var_list <- c("Race_or_Ethnicity", 'Gender', 'CVDPOINTS', 'CVDHISTORY')
  
  dat <- dat %>% select(unlist(var_list))
  
  
  #Age_Group = recode(Age_Group, '40-59'=0, '59+'=1),
  dat2 <- dat %>% mutate(Race_or_Ethnicity = recode(Race_or_Ethnicity,
                                                    'NH White' = 0,
                                                    'NH Black' = 1,
                                                    'NH Asian' = 2,
                                                    'Hispanic' = 3,
                                                    'Other' = 4))
  
  dat2 <- as_tibble(dat2)
  result <- ipu(dat2, targets, max_iterations = maxIter)
  weights <- result$weight_tbl$weight
  
  return(weights)
  
}

get_LDM_from_pop <- function(trialPopulation, targetPopulation){
  ev_list <- c("Gender", "Age_Group", "Race_or_Ethnicity")
  
  ev_level_list <- list(list('Female','Male'), list('40-59','59+'), 
                        c("NH White", "NH Black", "NH Asian", "Hispanic", "Other"))
  
  ev_counter <- 1
  df_colnames <- c("VAR", "Level", "Trial_Grp", "Trial_Non_Grp", 
                   "Target_Grp", "Target_Non_Grp", "LDM")
  ldm_df <- data.frame(matrix(ncol=7,nrow=0, dimnames=list(NULL, df_colnames)))
  for (ev in ev_list){
    ev_levels <- unlist(ev_level_list[ev_counter])
    for (i in 1:length(ev_levels)){
      level <- ev_levels[i]
      trial_grp_mmbr <- trialPopulation %>% filter((!!sym(ev))==level) %>% count()
      trial_non_grp_mmbr <- trialPopulation %>% filter((!!sym(ev))!=level) %>% count()
      target_grp_mmbr <- targetPopulation %>% filter((!!sym(ev))==level) %>% count()
      target_non_grp_mmbr <- targetPopulation %>% filter((!!sym(ev))!=level) %>% count()
      
      bg_rate <- Rate_Calculation(target_grp_mmbr, target_grp_mmbr + target_non_grp_mmbr)
      user_rate <- Rate_Calculation(trial_grp_mmbr, trial_grp_mmbr + trial_non_grp_mmbr)
      LDM <- Log_Disparate_Impact(as.numeric(bg_rate),as.numeric(user_rate), for_plot = FALSE)
      ldm_df[nrow(ldm_df)+1, ] <- c(ev, level, trial_grp_mmbr, trial_non_grp_mmbr,
                                    target_grp_mmbr, target_non_grp_mmbr, LDM)
    }
    ev_counter <- ev_counter + 1
  }
  return(ldm_df)
}

get_LDM_from_count <- function(trialPopulation) {
  
  ev_list <- c("Gender", "Age_Group", "Race_or_Ethnicity")
  
  ev_level_list <- list(list('Female','Male'), list('40-59','59+'), 
                        c("NH White", "NH Black", "NH Asian", "Hispanic", "Other"))
  ev_counts_list <- list(list(30505435, 24594942), #F,M
                         list(17182561, 37917816), 
                         list(38165591, 6624866, 2155632, 5529827, 2624461)) #W,B,A,H,O
  
  ev_counter <- 1
  df_colnames <- c("VAR", "Level", "Trial_Grp", "Trial_Non_Grp", 
                   "Target_Grp", "Target_Non_Grp", "LDM")
  ldm_df <- data.frame(matrix(ncol=7,nrow=0, dimnames=list(NULL, df_colnames)))
  for (ev in ev_list){
    ev_levels <- unlist(ev_level_list[ev_counter])
    ev_counts <- unlist(ev_counts_list[ev_counter])
    for (i in 1:length(ev_levels)){ #level in unlist(ev_levels)
      level <- ev_levels[i]
      trial_grp_mmbr <- trialPopulation %>% filter((!!sym(ev))==level) %>% count()
      trial_non_grp_mmbr <- trialPopulation %>% filter((!!sym(ev))!=level) %>% count()
      target_grp_mmbr <- ev_counts[i]
      target_non_grp_mmbr <- sum(ev_counts) - target_grp_mmbr

      bg_rate <- Rate_Calculation(target_grp_mmbr, target_grp_mmbr + target_non_grp_mmbr)
      user_rate <- Rate_Calculation(trial_grp_mmbr, trial_grp_mmbr + trial_non_grp_mmbr)
      LDM <- Log_Disparate_Impact(as.numeric(bg_rate),as.numeric(user_rate), for_plot = FALSE)

      ldm_df[nrow(ldm_df)+1, ] <- c(ev, level, trial_grp_mmbr, trial_non_grp_mmbr,
                                    target_grp_mmbr, target_non_grp_mmbr, LDM)
    }
    ev_counter <- ev_counter + 1
  }
  return (ldm_df)
}

check_subgroup_counts <- function(gender_sbgrps, race_sbgrps, cvd_sbgrps, frs_sbgrps, TA, HC) {
  flag <- 1
  if ((length(unique(TA$CVDHISTORY)) != cvd_sbgrps) |
      (length(unique(TA$CVDPOINTS)) != frs_sbgrps) |
      (length(unique(TA$Gender)) != gender_sbgrps) |
      (length(unique(TA$Race_or_Ethnicity)) != race_sbgrps)){
    flag <- 0
    return(flag)
  }
  else if ((length(unique(HC$CVDHISTORY)) != cvd_sbgrps) |
           (length(unique(HC$CVDPOINTS)) != frs_sbgrps) |
           (length(unique(HC$Gender)) != gender_sbgrps) |
           (length(unique(HC$Race_or_Ethnicity)) != race_sbgrps)){
    flag <- 0
    return(flag)
  }
  return(flag)
}

get_LDM_summary<- function(TA, CC, EC_candidate, HC, IPF_TA, IPF_HC){
  LDM_TA <- get_LDM_from_count(TA) %>% select(VAR, Level, LDM) %>% rename(LDM_TA = LDM)
  LDM_CC <- get_LDM_from_count(CC) %>% select(LDM) %>% rename(LDM_CC = LDM)
  LDM_EC_Sample <- get_LDM_from_count(EC_candidate) %>% select(LDM) %>% rename(LDM_EC_candidate = LDM)
  #LDM_EC <- get_LDM_from_count(EC) %>% select(LDM) %>% rename(LDM_EC = LDM)
  LDM_HC <- get_LDM_from_count(HC) %>% select(LDM) %>% rename(LDM_HC = LDM)
  LDM_IPF_TA <- get_LDM_from_count(IPF_TA) %>% select(LDM) %>% rename(LDM_IPF_TA = LDM)
  LDM_IPF_HC <- get_LDM_from_count(IPF_HC) %>% select(LDM) %>% rename(LDM_IPF_HC = LDM)
  
  LDM_summary <- do.call("cbind", list(LDM_TA, LDM_CC, LDM_EC_Sample,  
                                       LDM_HC, LDM_IPF_TA, LDM_IPF_HC))
  
  LDM_summary <- as.data.frame(LDM_summary)
  
  return(LDM_summary)
}

nhanes <- read.csv("./Data/Processed/nhanes_hypertension.csv")
#monthly_SPRINT_data <- read.csv("./Data/Processed/SPRINT_monthly_bp.csv")
baseline_SPRINT_data <- read.csv("./Data/Processed/SPRINT_example.csv")
baseline_SPRINT_data <- na.omit(baseline_SPRINT_data)
primary_outcome_data <- read.csv("./Data/Processed/primarysurv.csv")

data <- baseline_SPRINT_data %>%
  inner_join(primary_outcome_data, by="MASKID") %>%
  mutate(CVDHISTORY  = factor(CVDHISTORY),
         CVDPOINTS   = if_else(CVDPOINTS<=19,'Moderate', 'High'),
         GFRESTIMATE = if_else(GFRESTIMATE>=60, 'Normal', 'Disease')) %>%
  select(MASKID:GFRESTIMATE, EVENTDAYS, EVENTDAYS_POSTI, EVENT)

datap <- data %>% mutate(EVENT = ifelse(EVENTDAYS<=1190, EVENT, 0))

#https://www.niddk.nih.gov/health-information/professionals/advanced-search/explain-kidney-test-results

# month <- '18M'
# month_df <- monthly_SPRINT_data %>% 
#   filter(VISITCODE==month) %>%
#   select(MASKID, SEATSYS)
# data <- baseline_SPRINT_data %>% 
#   inner_join(month_df, by='MASKID') %>% 
#   select(MASKID:SEATSYS)
# data <- na.omit(data)
# data <- data %>% mutate(CVDPOINTS = if_else((Gender=='Male' & CVDPOINTS>=13)|(Gender=='Female' & CVDPOINTS>=16), 'High', 'Low')
#                         ,GFRESTIMATE = ifelse(GFRESTIMATE>=60,'Normal','Disease')) %>%
#   select(MASKID:SEATSYS)

