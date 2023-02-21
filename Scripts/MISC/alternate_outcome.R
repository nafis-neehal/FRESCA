nhanes <- read.csv("./Data/NHANES/NHANES_age_adj/nhanes_hypertension.csv")
monthly_SPRINT_data <- read.csv("./Data/Processed/SPRINT_monthly_bp.csv")
baseline_SPRINT_data <- read.csv("./Data/Processed/SPRINT_example.csv")
#baseline_SPRINT_data <- na.omit(baseline_SPRINT_data)
primary_outcome_data <- read.csv("./Data/SPRINT_data/primarysurv.csv")
bp_base_data <- read.csv("./Data/SPRINT_data/bp_manage_base.csv")

# month <- '18M'
# month_df <- monthly_SPRINT_data %>% 
#   filter(VISITCODE==month) %>%
#   select(MASKID, SEATSYS, FORMDAYS)
# data <- baseline_SPRINT_data %>% 
#   inner_join(month_df, by='MASKID') %>% 
#   select(MASKID:SEATSYS, FORMDAYS)
# data <- na.omit(data)
# 
# outcome_data <- data %>%
#   inner_join(primary_outcome_data, by="MASKID") %>%
#   select(MASKID:SEATSYS, FORMDAYS, EVENTDAYS, EVENTDAYS_POSTI, EVENT)
# 
# outcome_data_hasevent <- outcome_data %>%
#   filter(FORMDAYS < EVENTDAYS)

bp <- bp_base_data %>% inner_join(baseline_SPRINT_data, by="MASKID")
bp_filter <- bp %>% select(MASKID, SEATSYS, RANDASSIGN)

outcome_data <- bp_filter %>%
  inner_join(primary_outcome_data, by="MASKID") %>%
  select(SEATSYS:EVENTDAYS, EVENT)

library(survival)
library(survminer)

sfit <- survfit(Surv(EVENTDAYS, EVENT)~RANDASSIGN, 
               data=outcome_data)
summary(sfit)
ggsurvplot(sfit, pval=TRUE, risk.table = TRUE)

cfit <- coxph(Surv(EVENTDAYS, EVENT)~RANDASSIGN, data=outcome_data)
cfit$coefficient

