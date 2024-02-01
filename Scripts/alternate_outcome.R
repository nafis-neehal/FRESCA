library(survival)
library(survminer)
library(dplyr)
library(gtsummary)

setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
nhanes <- read.csv("./Data/Processed/nhanes_hypertension.csv")
monthly_SPRINT_data <- read.csv("./Data/Processed/SPRINT_monthly_bp.csv")
baseline_SPRINT_data <- read.csv("./Data/Processed/SPRINT_example.csv")
#baseline_SPRINT_data <- na.omit(baseline_SPRINT_data)
primary_outcome_data <- read.csv("./Data/Processed/primarysurv.csv")
primary_outcome_data_postin <- primary_outcome_data %>% filter(EVENTDAYS_POSTI==1)

outcome_data <- baseline_SPRINT_data %>%
  inner_join(primary_outcome_data, by="MASKID") %>%
  select(MASKID:GFRESTIMATE, EVENTDAYS, EVENTDAYS_POSTI, EVENT)

# sfit <- survfit(Surv(EVENTDAYS, EVENT)~RANDASSIGN, 
#                data=outcome_data)
# ggsurvplot(sfit, pval=TRUE, risk.table = TRUE)

get_effect <- function(population, weights){
  if (!is.null(weights)){
    fit <- coxph(Surv(EVENTDAYS, EVENT)~RANDASSIGN, data=population, weights = weights)
  }
  else{
    fit <- coxph(Surv(EVENTDAYS, EVENT)~RANDASSIGN, data=population)
  }
  return (coef(fit))
}
