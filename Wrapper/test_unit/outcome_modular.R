#### SOURCES ####
setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")

allhat_data <- read.csv("./Data/ALLHAT_data/ALLHAT_processed_secondary.csv")
allhat_data$X <- NULL

#create the treated populations
ctrl <- allhat_data %>% filter(RANDASSIGN==2) #Ch
tr1 <- allhat_data %>% filter(RANDASSIGN==3) #Am
tr2 <- allhat_data %>% filter(RANDASSIGN==4) #Li


## Split the Trial population into TA, CC and EC
## Change RANDASSIGN to 0 for controls and 1 for treated
g1 <- rbind(tr1, ctrl)
g1 <- g1 %>% mutate(RANDASSIGN = ifelse(RANDASSIGN==2, 0, 1))
g2 <- rbind(tr2, ctrl)
g2 <- g2 %>% mutate(RANDASSIGN = ifelse(RANDASSIGN==2, 0, 1))

get_cox_effect <- function(population, outcome_days, outcome_event, weights){
  
  # Extract survival time and status vectors from the dataframe
  survival_time <- population[[outcome_days]]
  survival_status <- population[[outcome_event]]
  
  # Create a survival object
  survival_data <- Surv(time = survival_time, event = survival_status)
  
  if (!is.null(weights)){
    fit <- coxph(survival_data~RANDASSIGN, data=population, weights = weights)
  }
  else{
    fit <- coxph(survival_data~RANDASSIGN, data=population)
  }
  
  estimate <- coef(fit)['RANDASSIGN']
  return (estimate)
}

exp(get_cox_effect(g1, "EVENTDAYS", "EVENT", NULL))



