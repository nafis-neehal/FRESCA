setwd("/home/neehan/data/Nafis/ESCA")
source("./Scripts/common.R")
source("./Modules/p_val_add.R")
source("./Modules/equity_metrics_miao.R")

get_ATT_val <- function(population, weights){
  
  #run linear model to get the treatment effect
  if (!is.null(weights)){
    fit <- lm(SEATSYS ~ RANDASSIGN, data = population, weights = weights)
  }
  else{
    fit <- lm(SEATSYS ~ RANDASSIGN, data = population)
  }
  
  estimate <- coef(fit)['RANDASSIGN']
  
  return(estimate)
}

TA_World <- data %>% filter(RANDASSIGN==1)
TA_size <- 1000
total_seeds <- 20
IPF_maxiter <- 50

seed_value <- 1

set.seed(seed_value)
TA <- sample_n(TA_World %>% filter(RANDASSIGN==1), TA_size)
W_TA_ipf <- get_IPF_weights(dat = TA, maxIter = IPF_maxiter)
TA$W_IPF <- W_TA_ipf

set.seed(seed_value)
CC_World <- sample_n(data %>% filter(RANDASSIGN==0), TA_size)

external_population <- setdiff(data, rbind(TA_World, CC_World))

EC_world <- external_population %>% filter(RANDASSIGN==0)
W_EC_ipf <- get_ECWorld_Bias_weights(dat = EC_world, maxIter = IPF_maxiter)
EC_world$W_IPF <- W_EC_ipf

#TA and EC with weights, without IPF resamples
get_ATT_val(population=rbind(TA %>% select(MASKID:SEATSYS), EC_world %>% select(MASKID:SEATSYS)), weights = c(TA$W_IPF, EC_world$W_IPF))

#TA and EC without weights, without IPF resamples
get_ATT_val(population=rbind(TA %>% select(MASKID:SEATSYS), EC_world %>% select(MASKID:SEATSYS)), weights = NULL)

#TA and EC without weights, with IPF resamples
set.seed(seed_value)
TA_IPF_sample <- TA[sample(seq_len(nrow(TA)), size = nrow(TA), replace = TRUE, prob = TA$W_IPF),]
set.seed(seed_value)
EC_IPF_sample <- EC_world[sample(seq_len(nrow(EC_world)), size = nrow(TA), replace = TRUE, prob = EC_world$W_IPF),]

get_ATT_val(population=rbind(TA_IPF_sample %>% select(MASKID:SEATSYS), EC_IPF_sample %>% select(MASKID:SEATSYS)), weights = NULL)

#TA and EC with weights, with IPF resamples
get_ATT_val(population=rbind(TA_IPF_sample %>% select(MASKID:SEATSYS), EC_IPF_sample %>% select(MASKID:SEATSYS)), 
            weights = c(TA_IPF_sample$W_IPF, EC_IPF_sample$W_IPF))





