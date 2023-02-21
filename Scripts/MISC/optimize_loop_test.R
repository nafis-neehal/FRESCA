setwd("/Users/nafisneehal/ESCA")
source("./Scripts/common.R")
source("./Modules/p_val_add.R")
source("./Modules/equity_metrics_miao.R")

library(doSNOW)
library(foreach)

cores=detectCores(logical = FALSE)
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoSNOW(cl)

iterations <- 20
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

seed_value <- 1
total_seeds <- 1

get_ipf <- function(seed){
  set.seed(seed_value)
  TA <- sample_n(data %>% filter(RANDASSIGN==1), 1000)
  W_TA_ipf <- get_IPF_weights(dat = TA, maxIter = 10)
  return(W_TA_ipf)
}

w_par <- foreach(seed_value=seq(0,1000,50), .combine = 'cbind', .packages = c("dplyr", "ipfr"),
                 .options.snow = opts) %dopar% {
  get_ipf(seed_value)
}
close(pb)
stopCluster(cl)
