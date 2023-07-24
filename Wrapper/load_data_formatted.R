

########## HELPER FUNCTIONS ###########
#define RCT TA and CC
get_rct_biased_ec_data <- function(seed_value, data, TA_Size, CC_Size, IPF_maxiter){
  
  TA_World <- data %>% filter(RANDASSIGN==1)
  
  #### RCT Population ####
  set.seed(seed_value)
  TA_RCT <- sample_n(TA_World %>% filter(RANDASSIGN==1), TA_Size)
  
  #set.seed(seed_value)
  CC_World <- sample_n(data %>% filter(RANDASSIGN==0), TA_Size)
  
  #set.seed(seed_value)
  CC_RCT <- sample_n(CC_World, CC_Size)
  
  #### External Population ####
  external_population <- setdiff(data, rbind(TA_World, CC_World)) 
  EC_World <- external_population %>% filter(RANDASSIGN==0)
  
  #W_EC_bias_ipf <- get_extraBiasSprint(dat = EC_World, maxIter = IPF_maxiter)
  W_EC_bias_ipf <- get_BiasAllhat(dat = EC_World, maxIter = IPF_maxiter)
  Biased_EC_World <- EC_World[sample(seq_len(nrow(EC_World)),
                                     size = nrow(EC_World),
                                     replace = TRUE, prob = W_EC_bias_ipf),]
  
  return(list(rbind(TA_RCT, CC_RCT), Biased_EC_World))
  
}


#monthly_SPRINT_data <- read.csv("./Data/Processed/SPRINT_monthly_bp.csv")
#datap <- data %>% mutate(EVENT = ifelse(EVENTDAYS<=1190, EVENT, 0))