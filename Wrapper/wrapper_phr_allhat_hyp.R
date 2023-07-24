#### SOURCES ####
setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
source("./Wrapper/load_packages.R")
source("./Wrapper/wrapper_modules.R")
source("./Wrapper/load_data_formatted.R")

allhat_data <- read.csv("./Data/ALLHAT_data/ALLHAT_processed.csv")
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


#### STUDY PARAMETERS ####
IPF_maxiter <- 100
num_seed <- 50
scale_to <- 100
weight_method <- "ipf" #"gen" or "ipf"


#### Run Main Simulation #####
run_main_simulation <- function(num_seed, pop){
  
  phr_estimate_list <- c()
  
  for (i in 1:num_seed){
    
    cat("\nIteration:", i, "\n")
    
    #Step 1: create data
    all_data <- partition_data_generic(pop, partition_seed=i)
    target_data <- data.frame(all_data[2])
    RCT <- data.frame(all_data[[1]][[1]])
    BIASED_EC <- data.frame(all_data[[1]][[2]])
    
    remove(all_data)
    cat("Data partitioned...\n")
    
    
    #Step 2: Borrow Synthetic Controls
    synthetic_controls <- get_matched_controls(method_name = "simple_matchit",
                                               RCT = RCT, BIASED_EC = BIASED_EC,
                                               N = SC_Size)
    sc_weights <- synthetic_controls$weights
    synthetic_controls$distance <- NULL
    synthetic_controls$weighs <- NULL
    synthetic_controls$RANDASSIGN <- 0
    synthetic_controls <- synthetic_controls %>% select(RANDASSIGN, Age_Group:EVENTDAYS)
    cat("Synthetic Controls borrowed...\n")

    hybrid_RCT <- rbind(RCT, synthetic_controls)
    
    #Step 3: Measure adjusted Treatment Effect
    #Use lr/rf. Lasso takes a lot of time - when using IPTW

    adjusted_weights <- get_adjusted_weights(hybrid_RCT, target_data,
                                             proba_estimate_method="lr", weight_method=weight_method)

    if(sum(unlist(adjusted_weights))!=scale_to){
      adjusted_weights <- scale_vector(unlist(adjusted_weights), scale_to)
    }

    cat("Adjusted weights generated...\n")

    #Step 4: Measure Adjusted Effects
    
    #phr_estimate <- get_adjusted_effect(hybrid_RCT, adjusted_weights)
    phr_estimate <- get_adjusted_effect(hybrid_RCT, NULL)
    
    cat("Adjusted effect estimated...\n")
    
    #append
    phr_estimate_list <- c(phr_estimate_list, exp(phr_estimate))
    
  }
    
  return (phr_estimate_list)
}

################################## Trial 1: AM vs CH ####################################

TA_Size <- nrow(tr1)
CC_Size <- 4000
SC_Size <- TA_Size - CC_Size

start <- Sys.time()
result <- run_main_simulation(num_seed=num_seed, pop=g1)
end <- Sys.time()

paste("Total time elapsed:", end-start)
meanCI <- CI(result)
paste("Mean:", meanCI["mean"])
paste("Lower:", meanCI["lower"])
paste("Upper:", meanCI["upper"])

################################## Trial 2: Li vs CH ####################################

TA_Size <- nrow(tr2)
CC_Size <- 4000
SC_Size <- TA_Size - CC_Size

start <- Sys.time()
result <- run_main_simulation(num_seed=num_seed, pop=g2)
end <- Sys.time()

paste("Total time elapsed:", end-start)
meanCI <- CI(result)
paste("Mean:", meanCI["mean"])
paste("Lower:", meanCI["lower"])
paste("Upper:", meanCI["upper"])


