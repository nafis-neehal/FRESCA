##### Note #####
# Top-N-Prop + IPTW Re-weighting ~ 8 minutes (100 runs)

#### SOURCES ####
setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
source("./Wrapper/load_packages.R")
source("./Wrapper/wrapper_modules.R")
source("./Wrapper/load_data_formatted.R")

#### STUDY PARAMETERS ####
TA_Size <- 2000
CC_Size <- 1000
SC_Size <- 1000
IPF_maxiter <- 100
num_seed <- 50
scale_to <- 100
weight_method <- "gen" #"gen" or "ipf"

run_main_simulation <- function(num_seed){

  phr_estimate_list <- c()
  
  for (i in 1:num_seed){
    
    cat("Iteration:", i, "\n")
    
    #Step 1: create data
    all_data <- read_data(partition_seed=i)
    target_data <- data.frame(all_data[2])
    RCT <- data.frame(all_data[[1]][[1]])
    BIASED_EC <- data.frame(all_data[[1]][[2]])
    
    remove(all_data)
    
    cat("Data partitioned...\n")
    
    #phr_estimate <- get_cox_effect(RCT, weights=NULL)
    
    #Step 2: Borrow Synthetic Controls
    synthetic_controls <- get_matched_controls(method_name = "one_to_one_with_replace_caliper",
                                               RCT = RCT, BIASED_EC = BIASED_EC,
                                               N = SC_Size)

    hybrid_RCT <- rbind(RCT, synthetic_controls)
    cat("Synthetic Controls borrowed...\n")

    # phr_estimate <- get_adjusted_effect(hybrid_RCT, NULL)
    
    #hybrid_RCT <- RCT
    
    #phr_estimate <- get_adjusted_effect(hybrid_RCT, NULL)
    
    #Step 3: Measure adjusted Treatment Effect
    #Use lr/rf. Lasso takes a lot of time
    
    adjusted_weights <- get_adjusted_weights(hybrid_RCT, target_data,
                                             proba_estimate_method="lr", weight_method=weight_method)

    cat("Adjusted weights generated...\n")

    if(sum(unlist(adjusted_weights))!=scale_to){
      adjusted_weights <- scale_vector(unlist(adjusted_weights), scale_to)
    }
    
    phr_estimate <- get_adjusted_effect(hybrid_RCT, adjusted_weights)
    
    cat("Adjusted effect estimated...\n")
    
    #append
    phr_estimate_list <- c(phr_estimate_list, exp(phr_estimate))
    
  }
  
  return (phr_estimate_list)
  
}

start <- Sys.time()
result <- run_main_simulation(num_seed=num_seed)
end <- Sys.time()

paste("Total time elapsed:", end-start)
meanCI <- CI(result)
paste("Mean:", meanCI["mean"])
paste("Lower:", meanCI["lower"])
paste("Upper:", meanCI["upper"])



  