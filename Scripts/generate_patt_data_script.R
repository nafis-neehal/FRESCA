setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
source("./Scripts/common.R")
source("./Scripts/run_patt_scenario.R")
source("./Modules/p_val_add.R")
source("./Modules/equity_metrics_miao.R")

TA_World <- data %>% filter(RANDASSIGN==1)
TA_size <- 2000
total_seeds <- 1 #number of total scenarios
total_success <- 5 #number of inner bootstraps
IPF_maxiter <- 100
step_size <- 500
randomization_ratio <- 1
num_col <- 7
CC_size_list <- seq(TA_size, 0, -step_size)
EC_size_list <- randomization_ratio*TA_size - CC_size_list

#final dataframe that will contain the result
Final_PATT_Summary <- data.frame(matrix(ncol = num_col, nrow = 0))
colsPATT <- c("Seed", "CC_Size", "TA_CC", "TA_HC", "TA_HC_Prop", 
              "TA_HC_IPF", "TA_HC_Both")
colnames(Final_PATT_Summary) <- colsPATT

#settings for cluster initiation and parallelization
cores=detectCores(logical = TRUE) #check number of cores
cl <- makeCluster(cores[1]) #using all but one nodes
registerDoSNOW(cl) #register the cluster

#packages to copy in each node in the cluster
packages <- c("dplyr", "ipfr", "MatchIt", "marginaleffects", "tableone", "survival")

#nested parallelization starts here
pb <- txtProgressBar(min=0, 
                     max = total_seeds * length(CC_size_list), 
                     style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

start <- Sys.time()
Final_PATT_Summary <- foreach(seed_value=seq(1,total_seeds, 1), .combine = 'rbind') %:%
  foreach (size=CC_size_list, .combine = "rbind", .packages = packages, 
           .options.snow = opts) %dopar%{
             run_one_patt_scenario(seed_value, size)
           }
end <- Sys.time()
paste("Total time took:", end-start)

#close progress bar and turn off + release cluster
close(pb)
stopCluster(cl)
stopImplicitCluster()
 
#save the data
setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
directory <- "./Data/Results/M6/"
filename <- "PATT_Summary_extrabias_2k_amia.csv"
write.csv(Final_PATT_Summary, paste(directory, filename, sep = ""), row.names = FALSE)




