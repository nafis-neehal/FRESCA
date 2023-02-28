setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
source("./Scripts/common.R")
source("./Scripts/run_ldm_scenario.R")
source("./Modules/p_val_add.R")
source("./Modules/equity_metrics_miao.R")

TA_World <- data %>% filter(RANDASSIGN==1)
TA_size <- 1000
total_seeds <- 15
total_success <- 50 #number of inner bootstraps
IPF_maxiter <- 100
step_size <- 250
randomization_ratio <- 1
num_col <- 10
CC_size_list <- seq(TA_size, 0, -step_size)
EC_size_list <- randomization_ratio*TA_size - CC_size_list

Final_LDM_Summary <- data.frame(matrix(ncol = 10, nrow = 0))
colsLDM <- c("Seed", "CC_Size", "Var", "Level", "LDM_TA", "LDM_CC", "LDM_HC_Prop", "LDM_HC_IPF", 
             "LDM_TA_Prop_IPF", "LDM_HC_Prop_IPF" )
colnames(Final_LDM_Summary) <- colsLDM

cores=detectCores(logical = TRUE)
cl <- makeCluster(cores[1]) #not to overload your computer
registerDoSNOW(cl)

packages <- c("dplyr", "ipfr", "MatchIt", "marginaleffects")

pb <- txtProgressBar(min = 0,
                     max = total_seeds * length(CC_size_list), 
                     style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

start <- Sys.time()
Final_LDM_Summary <- foreach(seed_value=seq(1,total_seeds,1), .combine = 'rbind', .inorder = TRUE)%:%
  foreach (size=CC_size_list, .combine = "rbind", .packages = packages,
           .options.snow = opts, .inorder = TRUE) %dopar%{
    run_one_ldm_scenario(seed_value, size)
  }

end <- Sys.time()
paste("Total time took:", end-start)

close(pb)
stopCluster(cl)
stopImplicitCluster()

setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
directory <- "./Data/Results/M6/"
filename <- "LDM_Summary_extrabias.csv"
write.csv(Final_LDM_Summary, paste(directory, filename, sep = ""), row.names = FALSE)



