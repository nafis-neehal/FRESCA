setwd("/Users/nafisneehal/ESCA")
source("./Scripts/common.R")

#directory <- "./Data/Results/PATT_m6/CC_EC_Vary_Size/summary_seed_results/matchCCEC/"
directory <- "./Data/Results/PATT_m6/CC_EC_Vary_Size/"
filename <- "1000TA_vary_CC_EC_matchCCEC_seed5.csv"
dat<- read.csv(paste(directory, filename, sep = ""))
gt <- read.csv("./Data/Results/PATT_ground.csv")

group1_data <- dat %>% filter(EC_Sample_Size==0) %>% select("Mean_Patt")
group2_data <- gt %>% select("PATT_ground")

# result <- t.test(group1_data, group2_data, var.equal = FALSE)
# result$p.value

var.test(unlist(group1_data, use.names = FALSE), unlist(group2_data, use.names = FALSE))
